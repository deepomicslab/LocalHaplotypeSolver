from Util.Utilities import *

__author__ = 'XU Chang'


def _reachability(g, startV, source=True, sink=True, targetGroup=None):
    r = False, False
    start_v = startV
    visited_stack = [startV]
    while 1:
        unvisited_next_edges = [next_e for next_e in start_v.mNextEdges if
                                not next_e.mVisited]
        if not unvisited_next_edges:
            visited_stack.remove(start_v)
            if not visited_stack:
                break
            else:
                start_v = visited_stack[-1]
                continue
        else:
            next_e = unvisited_next_edges[0]
            next_v = next_e.mTargetV
            next_seg = g.getSegment(next_v.mType, next_v.mIdx)
            if source and g.isSource(next_seg):
                if next_seg.getGroup() == targetGroup:
                    r = True, False
                    break
                else:
                    next_e.mVisited = True
                    visited_stack.append(next_v)
                    start_v = next_v
            elif sink and g.isSink(next_seg):
                if next_seg.getGroup() == targetGroup:
                    r = True, True
                    break
                else:
                    next_e.mVisited = True
                    visited_stack.append(next_v)
                    start_v = next_v
            else:
                next_e.mVisited = True
                visited_stack.append(next_v)
                start_v = next_v
        # print 'chose start v', start_v
    g.resetVisitedFlag()
    return r


def reachability_wrapper(g):
    try:
        reachability(g)
    except NotReachableException as e:
        stderr(e)
        add_normal_junctions(g, force=True)
        g.printGraphInfo()
        reachability(g)


def add_normal_junctions(g, force=True):
    if force:
        # Just do aggressive adding.
        for seg in g.mSegments:
            if g.isSink(seg):
                continue
            s_v = seg.mPosV
            n_v = chooseNextVertex(g, s_v)
            try:
                # g.addJunction(s_v.mType, s_v.mIdx, s_v.mDir,
                #               n_v.mType, n_v.mIdx, n_v.mDir,
                #               assignCov(s_v, n_v, g),
                #               cred=assignCred(s_v, n_v, g), inferred=True)
                g.addJunction(s_v.mType, s_v.mIdx, s_v.mDir,
                              n_v.mType, n_v.mIdx, n_v.mDir,
                              assignCov(s_v, n_v, g),
                              cred=SysSettings.INFERRED_JUNCTION_WEIGHT, inferred=True)
            except DuplicateJunctionException as e:
                continue
    else:
        for seg in g.mSegments:
            if g.isSink(seg):
                continue
            s_v = seg.mPosV
            n_v = chooseNextVertex(g, s_v)
            try:
                if assignCov(s_v, n_v, g) > 0.5 * 0.5 * g.mAvgCov:
                    j = g.addJunction(s_v.mType, s_v.mIdx, s_v.mDir,
                                      n_v.mType, n_v.mIdx, n_v.mDir,
                                      assignCov(s_v, n_v, g),
                                      cred=SysSettings.INFERRED_JUNCTION_WEIGHT, inferred=True)
            except DuplicateJunctionException as e:
                continue


def virusOriginallyReachesGroup(g, vtx):
    """Reach all groups that a virus segment can reach,
    So the input seg must be a virus segment
    """
    # reachable_group = []
    # for e in vtx.mNextEdges:
    #     if e.mTargetV.mType == 'H' and e.mTargetV.getGroup() not in reachable_group:
    #         reachable_group.append(e.mTargetV.getGroup())
    start_v = vtx
    visited_stack = [start_v]
    reachable_group = []
    while 1:
        unvisited_next_edges = []
        for next_e in start_v.mNextEdges:
            if not next_e.mTargetV.mType == 'H':
                if not next_e.mVisited:
                    unvisited_next_edges.append(next_e)
            else:
                if next_e.mTargetV.getGroup() not in reachable_group:
                    reachable_group.append(next_e.mTargetV.getGroup())

        if not unvisited_next_edges:
            visited_stack.remove(start_v)
            if not visited_stack:
                break
            else:
                start_v = visited_stack[-1]
                continue
        else:
            next_e = unvisited_next_edges[0]
            next_v = next_e.mTargetV
            next_e.mVisited = True
            visited_stack.append(next_v)
            start_v = next_v
    g.resetVisitedFlag()
    return reachable_group


def reachability(g):
    """If there are dead end segments, add necessary junctions.

    While it is required to input all segments, user may not list all
    junctions. Hence the need to ensure all segments are reachable to the
    source AND sink segments, i.e. no segment is a dead-end.

    > Orphan:
    Orphan segments(segments without in AND out edges) are not considered
    unless they are added into the graph while trying to ensure reachability
    of other non-orphan segments.

    > Choosing the next segment:
    We assume that the user has provided all abnormal junctions (
    junctions different from reference) and try to add edges according to
    the reference sequence. However, chances are, some important junctions
    are not provided and as a result, no next vertex can be found by our
    conservative method. User can set "exclude_incompatible_edges" as TRUE
    to let the program decide which edge to discard. (see below)
    **TODO**: In the future, we may recommend some possible junctions.

    > DuplicateJunctionException
    The exception occurs when we cannot safely infer a junction to the next
    segment. Since host segments is guaranteed to be able to connect to
    source or sink given the way our code is written, this exception only
    happens to virus segments. Under this situation, we try to remove
    junctions related to the segment. By allowing the junctions to have 0
    copy, it is up to Linear Programming to decide which part of the graph to
    dump.

    > Why choose THAT segment as the possible next segment?
    In the TestCases folder, file "reachability" best describes the issue.
    Why not connect V1+ => V2+ only and discard H2+ => H3+ => H4+ ? How do
    we keep artifacts like that under control? Check out details at function:
    assignCred(). Notice that if segment is originally orphan, its weight
    lower bound can be set as 0 during LP. (see BalanceGraph.py for detail)


    Added edge will have a credibility of the average of the two segments
    it bridges. And at least it will be assigned with half of the genome-wide
    average coverage, which stands for the expected coverage of a single copy.

    """
    v_needs_reachability = None

    for seg in g.mSegments:
        addedJunctions = []
        if g.isSource(seg) or g.isSink(seg):
            continue
        s = seg.mPosV
        if not s.mNextEdges and not s.mPrevEdges:  # in case of orphans
            if s.mType == 'H' and s.getCov() >= 0.25 * g.mAvgCov:
                # 20170914 update, some heuristic rules to keep a human segment in the graph.
                pass
            else:
                continue

        # update 20170913, check the groups a virus seg connects to.
        if seg.mType == 'V':
            l1 = virusOriginallyReachesGroup(g, seg.mPosV)
            l2 = virusOriginallyReachesGroup(g, seg.mNegV)
            groupsToReach = list(set(l1 + l2))
        else:
            groupsToReach = [seg.getGroup()]

        pos_to_sink = None  # positive vertex is reachable to sink
        for gp in groupsToReach:
            set_to_source, set_to_sink = True, True
            while 1:
                _do_not_break = False
                while True:
                    pos_r, pos_to_sink = _reachability(g, seg.mPosV,
                                                       source=set_to_source, sink=set_to_sink,
                                                       targetGroup=gp)
                    if pos_r:
                        break  # reachable to either source or sink.
                    else:
                        n = chooseNextVertex(g, s)
                        try:
                            # j = g.addJunction(s.mType, s.mIdx, s.mDir,
                            #                   n.mType, n.mIdx, n.mDir,
                            #                   assignCov(s, n, g),
                            #                   cred=assignCred(s, n, g), inferred=True)
                            j = g.addJunction(s.mType, s.mIdx, s.mDir,
                                              n.mType, n.mIdx, n.mDir,
                                              assignCov(s, n, g),
                                              cred=SysSettings.INFERRED_JUNCTION_WEIGHT, inferred=True)
                            addedJunctions.append(j)
                            '''
                            20170919 update, inferred edge should have no lower bound,
                            Below code shouldn't be here but it is, as a reminder.
                            In cases of multiple integration sites (HBV_262T_chr4_chr13), we 
                            leverage this to remove the inferred edge that should not exist.
                            (edges such as V4n->V3n are added since V4+ connects to chr4 but V4-
                            does not, while actually V4 should only connnect to chr13.)
                            '''
                            j.mLowerBoundLimit = False
                        except DuplicateJunctionException as e:
                            # > For virus not connecting back to host or user declared
                            #    a normal edge.
                            # > save the start vertex and throw the problem to the next
                            #    edge, if the problem returns back to the start vertex,
                            #    then there is a big problem
                            v_needs_reachability = seg.mPosV
                            if v_needs_reachability != n:
                                s = n
                                continue
                            else:
                                reachabilityErrorHandler(e, seg)
                                lowerSegAndJunCred(seg, g)
                                break
                        s = n
                s = seg.mNegV
                while True:
                    if pos_to_sink:  # check if its neg seg reaches source too
                        neg_r, neg_to_sink = _reachability(g, seg.mNegV, sink=False, targetGroup=gp)
                    else:
                        neg_r, neg_to_sink = _reachability(g, seg.mNegV, source=False, targetGroup=gp)
                    if neg_r and pos_to_sink != neg_to_sink:
                        break
                    else:  # pos and neg segs should not be reachable to the same end.
                        """There's an issue here, if pos reaches sink-pos, then neg
                        should reach source-neg, if pos reaches sink-neg, then neg
                        should reach source-pos, not source-neg. I did not perform
                        the check here, since LP will do that for me. However,
                        there are times when it goes wrong and triggers
                        NotReachableException, It's all because of the input is
                        not sufficient to construct a graph. Not trying to fix the
                        problem here."""
                        try:
                            n = chooseNextVertex(g, s)
                            try:
                                # j = g.addJunction(s.mType, s.mIdx, s.mDir,
                                #                   n.mType, n.mIdx, n.mDir,
                                #                   assignCov(s, n, g),
                                #                   cred=assignCred(s, n, g), inferred=True)
                                j = g.addJunction(s.mType, s.mIdx, s.mDir,
                                                  n.mType, n.mIdx, n.mDir,
                                                  assignCov(s, n, g),
                                                  cred=SysSettings.INFERRED_JUNCTION_WEIGHT, inferred=True)
                                addedJunctions.append(j)
                                j.mLowerBoundLimit = False
                                # print 'added a junc', j, j.mCov.mCov
                            except DuplicateJunctionException as e:
                                v_needs_reachability = seg.mNegV
                                if v_needs_reachability != n:
                                    s = n
                                    continue
                                else:
                                    reachabilityErrorHandler(e, seg)
                                    lowerSegAndJunCred(seg, g)
                                    break
                        except NotReachableException:
                            '''
                            Try to solve the case where H2+ connects to source, so that 
                            H2- needs to connect to sink, which is wrong.
                            '''
                            s = seg.mPosV
                            g.delJunctions(*addedJunctions)
                            if not pos_to_sink:
                                # This is usually the case, pos connects to source.
                                set_to_source, set_to_sink = False, True
                                _do_not_break = True
                                break
                        _do_not_break = False
                    s = n
                    set_to_source, set_to_sink = True, True
                if not _do_not_break:
                    break  # else then continue the outer while loop
    # return True


def isNeighbor(g, vertex1, vertex2):
    """When a normal connection is declared in input but we try to ensure
    reachability in the graph. Duplicate junction exception may be raised,
    if the junction to add fits below standard, then there won't be error.

    vertex1 is the source of a junction and vertex2 is the target

    Only when two vertices are from host segment,
    """
    if vertex1.mType == vertex2.mType:
        if vertex1.mType == 'V':
            return False
        if vertex1.mDir == vertex2.mDir == '+':
            if vertex1.mIdx + 1 == vertex2.mIdx:
                return True
            # if vertex1.mType == 'V' and SysSettings.CIRCULAR_VIRUS:
            #     if vertex1.mIdx == g.mMaxVIdx and vertex2.mIdx == g.mMinVIdx:
            #         return True
        if vertex1.mDir == vertex2.mDir == '-':
            if vertex1.mIdx - 1 == vertex2.mIdx:
                return True
            # if vertex1.mType == 'V' and SysSettings.CIRCULAR_VIRUS:
            #     if vertex2.mIdx == g.mMaxVIdx and vertex1.mIdx == g.mMinVIdx:
            #         return True
    return False


def weightedCred(vertex, avg_cov, is_start=True):
    """Return weighted credibility for a vertex during junction inference.

    average of (own_weight * overall_weight * (own_cov - in/out_cov) / avg_cov)

    surely, the difference between segment cov and in/out junction cov can be
    big, caused by high amount of reads supporting normal segment, making the
    credibility returned a high number.
    """
    if vertex.mType == 'H':
        coef = SysSettings.HOST_SEG_LP_COE * vertex.mCredibility
    else:
        coef = SysSettings.VIRUS_SEG_LP_COE * vertex.mCredibility
    if is_start:
        cred = coef * max(0, vertex.getCov() - vertex.getOutCov()) / (avg_cov / 2)
        return cred
    else:
        cred = coef * max(0, vertex.getCov() - vertex.getInCov()) / (avg_cov / 2)
        return cred


def assignCred(s_v, t_v, g):
    """Assign credibility for inferred junctions.

    > Is it affected by the weight of their surrounding segments?
    1. If a junction is user-defined, then no.
    2. If the junction is inferred, then we need to consider both segments'
    weight and their credibility settings. Consider "reachability" test case.

    the value is capped to 1.
    """
    c = (weightedCred(s_v, g.mAvgCov, is_start=True) +
         weightedCred(t_v, g.mAvgCov, is_start=False)) / 2
    return min(1, c)
    # return (s_v.mCredibility + t_v.mCredibility) / 2.0


def assignCov(s_v, t_v, g):
    in_cov_to_fill = max(0, t_v.getCov() - t_v.getInCov())
    out_cov_to_fill = max(0, s_v.getCov() - s_v.getOutCov())
    # return max(0.01, (in_cov_to_fill + out_cov_to_fill) / 2.0) * (a + b) / 2
    # return max(0, (in_cov_to_fill + out_cov_to_fill) * SysSettings.INFER_COV_COE / 2.0)
    return min(in_cov_to_fill, out_cov_to_fill)


def chooseNextVertex(g, startV):
    """Returns the vertex that a normal edge starting from startV points to.

    Human segs are not circular.
    """
    nextV = None
    try:
        if startV.mType == 'H':
            if startV.mDir == '+':
                nextV = g.getSegment('H', startV.mIdx + 1).mPosV
            else:
                nextV = g.getSegment('H', startV.mIdx - 1).mNegV
            if not g.inSameGroup(startV, nextV):
                """
                Requires that human segs in same group are in consecutive order.
                Since source and sink are not considered in this function, therefore
                the choice of nextV will be within the group.
                """
                raise NotInSameReferenceGroupException('from {} to {}'.format(startV, nextV))
        if startV.mType == 'V':
            if startV.mDir == '+':
                if startV.mIdx == g.mMaxVIdx:
                    if SysSettings.CIRCULAR_VIRUS:
                        nextV = g.getSegment('V', g.mMinVIdx).mPosV
                else:
                    nextV = g.getSegment('V', startV.mIdx + 1).mPosV
            else:
                if startV.mIdx == g.mMinVIdx:
                    if SysSettings.CIRCULAR_VIRUS:
                        nextV = g.getSegment('V', g.mMaxVIdx).mNegV
                else:
                    nextV = g.getSegment('V', startV.mIdx - 1).mNegV
    except AttributeError as ae:
        import traceback
        raise NotReachableException('{}\nCannot choose next V for {}.\n'
                                    'Make sure the numbering of segments from '
                                    'the same group are consecutive.'
                                    .format(traceback.format_exc(), startV))
    if not nextV:  # I guess this error should not occur for circular virus.
        sys.exit('No next vertex are found suitable while trying to ensure'
                 'reachability of segments.\nAt vertex {}'.format(startV))
    return nextV


reach_error_msg = """
> Error while trying to add necessary edges for segment [ {} ].
This implies that the segment is not compatible in the graph with
supplied breakpoint information. For example, there should, at least,
exist an edge that starts from host and ends at virus, and an edge
starts from virus and ends at host. The program will not assume any
abnormal junction that bridges host and virus. If you have turned on
"exclude_incompatible_edges", The program will try to remove all edges
with the segment involved which will, in effect, orphan the segment.

"""
reach_error_msg_short = "> Error while trying to add necessary edges for segment [ {} ].\n"
ReachabilityErrorRaised = False


def reachabilityErrorHandler(exception, segment):
    global reach_error_msg, ReachabilityErrorRaised
    if not ReachabilityErrorRaised:
        error_msg = str(exception) + reach_error_msg.format(segment.getAbbr())
        ReachabilityErrorRaised = True
    else:
        error_msg = str(exception) + reach_error_msg_short.format(segment.getAbbr())
    if SysSettings.DELETE_INCOMPATIBLE_EDGES:
        sys.stderr.write(error_msg)
    else:
        sys.exit(error_msg)


def lowerSegAndJunCred(seg, g):
    seg.setCred(0)
    seg.mLowerBoundLimit = False
    for e in seg.mPosV.mNextEdges:
        j = g.mJunctions[g.getJunctionIdx(e)]
        j.setCred(0)
        j.mLowerBoundLimit = False

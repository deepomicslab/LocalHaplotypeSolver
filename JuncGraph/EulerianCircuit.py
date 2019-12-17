import SysSettings
import random
from Util.Utilities import *
from JuncGraph.Graph import Graph
import time


__author__ = 'XU Chang'

FORKING_CHOICE = {}

def findCycle_recur(g, startV):
    """Deprecated, recursion is no good"""
    cycleFound = [startV]
    def _findCycle(aFrom):
        next_edge = selectNextEdge(aFrom, g, startV)
        aFrom.mLastEdgeChoice = next_edge
        next_edge.traverse()
        cycleFound.append(next_edge.mTargetV)
        if next_edge.mTargetV == startV:
            return
        else:
            _findCycle(next_edge.mTargetV)

    _findCycle(startV)
    return cycleFound


def findCycle(g, startV):
    """Find a cycle that begins from startV and ends at startV

    Returns the cycle, and its count.

    Important: only traversal with memory can do the trick marked by (1).
      If there is no cycle within the startV -> startV cycle, then by
      the definition of Eulerian circuit, all edges in the cycle can
      be further traversed with cn=min(traversed_edges_remaining_weight)

      If there are cycles within the startV -> startV cycle, then by
      "traversal with memory" method, all inner cycles are drained,
      causing the min() function returning 0, no problems there.

      However, if we are not traversing with memory, we will have below
      situation:

      startV--e1-->CYC1-->CYC1--e2-->startV,

      after initial traversal, we may have
      min_weight_left(e1)=1, min_weight_left(e2)=1,
      min_weight_left(edges_in_CYC1)=1

      However, the simplified deduction will subtract 2 from each edge
      in CYC1, because there are 2 copies of CYC1 during the initial
      traversal, which results in their weight < 0, which is wrong.


    """
    cycleFound = [startV]

    aFrom = startV
    remaining_cn = 0
    _min_remain_w = float('inf')
    traversed_edges = []
    # cnt = 0
    while 1:
        # cnt += 1
        next_edge = selectNextEdge(aFrom, g, startV)
        # print
        next_v = next_edge.mTargetV
        aFrom.mLastEdgeChoice = next_edge
        traversed_edges.append(next_edge)
        _min_remain_w = min(_min_remain_w, next_edge.traverse())
        # if cnt > 100:
        #     time.sleep(0.1)
        #     print next_v
        if _min_remain_w < 0:
            # Bug fix, if this is smaller than 0, traverse could never end.
            # Treat this as invalid cycle and ignore! This is caused by
            # breakpoint choice that ends with a vertex being traversed
            # both positively and negatively in a cycle.
            raise TooManyTraversalException()
        cycleFound.append(next_v)
        if next_v == startV:
            if (SysSettings.EDGE_SELECTION == 'try_all' or
                    SysSettings.EDGE_SELECTION == 'smart_random'):
                # (1)
                remaining_cn = _min_remain_w
                if remaining_cn:
                    [e.traverse(cn=remaining_cn) for e in traversed_edges]
            break
        else:
            aFrom = next_v

        if next_edge.isImaginaryJunction():
            """
            e.g. during [H2-] => [H1-] => [H(sink)-], the above "break"
             will not work
            """
            raise ImaginaryJunctionReachedException(
                'Cannot find cycle. Please try disabling "preferred" '
                'junctions, or try disable --lp_with_original_ploidy')
    '''
    20170918 update, to solve situations like < V5+ V1+ H4- H6- V4+ V5+ > 
     where H4 is actually the source.
    '''
    # print 'cycle', [v.getAbbr() for v in cycleFound]
    # for i, v in enumerate(cycleFound[1:-1]):
    #     if g.isSource(v.mSeg):
    #         # change < H2+ H3+ H1+ H2+ > to < H1+ H2+ H3+ H1+ >
    #         cycleFound = cycleFound[i:] + cycleFound[1:i] + [v]
    #         if v.mDir == '-':
    #             tmp = []
    #             for _v in reversed(cycleFound):
    #                 tmp.append(_v.getSiblingVtx())
    #             cycleFound = tmp
    #         break
    # print [v.getAbbr() for v in cycleFound]
    # print 'after', startV, startV.getWeight()
    return cycleFound, int(1 + remaining_cn)


def selectNextEdge(aFrom, g, cycleStartV):
    if SysSettings.EDGE_SELECTION == 'default':
        nextEdges = aFrom.mNextEdges
        group = aFrom.getGroup()
        if cycleStartV == g.mSource[group].mPosV:
            for e in nextEdges:
                if (e.hasWeight() and e.mTargetV.mType == aFrom.mType and
                            e.mTargetV.mIdx == aFrom.mIdx + 1):
                    return e
            for e in nextEdges:
                if (e.hasWeight() and e.mTargetV.mType == aFrom.mType and
                            e.mTargetV.mIdx > aFrom.mIdx):
                    """This solves the issue described in detail below.
                    Always choose the larger index."""
                    return e
        for e in nextEdges:
            if e.hasWeight() and e.mTargetV == cycleStartV:
                return e
        for e in nextEdges:
            if e.hasWeight():
                return e
    elif SysSettings.EDGE_SELECTION == 'random':
        nextEdges = [e for e in aFrom.mNextEdges if e.hasWeight()]
        index = random.randint(0, len(nextEdges) - 1)
        return nextEdges[index]
    elif SysSettings.EDGE_SELECTION == 'smart_random':
        previous_choice = aFrom.mLastEdgeChoice
        if previous_choice and previous_choice.hasWeight():
            return previous_choice
        else:
            nextEdges = [e for e in aFrom.mNextEdges if e.hasWeight()]
            for e in nextEdges:
                if e.mPreferred:
                    return e
            """If the cycle being traversed is in minus strand of human
            reference, we try to bring it back to the plus strand. if not
            successful, the path will reach H(source)- and continues to
            H(sink)-, which is wrong since H(source) <=> H(sink) junction is
            an imaginary junction. Notice that this error DOES NOT imply that
            the input file does not provide enough information in constructing
            a proper graph (since reachability and degree balancing has
            been ensured in previous operations). Therefore, with below
            code, such error will not occur."""
            if aFrom.mType == 'H' and aFrom.mDir == '-':
                # filtered_nextEdges = [e for e in nextEdges if
                #                       not e.connectsSameDir() or
                #                       (e.connectsSameDir() and
                #                        e.mTargetV.mIdx > aFrom.mIdx)]

                # Try to return back to plus strand.
                filtered_nextEdges = [e for e in nextEdges if not e.connectsSameDir()]
                if filtered_nextEdges:
                    nextEdges = filtered_nextEdges
            # if len(nextEdges) == 0:
            #     raise TraversalException(
            #         'Vertex {} has no available next edges left. '
            #         'Please try disabling "preferred" junctions, '
            #         'or try disable --lp_with_original_ploidy'.format(aFrom))
            index = random.randint(0, len(nextEdges) - 1)
            return nextEdges[index]
    elif SysSettings.EDGE_SELECTION == 'try_all':
        global FORKING_CHOICE
        # print [a.getAbbr() for a in aFrom.mNextEdges]
        if len(aFrom.mNextEdges) == 1:
            return aFrom.mNextEdges[0]
        elif len(aFrom.mSafeNextEdges) == 1:
            return aFrom.mSafeNextEdges[0]
        else:
            if not aFrom.mSafeNextEdges:
                # Then it might have only one effective edge, check it here.
                forks = []
                for e in aFrom.mNextEdges:
                    if not g.edgeReachesImaginaryJunction(e):
                        forks.append(e)
                        # Update this so that edgeReachesImaginaryJunction() is less called.
                        aFrom.mSafeNextEdges.append(forks[0])
                if len(forks) == 1:
                    # Still do not need to be a fork point
                    return forks[0]
            # More than one edge without imaginary junction.
            try:
                for next_e_idx in FORKING_CHOICE[aFrom]:
                    e = aFrom.mNextEdges[next_e_idx]
                    if e.hasWeight():
                        # if g.edgeReachesImaginaryJunction(e):
                        #     raise ImaginaryJunctionReachedException()
                        return e
                # print 'i was also in ', aFrom
            except KeyError:  # happens at FORKING_CHOICE[aFrom]
                # This should never happen if all forks are used for permutation and combination.
                # stderr(ex)
                # return [e for e in aFrom.mNextEdges if e.hasWeight()][0]
                #raise ImaginaryJunctionReachedException('Key error at FORKING_CHOICE[aFrom] for {}'.format(
                # aFrom.getAbbr()))
                x = BreakPointNotEnumeratedException()
                x.breakPoint = aFrom
                raise x
            """
            The graph is balanced, so there is no way below code can be reached 
             under normal circumstances. However, We excluded edge that reaches
             imaginary junction from the list of fork points. Therefore this 
             part can be reached when the path is going into the imaginary 
             junction (e.g. when aFrom=H2- and H2- connects to H1-)
            """

            raise ImaginaryJunctionReachedException()
    else:
        sys.exit('Please select a valid mode for choosing next edges at '
                 'traversal.')


def findLGM(g, fork_pts=None):
    # type: (Graph, list) -> (list, list)
    # Returns the path that uses up all weights.
    # LGM=[[v1,v2,v3,...] [v2,v4,...]]

    # s = start_v if start_v else g.getStartVertex() # get source positive vertex
    ALL_LGMs = []
    ALL_LGMs_cnt = []

    # First we drain source and sink from all groups.
    while 1:
        st_vtx = g.getStartVertex()
        if not st_vtx:
            break
        cycle, cycle_cnt = findCycle(g, st_vtx)
        ALL_LGMs.append(cycle)
        ALL_LGMs_cnt.append(cycle_cnt)
    # Then we check if any vertex in the lgm is not drained.
    br = False
    while not br:
        _lgm_flat = [v for c in ALL_LGMs for v in c]  # [[a,b], [c,d]] => [a,b,c,d]
        for idx, v in enumerate(_lgm_flat):
            if v.hasWeight():
                st_vtx = v
                cycle, cycle_cnt = findCycle(g, st_vtx)
                ALL_LGMs.append(cycle)
                ALL_LGMs_cnt.append(cycle_cnt)
                break
            elif idx == len(_lgm_flat) - 1:
                # Reached the last vertex of the lgm
                br = True

    return ALL_LGMs, ALL_LGMs_cnt
    # s = g.getStartVertex()
    # print 'got source', s
    # if not s: # Until all sources are drained.
    #     break
    # br = False
    # LGM = []
    # while not br:
    #     cycle = findCycle(g, s, group=s.getGroup())
    #     print 'raw cycle', [v.getAbbr() for v in cycle]
    #     LGM.append(cycle)
    #     _lgm_flat = [v for c in LGM for v in c]  # [[a,b], [c,d]] => [a,b,c,d]
    #     for idx, v in enumerate(_lgm_flat):
    #         if v.hasWeight():
    #             s = v
    #             break
    #         elif idx == len(_lgm_flat) - 1:
    #             br = True
    #
    # printLGM(ALL_LGMs)
    # ALL_LGMs += LGM




def mergeCycles(cycles):
    """Merge cycles into an LGM, preferably with some rule.
        (not implementing it for now).

        :param cycles: [ [cyc1], [cyc2] ], expect cyc1[0] = source

        We might not be able to merge the cycles into one composite cycle.
        So this function tries to merge all cycles. If there are cycles
        left alone, return them.

        # If the lgm is an legal LGM, then if we take the cycle starting from
        # source (naturally it means the cycle ends with the sink vertex), we
        # must be able to merge the cycles into one LGM. If the lgm param
        # provided cannot merge into one LGM, then the input is somehow
        # problematic. Maybe the sequencing coverages does not reflect the
        # ploidy value provided.
    """
    merged_lgm = cycles[0]
    lgm_dup = [c for c in cycles[1:]]  # a replica of lgm
    while True:
        _new_merge = False
        _lgm_del = []
        for c in lgm_dup:  # by default, we merge cycles[1:] into cycles[0]
            _ccycle = _mergeCycles(merged_lgm, c)  # composite cycle
            if _ccycle:  # success merge
                _new_merge = True
                merged_lgm = _ccycle
                _lgm_del.append(c)
        for c in _lgm_del:
            lgm_dup.remove(c)  # delete cycle if already been merged into cycle[0]
        if not _new_merge:
            break

    return merged_lgm, lgm_dup
    # if len(lgm_dup) > 1:
    #     return
    # if len(cycles) > 1 and merged_lgm == cycles[0]:
    #     return merged_lgm, lgm_dup  # TODO: delete first
    # return merged_lgm


def _mergeCycles(cyc1, cyc2):
    merged_cycle = [v for v in cyc1]  # take cyc1 as template, avoid list change.
    if cyc2[0] in merged_cycle:
        i = merged_cycle.index(cyc2[0]) + 1
        merged_cycle[i:i] = cyc2[1:]
        return merged_cycle
    return


def getBreakPoints(cycle):
    breakPoints = []
    for v in cycle:
        if len(v.mNextEdges) > 1:
            breakPoints.append(v)
    return breakPoints


def breakdownLGM(merged_lgm, break_point_v_list):
    """try to break down big cycles into unit cycles

        :param merged_lgm:  one composite cycle in list form.
                            Not necessarily the complete lgm. [v1, v2, ..., v1]

        :param break_point_v_list:  a list of breakpoints. [v1, v2, ...]
                                    (vertices with > 1 junctions going out.)
    """
    cycles = [[v for v in merged_lgm]]
    for v in break_point_v_list:
        old_cycles_len = len(cycles)  # remember the length of cycles. 1 for the first for-loop.
        while 1:
            cycles = _breakdown_q(cycles, v)  # cycles is updated in every while-loop.

            if len(cycles) == old_cycles_len:  # if length not changed, go check next breakpoint.
                break
            old_cycles_len = len(cycles)
    return cycles


def _breakdown(cycle_list, break_point_v):
    """
        :param cycle_list: composite cycle to be broken down. [[v1, v2, ..., v1]]
        :param break_point_v: break cycle_list with this vertex.

    Using an offset to lower the time cost.

    When appending new element, while loop in parent function continues."""
    for cyc_idx in range(len(cycle_list)):
        c = cycle_list[cyc_idx]
        for idx1 in range(len(c) - 1):
            # find idx1 from left to right, find idx2 from right to left.
            if c[idx1] == break_point_v:
                for idx2 in reversed(range(len(c) - 1)):  # len(c) - 1 since we ignore first or last element.
                    v2 = c[idx2]
                    if idx2 == idx1:
                        break
                    if v2 == break_point_v:
                        # print idx1, idx2
                        new_cycle = c[idx1: idx2]
                        new_cycle.append(v2)
                        del c[idx1: idx2]
                        cycle_list.append(new_cycle)
                        return cycle_list
    return cycle_list


def _breakdown_q(cycle_list, break_point_v):
    new_cycle_list = []
    for c in cycle_list:
        v_occur = [_i for _i in xrange(len(c) - 1) if c[_i] == break_point_v]
        if len(v_occur) >= 2:
            if v_occur[0] != 0:
                # source AND sink are not break_point_v
                new_cycle_list.append(c[:v_occur[0]] + c[v_occur[-1]:])

            new_cycle_list.extend([c[v_occur[_]: v_occur[_ + 1] + 1] for _ in xrange(len(v_occur) - 1)])
            # Same effect as below
            # for _i, o in enumerate(v_occur[:-1]):
            #     new_cycle = c[o: v_occur[_i+1] + 1]
            #     new_cycle_list.append(new_cycle)
        else:
            new_cycle_list.append(c)
    return new_cycle_list


def _breakdown_qq(cycle_list, break_point_v):
    new_cycle_list = []
    for cyc in cycle_list:
        reached_first_bp = False
        new_cyc = []
        base_cyc = []
        for v in cyc:
            if v == break_point_v:
                if reached_first_bp:
                    # End of a new cycle
                    new_cyc.append(v)
                    new_cycle_list.append(new_cyc)
                    new_cyc = [v]
                else:
                    # In base cycle
                    reached_first_bp = True
                    new_cyc = [v]
            else:
                if reached_first_bp:
                    new_cyc.append(v)
                else:
                    base_cyc.append(v)
        if v != break_point_v:
            # need to fill base cyc
            new_cycle_list.append(base_cyc + new_cyc)

    return new_cycle_list

def findWeightedCycle(g, cn_cycle_dict, target_cn, start_vtx, last_try=False):
    """find a composite cycle with weight specified

        :param g: graph instance
        :param cn_cycle_dict: {copy_number: [ [cyc1], [cyc2] ]}
        :param target_cn: cut-off number for selecting cycles

        Well we do assume that if run with target_cn = 1, all cycles cna be
        merged into one lgm.
    """
    cycles_to_merge = []
    br = False
    # Find a cycle as template (the cycle with start vertex)
    for cn in cn_cycle_dict:
        if cn >= target_cn:  # major or minor allele frequency
            if SysSettings.EDGE_SELECTION == 'try_all':
                for cyc in cn_cycle_dict[cn]:
                    if not g.isNormalAllele(cyc[:-1]) and cyc[0] == start_vtx:
                        cycles_to_merge.append(cyc)
                        br = True
                        break
            if not br:
                for cyc in reversed(cn_cycle_dict[cn]):
                    # backward, leave normal map to lower cn allele, only for normal way of traversing,
                    # does not work for brute-forcing.
                    # if cyc[0] == g.getStartVertex():
                    if cyc[0] == start_vtx:
                        cycles_to_merge.append(cyc)  # cycles_to_merge contains only one cycle with start vertex.
                        br = True
                        break
        if br:
            break
    merged = []
    # Two possibilities: 1. cycles_to_merge = [], 2. this function cannot use up all cycles according to target_cn.
    if not cycles_to_merge:
        return merged  # []
    cycles_to_lower_cn = [cycles_to_merge[0]]
    br = False
    while True:  # drain the cn_cycle_dict!
        if target_cn < 1:
            break
        # Find all cycles with copy number >= cut-off allele frequency.
        for cn in cn_cycle_dict:
            if cn >= target_cn:
                for cyc in cn_cycle_dict[cn]:
                    # if cyc[0] != g.getStartVertex():
                    if cyc[0] != start_vtx:
                        cycles_to_merge.append(cyc)
                        cycles_to_lower_cn.append(cyc)  # add all, delete not merged later.
        if len(cycles_to_merge) > 1:
            _merged_tmp, not_merged = mergeCycles(cycles_to_merge)
            merged = _merged_tmp
            if _merged_tmp == cycles_to_merge[0]:
                br = True  # if merge is not successful. break, try next template.
            for cyc in not_merged:
                cycles_to_lower_cn.remove(cyc)
        else:
            merged = cycles_to_merge[0]
            br = True
        # Lower copy number according to cycles_to_lower_cn
        for cyc in cycles_to_lower_cn:
            for cn in cn_cycle_dict:
                if cyc in cn_cycle_dict[cn]:
                    cn_cycle_dict[cn].remove(cyc)
                    result_cn = cn - target_cn
                    if result_cn > 0:
                        if result_cn in cn_cycle_dict:
                            cn_cycle_dict[result_cn].append(cyc)
                        else:
                            cn_cycle_dict.update({result_cn: [cyc]})
                    break
        if br:
            break
        cycles_to_merge = [merged]

        cycles_to_lower_cn = []
    return merged

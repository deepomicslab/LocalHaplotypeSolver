from pulp import *
import SysSettings
from Util.Utilities import *

__author__ = 'XU Chang'


def solveLp(g, ploidy_dict, low_bound=1, rm_inferred=False):
    """Ensure balance of the graph with a Linear Programming approach.

     Notice that junctions with below properties can be discarded (lowbound=0)
     1. Inferred
     2. Situation in the "deadend" test case.

     Segments with below properties can be discarded:
     1. Originally an orphan, i.e. all junctions connecting to it are inferred.
     2. "deadend" scenario.
     """

    # Variable lists
    var_seg_tar = []  # target weight of segments
    var_jun_tar = []  # target weight of junctions [LpVar_jun1, LpVar_jun2,...]
    var_seg_eps = []  # difference between original and target weight
    var_jun_eps = []
    var_seg_insum = []  # Variable, segments' in-weight
    var_seg_outsum = []  # Variable, segments' out-weight

    v_seg_cnt = 0
    h_seg_cnt = 0
    # Expression lists
    seg_expr_v_tmp = []
    seg_expr_h_tmp = []
    jun_expr_tmp = []

    for idx, j in enumerate(g.mJunctions):
        # if (rm_inferred and j.mInferred) or not j.mLowerBoundLimit:
        if j.mInferred or not j.mLowerBoundLimit:
            _jun_tar = LpVariable('{0}_j'.format(j.getVariableName()),
                                  lowBound=0, cat='Integer')
            # print 'gave ', j, 'lowerbound0', j.mInferred, j.mLowerBoundLimit
        else:
            _jun_tar = LpVariable('{0}_j'.format(j.getVariableName()),
                                  lowBound=low_bound, cat='Integer')
            # print 'gave ', j, 'lowerbound1', j.mInferred, j.mLowerBoundLimit, low_bound
        _jun_eps = LpVariable('{0}_e'.format(j.getVariableName()),
                              lowBound=0, cat='Continuous')
        var_jun_tar.append(_jun_tar)
        var_jun_eps.append(_jun_eps)
        jun_expr_tmp.append((_jun_eps, j.mCredibility))

    for s in g.mSegments:
        if s.isDeadEnd():  # if s is an orphan.
            continue
        if s.mIsOrphan or not s.mLowerBoundLimit:  # orphan as user input
        # if not s.mLowerBoundLimit:  # orphan as user input
            stderr('Segment [ {} ] is orphan in user input.'.format(s.getAbbr()))
            _seg_tar = LpVariable('{0}{1}'.format(s.mType, s.mIdx),
                                  lowBound=0, cat='Integer')
        else:
            _seg_tar = LpVariable('{0}{1}'.format(s.mType, s.mIdx),
                                  lowBound=low_bound, cat='Integer')
        _seg_eps = LpVariable('{0}{1}_e'.format(s.mType, s.mIdx),
                              lowBound=0, cat='Continuous')
        var_seg_tar.append(_seg_tar)
        var_seg_eps.append(_seg_eps)

        _var_out_edges = []
        outsum_expr_tmp = []
        for e in s.mPosV.mNextEdges:
            _tar = var_jun_tar[g.getJunctionIdx(e)]
            _var_out_edges.append(_tar)
            _coe = 1 if not e.isPalindrome() else 2
            outsum_expr_tmp.append((_tar, _coe))
        var_seg_outsum.append(LpAffineExpression(outsum_expr_tmp))

        _var_in_edges = []
        insum_expr_tmp = []
        for e in s.mPosV.mPrevEdges:
            _tar = var_jun_tar[g.getJunctionIdx(e)]
            _var_in_edges.append(_tar)
            _coe = 1 if not e.isPalindrome() else 2
            insum_expr_tmp.append((_tar, _coe))
        var_seg_insum.append(LpAffineExpression(insum_expr_tmp))

        # old code, not considering palindrome.
        # _var_out_edges = [var_jun_tar[g.getJunctionIdx(e)] for e in
        #                   s.mPosV.mNextEdges]
        # _var_in_edges = [var_jun_tar[g.getJunctionIdx(e)] for e in
        #                  s.mPosV.mPrevEdges]
        # var_seg_outsum.append(
        #     LpAffineExpression([(x, 1) for x in _var_out_edges]))
        # var_seg_insum.append(
        #     LpAffineExpression([(x, 1) for x in _var_in_edges]))
        if s.mType == 'V':
            v_seg_cnt += 1
            seg_expr_v_tmp.append((_seg_eps, s.mCredibility))
        else:
            h_seg_cnt += 1
            seg_expr_h_tmp.append((_seg_eps, s.mCredibility))

    prob = LpProblem(g.mName, LpMinimize)

    seg_expr_v = LpAffineExpression(seg_expr_v_tmp)
    seg_expr_h = LpAffineExpression(seg_expr_h_tmp)
    jun_expr = LpAffineExpression(jun_expr_tmp)
    # Define objective function first
    prob += SysSettings.HOST_SEG_LP_COE * seg_expr_h + \
            SysSettings.VIRUS_SEG_LP_COE * seg_expr_v + \
            SysSettings.JUNCTION_LP_COE * jun_expr
    # Define constraints
    cnt = -1
    source_idx, sink_idx = [], []
    for s in g.mSegments:
        if s.isDeadEnd():
            stderr(s, 'is dead')
            continue
        cnt += 1
        if s in g.mSource.values():
            source_idx.append(cnt)
        if s in g.mSink.values():
            sink_idx.append(cnt)
        prob += s.getWeight() + var_seg_eps[cnt] >= var_seg_tar[cnt]
        prob += s.getWeight() - var_seg_eps[cnt] <= var_seg_tar[cnt]
        prob += var_seg_insum[cnt] == var_seg_tar[cnt]
        prob += var_seg_outsum[cnt] == var_seg_tar[cnt]
    for idx, j in enumerate(g.mJunctions):
        prob += j.getWeight() + var_jun_eps[idx] >= var_jun_tar[idx]
        prob += j.getWeight() - var_jun_eps[idx] <= var_jun_tar[idx]
    # fix copy number for source and sink segments.
    # prob += var_seg_tar[g.mSegments.index(g.mSource)] == source_sink_cn
    # prob += var_seg_tar[g.mSegments.index(g.mSink)] == source_sink_cn
    for si in source_idx + sink_idx:
        prob += var_seg_tar[si] == sum(ploidy_dict[g.mSegments[si].mGroup])

    if SysSettings.OUTPUT_LP:
        prob.writeLP("Model.lp")

    # prob.solve(CPLEX())
    prob.solve()

    # For logging purposes
    logs = '\n> Linear Programming Log\n'
    logs += '{0:>12} {1:>12} {2:>12} {4:>12} {3:>12} {5:>7}\n' \
        .format('SEG/JUNC', 'OLD_CN', 'NEW_CN', 'DIFF', 'CRED', 'DIFF%')

    # -- Write to Graph --
    cnt = -1
    for s in g.mSegments:
        if s.isDeadEnd():
            continue
        cnt += 1
        v = var_seg_tar[cnt]
        if not rm_inferred and (not float.is_integer(v.varValue) or v.varValue < 0):
            stdout('\n***Graph can not balanced, try to remove inferred junctions.\n')
            return solveLp(g, ploidy_dict, low_bound=low_bound, rm_inferred=True)
        if not float.is_integer(v.varValue) or v.varValue < 0:
            stdout('\n***Graph still can not balanced, try to discard lower bound\n')
            return solveLp(g, ploidy_dict, low_bound=0, rm_inferred=True)
        diff = v.varValue - s.getWeight()
        if s.mType == 'V':
            cred = SysSettings.VIRUS_SEG_LP_COE * s.mCredibility
        else:
            cred = SysSettings.HOST_SEG_LP_COE * s.mCredibility
        logs += '{0:>12} {1:>12.2f} {2:>12.2f} {3:>12.2f} {4:>12.2f} ' \
            .format(v.name, s.getWeight(), v.varValue, cred, diff)
        if s.getWeight() == 0:
            diff_pct = 'inf'
        else:
            diff_pct = 100 * diff / s.getWeight()
        logs += '{0:6.0f}%\n'.format(float(diff_pct))
        if var_seg_tar[cnt].varValue == 0:
            stderr('Segment [ {} ] has been discarded.'. format(s.getAbbr()))
        s.setLpStat(var_seg_tar[cnt].varValue, cred, float(diff_pct))
        s.setWeight(var_seg_tar[cnt].varValue)

    for idx, j in enumerate(g.mJunctions):
        v = var_jun_tar[idx]
        if not float.is_integer(v.varValue):
            stderr('\n***Graph can not balanced\n')
            return solveLp(g, ploidy_dict, 0)
        diff = v.varValue - j.getWeight()
        cred = SysSettings.JUNCTION_LP_COE * j.mCredibility

        logs += '{0:>12} {1:>12.2f} {2:>12.2f} {3:>12.2f} {4:>12.2f} '\
            .format(v.name, j.getWeight(), v.varValue, cred, diff)
        if j.getWeight() == 0:
            diff_pct = 'inf'
        else:
            diff_pct = 100 * diff / j.getWeight()
        logs += '{0:6.0f}%\n'.format(float(diff_pct))
        # if j.getWeight() == 0:
        #     logs += '{0:6.0f}%\n'.format(float('inf'))
        # else:
        #     logs += '{0:6.0f}%\n'.format(100 * diff / j.getWeight())
        j.setLpStat(var_jun_tar[idx].varValue, cred, float(diff_pct))
        j.setWeight(var_jun_tar[idx].varValue)

    logs += '\nObjective  = {0:.3f}\n'.format(value(prob.objective))
    # logs += 'Confidence = {0:.2f} / 1\n' \
    #     .format(max(0.0001, 1 - avg_diff) ** 0.5)

    # -- Write to Log --
    g.mLog['lp'] = logs

    # returning the total sum of epsilons and weight sum of epsilons
    return value(prob.objective)

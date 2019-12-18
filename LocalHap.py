#!/bin/env python
# -*- coding: utf-8 -*-
import JuncGraph.Graph
from JuncGraph import CompleteGraph, BalanceGraph
from JuncGraph import EulerianCircuit
import argparse
import itertools
from Util.Utilities import *
from collections import namedtuple

__version__ = '1.0.0'


def createGraph(path):
    input_lines = open(path, 'rb').readlines()
    g_tmp = JuncGraph.Graph.Graph()
    g_tmp.readGraph2(input_lines, verbose=True)
    raw_ploidy = g_tmp.getRawPloidy()
    if SysSettings.IF_DELETE_NORMAL_ALLELE:
        '''
        In case of Xm0, no normal copy needs to be removed.
        In case of Xm(0.5X), only major copy needs to be removed.
        In case of Xm(!=0.5X), try to remove minor & major separately.
        '''
        drop_setting = {}
        for group in raw_ploidy:
            major_minor_cn = raw_ploidy[group]
            if major_minor_cn[0] == major_minor_cn[1]:
                _drop = 'M'
            elif major_minor_cn[1] == 0:
                _drop = '-'
            else:
                _drop = SysSettings.DELETE_NORMAL.get(group, 'M')
            drop_setting.update({group: _drop})
        g, lp_obj = _createGraph(input_lines, drop=True, drop_setting=drop_setting)

        #
        # group_order, iter_setting = [], []
        # for group in raw_ploidy:
        #     group_order.append(group)
        #     major_minor_cn = raw_ploidy[group]
        #     '''
        #     In case of Xm0, no normal copy needs to be removed.
        #     In case of Xm(0.5X), only major copy needs to be removed.
        #     In case of Xm(!=0.5X), try to remove minor & major separately.
        #     '''
        #     if major_minor_cn[0] == major_minor_cn[1]:
        #         iter_setting.append(('M'))
        #     elif major_minor_cn[1] == 0:
        #         iter_setting.append(('-'))
        #     else:
        #         if group in SysSettings.DELETE_NORMAL:
        #             iter_setting.append((SysSettings.DELETE_NORMAL[group]))
        #         else:
        #             iter_setting.append(('M', 'm'))
        #
        # min_lp_obj = float('inf')
        # min_lp_obj_g = None
        # for flag_by_group in itertools.product(*iter_setting):
        #     '''
        #     flag_by_group: ('M', '-', 'm'), group_order: [gp1, gp2, gp3]
        #
        #     Meaning that gp1 removes major, do nothing to gp2 and gp3 removes
        #     minor allele as normal allele.
        #     '''
        #     drop_setting = {group_order[gp_idx]: flag for gp_idx, flag in enumerate(flag_by_group)}
        #     _g, _lp_obj = _createGraph(input_lines, drop=True, drop_setting=drop_setting)
        #     # stdout(_g.mLog['lp'])
        #     # stderr('--' * 30)
        #     print 'Mode: Dropping', drop_setting, 'Normalized Obj', _lp_obj
        #     if _lp_obj < min_lp_obj:
        #         min_lp_obj = _lp_obj
        #         min_lp_obj_g = _g
        # g, lp_obj = min_lp_obj_g, min_lp_obj
        print 'chose', lp_obj
    else:
        g, lp_obj = _createGraph(input_lines, drop=False)

    g.printGraphInfo()
    g.backupWeight()
    stdout(g.mLog['lp'])
    stderr('--' * 30)
    if SysSettings.OUTPUT_LP:
        with open(SysSettings.OUTPUT_LP, 'wb') as f:
            f.write(g.mLog['lp'])
    return g

    # if SysSettings.DELETE_NORMAL_ALLELE:
    #     g1, lp_obj1 = _createGraph(path,
    #                                drop_normal=True,
    #                                drop_major=False)
    #     stdout(g1.mLog['lp'])
    #     stderr('--' * 30)
    #     g2, lp_obj2 = _createGraph(path,
    #                                drop_normal=True,
    #                                drop_major=True)
    #     stdout(g2.mLog['lp'])
    #     stderr('--' * 30)
    #     stdout('obj1', lp_obj1, 'obj2', lp_obj2)
    #     g = g1 if lp_obj1 <= lp_obj2 else g2
    # else:
    #     g, lp_obj = _createGraph(path, drop_normal=False)
    # g.printGraphInfo()
    # g.backupWeight()
    # stdout(g.mLog['lp'])
    # stderr('--' * 30)
    # return g


def _createGraph(input_lines, drop=False, drop_setting=None):
    g = JuncGraph.Graph.Graph()
    g.readGraph2(input_lines)
    g.setOrphan()
    CompleteGraph.add_normal_junctions(g, force=False)
    CompleteGraph.reachability_wrapper(g)
    g.setLoop()
    cn_gcd = g.calculateWeight(drop=drop, drop_setting=drop_setting)
    g.setInferredJunctionCred()
    # g.printGraphInfo()
    raw_lp_obj = BalanceGraph.solveLp(g, g.mPloidy)
    return g, raw_lp_obj * cn_gcd


# def _cycle_dict_simplicity(cycles_dict):
#     """Contig dicts with less number of values are considered simpler.
#     """
#     return sum([len(_) for _ in cycles_dict.values()])

def getCycles(g):
    # =========== #
    #  Find LGMs  #
    # =========== #
    if SysSettings.EDGE_SELECTION == 'try_all':
        # 2017.8.29 update
        # Get some forking vertices, we might not be able to find all forks at the first attempt.
        '''

        '''
        _all_fork_pts, _total_num_forks = g.forkingVertices()
        stdout('> All possible forking points: ', *[p.getAbbr() for p in _all_fork_pts])
        # First try default mode to minimize randomness of the output.
        SysSettings.EDGE_SELECTION = 'default'
        try:
            lgm, lgm_cycle_cnt = EulerianCircuit.findLGM(g)
            stdout('> Using default fork search for basic LGM.')
        except:
            import traceback
            print traceback.format_exc()
            stdout('> Default fork search failed, using random search. ', *[p.getAbbr() for p in _all_fork_pts])
            g.refreshGraph()
            # If ImaginaryJunctionReachedException happened, then we need smart_random.
            SysSettings.EDGE_SELECTION = 'smart_random'
            lgm, lgm_cycle_cnt = None, None
            for _ in range(1000):
                # This is hard coded.. lgm might come across
                # ImaginaryJunctionReachedException, so we need to try many times to find
                # a correct LGM
                try:
                    lgm, lgm_cycle_cnt = EulerianCircuit.findLGM(g)
                    break
                except:
                    g.refreshGraph()
            stderr('Tried {} times for finding breakpoints.'.format(_ + 1))
        if not lgm:
            sys.exit('LGM not found, please try again.')
        SysSettings.EDGE_SELECTION = 'try_all'
        stdout('> Starting try_all mode...')
        _fork_pts = []
        for raw_cycle in lgm:
            _fork_pts += EulerianCircuit.getBreakPoints(raw_cycle) # the breakpoints at this lgm.

        # fork_pts = list(set(fork_pts))
        fork_pts = []
        [fork_pts.append(v) for v in _fork_pts if v not in fork_pts]  # ordered dedup
        stdout('> Forking points chosen (first attempt): ', *[p.getAbbr() for p in fork_pts])

        # fork_pts, num_forks = g.forkingVertices()
        # stdout('> Forking points: ', *zip([p.getAbbr() for p in fork_pts], num_forks))
        # if len(fork_pts) > 12 or sum(num_forks) > 30:
        #     stdout('[WARNING] Too many forks to process. Try lowering number of forks.')
        #     SysSettings.EDGE_SELECTION = 'smart_random'
        #     lgm = findLGM(g)
        #     SysSettings.EDGE_SELECTION = 'try_all'
        #     fork_pts = []
        #     for raw_cycle in lgm:
        #         fork_pts += getBreakPoints(raw_cycle)
        #     fork_pts = list(set(fork_pts))  # get it in order
        #     if len(fork_pts) > 12 or sum(num_forks) > 30:
        #         stdout('[WARNING] Still too many forks to process. Try lowering number of forks to 12.')
        #         for i in range(len(fork_pts) - 12):
        #             fork_pts.pop()
        #     stdout('> Forking points: ', *[p.getAbbr() for p in fork_pts])
    else:
        lgm, lgm_cycle_cnt = EulerianCircuit.findLGM(g)
        # Get fork points first. Since they don't change.
        # fork_pts = [g.getStartVertex()]  # 20170910 TODO: why start_vtx already in it?
        _fork_pts = []
        for raw_cycle in lgm:
            _fork_pts += EulerianCircuit.getBreakPoints(raw_cycle)
        fork_pts = []
        [fork_pts.append(v) for v in _fork_pts if v not in fork_pts]  # ordered dedup
        # fork_pts = list(set(fork_pts))  # get it in order
        stdout('> Forking points: ', *[p.getAbbr() for p in fork_pts])

    # ================================== #
    # Break down LGMs into Unique Cycles #
    # ================================== #
    if SysSettings.EDGE_SELECTION != 'try_all':
        cnt_cyc_dict = getCntCycDict(g, lgm, lgm_cycle_cnt, fork_pts)
        stdout('\n> Cycles:')
        print_cnt_cyc_dict(cnt_cyc_dict)
        return [cnt_cyc_dict]
    else:  # If in 'try_all' mode.
        # fork_pts = _all_fork_pts
        print [p.getAbbr() for p in fork_pts]
        # Contig dicts with less number of values are considered simpler.
        _cycle_dict_simplicity = lambda d : sum([len(_) for _ in d.values()])
        while 1:
            try_again = False
            duplication_cnt = 0
            failure_cnt = 0
            fork_choices = []
            for fp in fork_pts:
                # fork_choices = [[(A1,A2,A3), (A1,A3,A2), (A2,A1,A3)...], [(B1,B2), (B2,B1)]]

                fork_choices.append([p for p in itertools.permutations([_i for _i,_ne in enumerate(fp.mNextEdges)
                                                                        if not g.edgeReachesImaginaryJunction(_ne)])])
                EulerianCircuit.FORKING_CHOICE.update({fp: None})
            # Get the number of iterations needed.
            num_of_attempts = 1
            for choice in fork_choices:
                num_of_attempts *= len(choice)
            stdout('> Predicted number of enumerations: {}'.format(num_of_attempts))

            allUniqueCycles = []
            _cnt = 0

            for fork_comb in itertools.product(*fork_choices):
                _cnt += 1
                if _cnt > SysSettings.MAX_BRUTE_FORCE_LOOPS:
                    break
                stderr('\r> Trying combinations of forks: {}'.format(_cnt), end='')
                g.refreshGraph()
                for i, choice in enumerate(fork_comb):
                    EulerianCircuit.FORKING_CHOICE[fork_pts[i]] = choice
                try:
                    lgm, lgm_cycle_cnt = EulerianCircuit.findLGM(g, fork_pts)
                except ImaginaryJunctionReachedException:
                    stdout('Imaginary junction reached for current choice of fork.')
                    failure_cnt += 1
                    print '-' * 100
                    continue
                except TooManyTraversalException:
                    stdout('Too many traversals for some vertices for '
                           'current choice of fork.')
                    print '-' * 100
                    failure_cnt += 1
                    continue
                except BreakPointNotEnumeratedException as ex:
                    fork_pts.append(ex.breakPoint)
                    try_again = True
                    stdout('\n> Updated Forking points: ', *[p.getAbbr() for p in fork_pts])
                    break

                # When LGM is successfully found.
                try:
                    cnt_cyc_dict = getCntCycDict(g, lgm, lgm_cycle_cnt, fork_pts)  # break down LGM
                except MixedGroupException as e:
                    failure_cnt += 1
                    stdout(str(e))
                    continue
                is_new_lgm = True
                for cyc_d in allUniqueCycles:
                    if compareCycles(cnt_cyc_dict, cyc_d):
                        is_new_lgm = False
                        break
                if is_new_lgm:
                    if SysSettings.CHOOSE_SIMPLE_CONTIG:
                        # if _cnt % 10 == 0:
                        #     stderr('tried {} combinations of forks'.format(_cnt))
                        # Contig dicts with less number of values are considered simpler.
                        if not allUniqueCycles:
                            allUniqueCycles.append(cnt_cyc_dict)
                        elif _cycle_dict_simplicity(cnt_cyc_dict) < _cycle_dict_simplicity(allUniqueCycles[0]):
                            allUniqueCycles = [cnt_cyc_dict]
                        elif _cycle_dict_simplicity(cnt_cyc_dict) == _cycle_dict_simplicity(allUniqueCycles[0]):
                            if len(allUniqueCycles) < SysSettings.LIMIT_SIMPLE_CONTIG:
                                allUniqueCycles.append(cnt_cyc_dict)
                        else:
                            # Do nothing, since the new cnt_cyc_dict has more unique cycles.
                            pass
                    else:
                        stdout('\n> Cycles found:')
                        print_cnt_cyc_dict(cnt_cyc_dict)
                        allUniqueCycles.append(cnt_cyc_dict)
                else:
                    duplication_cnt += 1
            stderr()
            if try_again:
                continue
            else:
                if not SysSettings.CHOOSE_SIMPLE_CONTIG:
                    # allUniqueCycles might be empty in "choose_simple_contig" mode,
                    # causing below result erroneous. number unique cycle compositions
                    # will be higher.
                    stdout('\n> Number of failed attempt: {}'.format(failure_cnt))
                    stdout('> Number of duplicate cycle compositions: {}'.format(duplication_cnt))
                    stdout('> Number of unique cycle compositions: {}'.format(num_of_attempts - duplication_cnt))
                return allUniqueCycles


def getCntCycDict(g, lgm, cyc_cnt, fork_pts):
    breakdown = []
    breakdown_cnt = []
    for _lgm_idx, c in enumerate(lgm):
        for br_c in EulerianCircuit.breakdownLGM(c, fork_pts):
            '''
            Mixed group exception has to be caught here, since composite cycle with
             more than one group can be broken down into multiple cycles that do not
             cause conflict.
            '''
            if len(set([v.getGroup() for v in br_c if v.mType == 'H'])) > 1:
                raise MixedGroupException(str([v.getAbbr() for v in br_c]))
            # print [v.getAbbr() for v in br_c]
            for i, v in enumerate(br_c[1:-1]):
                if g.isSource(v.mSeg):
                    # change < H2+ H3+ H1+ H2+ > to < H1+ H2+ H3+ H1+ >
                    br_c = br_c[i:] + br_c[1:i] + [v]
                    break
            breakdown.append(br_c)
            breakdown_cnt.append(cyc_cnt[_lgm_idx])

    # Get the dict to show the copy number of each unique unit cycles.
    unique_cycles = []  # unique cycles in breakdown list
    unique_cycles_cnt = []
    for i, c in enumerate(breakdown):
        try:
            c_idx = unique_cycles.index(c)
            unique_cycles_cnt[c_idx] += breakdown_cnt[i]
        except ValueError:
            # not in unique_cycles
            unique_cycles.append(c)
            unique_cycles_cnt.append(breakdown_cnt[i])

    cnt_cyc_dict = {}
    for cnt, c in zip(unique_cycles_cnt, unique_cycles):
        if cnt in cnt_cyc_dict:
            cnt_cyc_dict[cnt].append(c)
        else:
            cnt_cyc_dict.update({cnt: [c]})

    return cnt_cyc_dict

    # [unique_cycles.append(c) for c in breakdown if c not in unique_cycles]  # cannot use set()
    # cnt_cyc_dict = {}
    # for c in unique_cycles:
    #     cycle_cnt = breakdown.count(c)
    #     if cycle_cnt in cnt_cyc_dict:
    #         cnt_cyc_dict[cycle_cnt].append(c)
    #     else:
    #         cnt_cyc_dict.update({cycle_cnt: [c]})
    # return cnt_cyc_dict


def compareCycles(cycles_dict1, cycles_dict2):
    """Receive two cycles dictionary, tell if they are the same"""
    # Pre-check
    if set(cycles_dict1.keys()) != set(cycles_dict2.keys()):
        return False
    for cnt in cycles_dict1:
        # if have different number of vertices per cycles per copy number.
        if set([len(c) for c in cycles_dict1[cnt]]) != set([len(c) for c in cycles_dict2[cnt]]):
            return False
        # if same cycle length are the same but comprised of different vertices
        cycle1_v_set = [set(c) for c in cycles_dict1[cnt]]
        for c in cycles_dict2[cnt]:
            if set(c) not in cycle1_v_set:
                return False
    # Compare
    for cnt in cycles_dict1:
        for cyc1 in cycles_dict1[cnt]:
            cyc1_extended = cyc1 + cyc1[1:-1]  # circular extension
            cyc1_extended_str = ''.join([v.getAbbr() for v in cyc1_extended])
            match = False
            for cyc2 in cycles_dict2[cnt]:
                if len(cyc2) == len(cyc1) and ''.join([v.getAbbr() for v in cyc2]) in cyc1_extended_str:
                    match = True
            if not match:
                return False
    return True


def consumedContigs(new_cycles, ori_cycles):
    used_contig = {}
    for o_cn in ori_cycles:
        for o_cyc in ori_cycles[o_cn]:
            found = False
            for cn in new_cycles:
                if o_cyc in new_cycles[cn]:
                    found = True
                    if o_cn > cn:
                        used_contig.update({tuple(o_cyc): o_cn - cn})
            if not found:
                used_contig.update({tuple(o_cyc): o_cn})
    return used_contig


CYCLES_PROCESSED_CNT = 0


def constructAlleles(g, cycles_dict, cycles_dict_backup):
    global CYCLES_PROCESSED_CNT
    CYCLES_PROCESSED_CNT += 1
    alleles, contig_per_allele = [], []
    stdout('\n> Reconstructing allele with specified copy number profile...')
    stdout('\n> Processing Cycles {}:'.format(CYCLES_PROCESSED_CNT))
    print_cnt_cyc_dict(cycles_dict)
    for gp in g.mGroupPriority:
        ploidy = g.mPloidy[gp]
        for idx, cn in enumerate(ploidy):
            prev_cycles_dict = {}
            [prev_cycles_dict.update({_cn: [_cyc for _cyc in cycles_dict[_cn]]}) for _cn in cycles_dict]
            merged_lgm = EulerianCircuit.findWeightedCycle(g, cycles_dict, cn, g.getGroupStartVertex(gp),
                                                           last_try=len(ploidy) == idx + 1)
            contig_per_allele.append(consumedContigs(cycles_dict, prev_cycles_dict))
            # Not adding if cn ==0, leaving contigs_per_allele with a useless {} ending.
            if cn > 0: alleles.append({cn: merged_lgm})
            stdout('\n  Cycles left:')
            print_cnt_cyc_dict(cycles_dict)
    if not any([cycles_dict[cn] for cn in cycles_dict]):
        stdout('\n  Result')
        return alleles, contig_per_allele

    stdout('\nFailed at provided copy number profile. Enter safe mode...')
    cycles_dict = cycles_dict_backup
    if SysSettings.EQUALLY_ASSIGN_UCYC_TO_ALLELES:
        alleles, contig_per_allele = \
            EulerianCircuit.equallyAssignCyclesToAllele(g, cycles_dict)
    else:
        alleles, contig_per_allele = [], []
        for gp in g.mGroupPriority:
            ploidy = g.mPloidy[gp]
            for i in range(sum(ploidy)):
                prev_cycles_dict = {}
                # just a deep copy
                [prev_cycles_dict.update(
                    {_cn: [_cyc for _cyc in cycles_dict[_cn]]}
                ) for _cn in cycles_dict]

                merged_lgm = EulerianCircuit.findWeightedCycle(
                    g, cycles_dict, 1, g.getGroupStartVertex(gp)
                )
                contig_per_allele.append(
                    consumedContigs(cycles_dict, prev_cycles_dict)
                )
                stdout('\n  Cycles left:')
                print_cnt_cyc_dict(cycles_dict)
                alleles.append({1: merged_lgm})
    stdout('\n  Result')
    return alleles, contig_per_allele

# G = None
# def foo():
#     return getCycles(G)

def main():
    params = readParam()
    sample_path = params.input_file
    output_path = params.output_file
    stderr('-' * 15, 'Processing sample %s' % sample_path, '-' * 15, '\n')
    # CONTIGS Format: {(A,B,C,A): serial_number}
    CONTIGS = {}  # Serial number for UCYC1~n naming, shared between solutions.
    RESULTS = []
    # FREE_VIRUS: if any contig in CONTIGS belongs to defined free virus.
    FREE_VIRUS = {}
    _CONTIG_CNT = 1

    # Create the graph object
    g = createGraph(sample_path)
    # Add dropped normal allele into contigs list.
    dropped_normal_allele = {} # {cyc_tuple: cn}
    for gp in g.mDeleteNormalCNByGroup:
        contig_tuple = tuple([_.mPosV for _ in g.mDeletedNormalSegList[gp]])
        if g.mDeleteNormalCNByGroup[gp] > 0:
            CONTIGS.update({contig_tuple: _CONTIG_CNT})
            _CONTIG_CNT += 1
            dropped_normal_allele.update({contig_tuple: g.mDeleteNormalCNByGroup[gp]})

    # global G
    # G = g
    # import cProfile
    # cProfile.run('foo()')
    # sys.exit()

    all_cycles_dict = getCycles(g) # [{cyc1_cnt: cyc1, cyc2_cnt:cyc2}, {...}]
    # For each solution
    for index, cycles_dict in enumerate(all_cycles_dict):
        # Back up all cycles
        cycles_dict_backup = {}
        [cycles_dict_backup.update({cn: [cyc for cyc in cycles_dict[cn]]}) for cn in cycles_dict]
        # Prepare for output HEADER and all_contigs column
        all_contigs = {}  # {serial_no: cn}
        for _t in dropped_normal_allele:
            all_contigs.update({CONTIGS[_t]: dropped_normal_allele[_t]})
        for cn in cycles_dict:
            contig_list = cycles_dict[cn]
            for contig in contig_list:
                contig_tuple = tuple(contig)
                if contig_tuple not in CONTIGS:
                    CONTIGS.update({contig_tuple: _CONTIG_CNT})
                    _CONTIG_CNT += 1
                # Include droppped normal copy count.
                _dropped_n_cn = all_contigs.get(CONTIGS[contig_tuple], 0)
                # This is for header, so need to multiply by g.mCNDividedBy
                all_contigs.update({CONTIGS[contig_tuple]: _dropped_n_cn + cn * g.mCNDividedBy})  # used to be only "cn"
        all_contigs_str = ','.join(['UCYC{}:{}'.format(k, all_contigs[k]) for k in sorted(all_contigs.keys())])

        # Try to merge unit cycles into alleles according to target CN profile.
        alleles, contig_per_allele = constructAlleles(g, cycles_dict, cycles_dict_backup)

        # Check if any contig is a free virus
        for contig in CONTIGS:
            if all(v.mType == 'V' for v in contig):
                _sorted_contig_str_1 = tuple(g.sortVirusGenome(','.join(v.getAbbr()
                                                                        for v in contig)))
                _sorted_contig_str_2 = tuple(g.sortVirusGenome(','.join(v.getComplementAbbr()
                                                                        for v in reversed(contig))))
                free_virus_ctg_id = None
                if _sorted_contig_str_1 in g.mFreeVirusGenomeTemplate:
                    free_virus_ctg_id = g.mFreeVirusGenomeTemplate[_sorted_contig_str_1]
                if _sorted_contig_str_2 in g.mFreeVirusGenomeTemplate:
                    free_virus_ctg_id = g.mFreeVirusGenomeTemplate[_sorted_contig_str_2]
                if free_virus_ctg_id != None:
                    FREE_VIRUS.update({contig: free_virus_ctg_id})

        # Prepare to dump out the results.
        _res = namedtuple('result', ['sol_no', 'group', 'sol_ucyc',
                                     'allele_no', 'allele_cn',
                                     'allele_ucyc', 'free_virus_ucyc',
                                     'notes'])

        for idx, lgm_dict in enumerate(alleles):
            stdout('\n  Allele {}, CN: {}'.format(idx + 1, lgm_dict.keys()[0]))
            printLGM([lgm_dict.values()[0]])
            allele_cn = lgm_dict.keys()[0]
            # Divide by allele_cn to get cn of cycles for a single allele.
            tmp = [(CONTIGS[cyc_t], contig_per_allele[idx][cyc_t] / allele_cn, cyc_t)
                   for cyc_t in contig_per_allele[idx]]
            # If allele is already divided by a gcd in certain settings.
            actual_allele_cn = allele_cn * g.mCNDividedBy
            # Sort cycle report order by its serial number.
            _free_v_tmp = []
            _ctg_by_allele_tmp = []
            for t in sorted(tmp, key=lambda x: x[0]):  # sort by serial number of UCYC
                _s = 'UCYC{}:{}'.format(t[0], t[1])
                if t[2] in FREE_VIRUS:
                    _free_v_tmp.append(_s)
                else:
                    _ctg_by_allele_tmp.append(_s)
            free_virus_str = ','.join(_free_v_tmp) if _free_v_tmp else 'N/A'
            contig_by_allele_str = ','.join(_ctg_by_allele_tmp)

            # Determine the group the allele belongs to
            if lgm_dict.values()[0]:
                _group = lgm_dict.values()[0][0].getGroup()  # group of the first vertex
            else:
                _group = 'N/A'  # empty allele, should not happen
            if not _group:
                _group = '1'

            # Organize output.
            RESULTS.append(_res(str(index+1),
                                _group,
                                all_contigs_str,
                                str(idx + 1),
                                str(actual_allele_cn),
                                contig_by_allele_str,
                                free_virus_str,
                                notes=()))
        for gp in g.mDeleteNormalCNByGroup:
            idx += 1
            removed_cn = g.mDeleteNormalCNByGroup[gp]
            if removed_cn > 0:
                ucyc_id = CONTIGS[tuple([_.mPosV for _ in g.mDeletedNormalSegList[gp]])]
                RESULTS.append(_res(str(index + 1),
                                    gp if gp else '1',
                                    all_contigs_str,
                                    str(idx + 1),
                                    str(removed_cn),
                                    'UCYC{}:1'.format(ucyc_id),
                                    'N/A',  # no virus is involved in normal copy.
                                    notes=('REF_ALLELE', )))

    if SysSettings.EQUALLY_ASSIGN_UCYC_TO_ALLELES:
        print("\n[IMPORTANT]\n"
              "You are using --equally_assign_cycles_to_alleles option, "
              "which is an experimental feature and the above 'Result' "
              "does not reflect the actual appearance of possible LGMs. "
              "Please refer to the result file at:"
              "\n\n{}\n".format(output_path))

    dumpResult(output_path, CONTIGS, RESULTS, FREE_VIRUS, g)


def dumpResult(path, contigs, results, free_virus, graph):
    with open(path, 'wb') as of:
        # Write general header
        _doc_template = '@DOC\t{}\n'
        _doc_fields = []
        _doc_fields.append('VERSION={}'.format(__version__))
        _doc_fields.append('SAMPLE_ID={}'.format(graph.mName))
        ploidy_strs = []
        for gp in graph.mPloidyRaw:
            gp_name = gp if gp else '1'
            ploidy_strs.append('GROUP{}:{}M:{}m'.format(gp_name, *graph.mPloidyRaw[gp]))
        _doc_fields.append('PLOIDY={}'.format(','.join(ploidy_strs)))
        _doc_fields.append('ARGS={}'.format(sys.argv[1:]))
        of.write(_doc_template.format(';'.join(_doc_fields)))

        # Write V and H segments in header
        h_segs_header, v_segs_header = [], []
        for s in graph.mSegments:
            _header_fields = []
            _header_fields.append('ID={}'.format(s.getAbbr()))
            if s.mType == 'H':
                _header_fields.append('GROUP={}'.format(s.mGroup if s.mGroup else '1'))

            _header_fields.append('CN_ORIG={:.2f}'.format(2.0 * s.mCov.mCovBackup / graph.mAvgCov))

            if s.mType == 'H':
                _weight_change = 2.0 * (s.mCov.mCovBackup - s.mCov.mCov) / graph.mAvgCov
                _header_fields.append('CN_REF={:.2f}'.format(_weight_change))
            lp_stat = s.getLpStat()
            _header_fields.append('CN_ALT={:.2f}'.format(lp_stat.orig_cn * graph.mCNDividedBy))

            if lp_stat.lp_cn == 'N/A':
                _alt_lp = 'N/A'
                _lp_w = 'N/A'
                _lp_offset = 'N/A'
            else:
                _alt_lp = '{:.2f}'.format(lp_stat.lp_cn * graph.mCNDividedBy)
                _lp_w = '{:.2f}'.format(lp_stat.lp_cred)
                _lp_offset = '{:.2f}%'.format(lp_stat.diff_pct)

            _header_fields.append('CN_ALT_LP={}'.format(_alt_lp))
            _header_fields.append('LP_WEIGHT={}'.format(_lp_w))
            _header_fields.append('LP_OFFSET%={}'.format(_lp_offset))

            if graph.isSource(s):
                _header_fields.append('NOTES=SOURCE')
            if graph.isSink(s):
                _header_fields.append('NOTES=SINK')
            for _tag in s.mTags:
                _header_fields.append('{}={}'.format(_tag, s.mTags[_tag]))
            _h_template = '@HSEG\t{}\n' if s.mType == 'H' else '@VSEG\t{}\n'
            header_str = _h_template.format(';'.join(_header_fields))
            if s.mType == 'H':
                h_segs_header.append(header_str)
            else:
                v_segs_header.append(header_str)
        [of.write(_) for _ in h_segs_header]
        [of.write(_) for _ in v_segs_header]
        # Write Junctions into header
        # @JUNC   ID=JUNC1;GROUP=1;SEGLINK:H1+,V2+;DROP_NORMAL_CN=xxx;
        # ORIG_CN=xxx;LP_CN=xxx;LP_WEIGHT=xxx;LP_OFFSET%=xxx%;TYPE=[REAL/INFER]
        junc_headers = []
        for i, j in enumerate(graph.mJunctions):
            if j.mInferred or j.isImaginaryJunction():
                continue
            _header_fields = []
            _header_fields.append('ID={}'.format(j.JUNCID))
            _header_fields.append('SEGLINK={},{}'.format(j.mEdgeA.mSourceV.getAbbr(),
                                                         j.mEdgeA.mTargetV.getAbbr()))
            _header_fields.append('CN_ORIG={:.2f}'.format(2.0 * j.mCov.mCovBackup / graph.mAvgCov))
            _weight_change = 2.0 * (j.mCov.mCovBackup - j.mCov.mCov) / graph.mAvgCov
            if _weight_change > 0:
                _header_fields.append('CN_REF={:.2f}'.format(_weight_change))

            lp_stat = j.getLpStat()
            _header_fields.append('CN_ALT={:.2f}'.format(lp_stat.orig_cn * graph.mCNDividedBy))

            if lp_stat.lp_cn == 'N/A':
                _alt_lp = 'N/A'
                _lp_w = 'N/A'
                _lp_offset = 'N/A'
            else:
                _alt_lp = '{:.2f}'.format(lp_stat.lp_cn * graph.mCNDividedBy)
                _lp_w = '{:.2f}'.format(lp_stat.lp_cred)
                _lp_offset = '{:.2f}%'.format(lp_stat.diff_pct)

            _header_fields.append('CN_ALT_LP={}'.format(_alt_lp))
            _header_fields.append('LP_WEIGHT={}'.format(_lp_w))
            _header_fields.append('LP_OFFSET%={}'.format(_lp_offset))

            # _header_fields.append('TYPE={}'.format('INFERRED' if j.mInferred else 'REAL'))
            for _tag in j.mTags:
                _header_fields.append('{}={}'.format(_tag, j.mTags[_tag]))
            _h_template = '@JUNC\t{}\n'
            header_str = _h_template.format(';'.join(_header_fields))
            junc_headers.append(header_str)
        [of.write(_) for _ in junc_headers]

        # Write UnitCycles into header
        tmp_list = []
        for contig_t in contigs:
            tmp_list.append((contigs[contig_t], contig_t))
        for contig_num, contig_t in sorted(tmp_list, key=lambda x: x[0]):
            _header_fields = []
            _header_fields.append('ID=UCYC{}'.format(contig_num))
            _header_fields.append('SEGSTR={}'.format(cycleToStr(contig_t)))
            if contig_t in free_virus:
                _header_fields.append('FVGM_ID={}'.format(free_virus[contig_t]))
            _h_template = '@UCYC\t{}\n'
            of.write(_h_template.format(';'.join(_header_fields)))

        # Write Solution
        _s_template = 'SOLUTION\t{}\n'
        for res in results:
            #'sol_no', 'group', 'sol_ucyc', 'allele_no', 'allele_cn', 'allele_ucyc'
            _s_fields = ['NO={}'.format(res.sol_no),
                         'GROUP={}'.format(res.group),
                         'SOLUTION_UCYC={}'.format(res.sol_ucyc),
                         'ALLELE_NO={}'.format(res.allele_no),
                         'ALLELE_CN={}'.format(res.allele_cn),
                         'ALLELE_UCYC={}'.format(res.allele_ucyc),
                         'ALLELE_FVGM_UCYC={}'.format(res.free_virus_ucyc),
                         ]
            if res.notes:
                _s_fields.append('NOTES={}'.format(';'.join(res.notes)))
            if int(res.allele_cn) > 0:
                of.write(_s_template.format(';'.join(_s_fields)))


def readParam():
    parser = argparse.ArgumentParser(
        description='Reconstructing the local haplotype surrounding '
                    'virus integrated regions \non human genome. '
                    '\n\nVersion 1.0.1 (9-JUN-2018)',
        epilog='Details of file formats and other info will be posted'
               'on our GitHub page.\n\n'
               'Developed by:  \n'
               '  Chang Xu (xuchang0310@gmail.com)   \n'
               '  Wenlong Jia (wenlongkxm@gmail.com) \n',
        formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, max_help_position=30))

    parser.add_argument(dest='input_file',
                        help='Input file. Format requirements can be found \n'
                             'on our GitHub page.')

    parser.add_argument('-o', dest='output_file', default=None,
                        help='Output result to file. Visit our GitHub page \n'
                             'for explanation of output formats. \n'
                             '[<input_file>.ucyc.txt]')
    # === Weight Settings ===
    # Notice that the weight settings through params are overwritten if set in input file.
    parser.add_argument('--host_seg_weight_factor',
                        dest='host_seg_weight_factor',
                        default=1, type=float,
                        help='Weight of host segments. Linear programming \n'
                             'coefficient of all host segments are scaled \n'
                             'with this factor. 1 to disable scaling. [1.0]',
                        metavar='')

    parser.add_argument('--virus_seg_weight_factor',
                        dest='virus_seg_weight_factor',
                        default=1, type=float,
                        help='Weight of virus segments. Linear programming \n'
                             'coefficient of all virus segments are scaled \n'
                             'with this factor. 1 to disable scaling. [1.0]',
                        metavar='')

    parser.add_argument('--junction_weight_factor',
                        dest='junction_weight_factor',
                        default=1, type=float,
                        help='Weight of junctions. Linear programming \n'
                             'coefficient of all junctions are scaled \n'
                             'with this factor. 1 to disable scaling. [1.0]',
                        metavar='')

    parser.add_argument('--cap_sv_weight',
                        dest='cap_sv_weight',
                        default=1, type=float,
                        help='Depending on the detection rate of upstream \n'
                             'software, you can cap the LP weight of \n'
                             'junctions of structural variations on same \n'
                             'reference. 1 to disable capping. [1.0]',
                        metavar='')

    parser.add_argument('--cap_integ_weight',
                        dest='cap_integ_weight',
                        default=1, type=float,
                        help='Depending on the detection rate of upstream \n'
                             'software, you can cap the LP weight of \n'
                             'junctions supporting virus integrations. \n'
                             '1 to disable capping. [1.0]',
                        metavar='')

    parser.add_argument('--inferred_jun_weight',
                        dest='inferred_junction_weight',
                        default=0.01, type=float,
                        help='Set the LP weight of inferred junctions. '
                             '[0.01]',
                        metavar='')
    # parser.add_argument('-isw', '--inferred_seg_weight', dest='inferred_seg_weight',
    #                     default=1, type=float,
    #                     help='Lower the coverage of inferred segments. [1]')

    parser.add_argument('--cap_jun_weight', dest='cap_jun_weight',
                        type=float, default=1.0,
                        help='Cap junctions\' weight to reflect the effect \n'
                             'of sequencing coverage fluctuations. [1.0]',
                        metavar='')

    parser.add_argument('--auto_weight', dest='auto_weight',
                        default=False, action='store_true',
                        help='If --auto_weight option is set,               \n'
                             '--weight_by_length and --weight_by_cn         \n'
                             'will become effective. All "W" tags in input  \n'
                             'file will be overridden. At the same time,    \n'
                             'you may use other parameters to further tune  \n'
                             'weights, when multiple weighing rules are     \n'
                             'imposed, the program takes the lowest one.    \n'
                             '[False]')

    parser.add_argument('--weight_by_length', dest='weight_by_length',
                        default='0-50:0.2,51-100:0.4,101-150:0.6,'
                                '151-350:0.8,351-:1',
                        help='Assign different LP weight for segments of    \n'
                             'different length. Longer segments should have \n'
                             'larger LP weight since the calculation of     \n'
                             'their weight is more reliable. This option    \n'
                             'also affects the weighing of junctions        \n'
                             'connecting segments of varying lengths.       \n'
                             'format: \n'
                             ' --weight_by_length 0-50:0.2,51-100:0.4,\n'
                             ' 101-150:0.6,151-350:0.8,351-:1',
                        metavar='')

    parser.add_argument('--weight_by_cn', dest='weight_by_cn',
                        default='-100:1,101-200:0.9,201-300:0.8,301-400:0.7,401-500:0.6,501-:0.5',
                        help='Assign different LP weight for segments or    \n'
                             'junctions of different copy number. Segments  \n'
                             'or junctions with high copy number are more   \n'
                             'likely to be changed by LP.                   \n'
                             'format: \n'
                             ' --weight_by_cn 0-100:1,101-200:0.9,201-300:  \n'
                             '0.8,301-400:0.7,401-500:0.6,501-:0.5',
                        metavar='')

    # parser.add_argument('--override_w_tag', dest='override_w_tag', default=False, action='store_true',
    #                     help='If "W" tag has been given in the input file, then you need to add '
    #                          '--override_weight argument to override the weight setting in input file '
    #                          'so that --weight_by_length and --weight_by_depth can take effect.\n')

    # === Graph Construction and Traversal ===
    parser.add_argument('--no_infer_normal', dest='no_infer_normal',
                        action='store_true', default=False,
                        help='By default, normal junctions are inferred     \n'
                             'based on copy number of existing junctions    \n'
                             'and surrounding segments. Use this option     \n'
                             'to disable the feature. [False]')

    parser.add_argument('--keep_ref_allele', dest='keep_ref_allele',
                        action='store_true', default=False,
                        help='Keep the allele which is considered a reference\n'
                             'allele by default. In the mean time, LP will   \n'
                             'be solved with original copy number. [False] \n')

    parser.add_argument('--lp_with_original_ploidy',
                        dest='lp_with_original_ploidy',
                        action='store_true', default=False,
                        help='Do linear programming with original ploidy  \n'
                             'setting instead of lowering it.             \n'
                             'i.e. do not change 2m0 into 1m0 for better  \n'
                             'LP performance. [False]')

    parser.add_argument('--equally_assign_cycles_to_alleles',
                        action='store_true', default=False,
                        help='(EXPERIMENTAL) \n'
                             'When multiple alleles are found for a sample  \n'
                             '(and the alleles are different), try to assign\n'
                             'unit cycles to each allele with equal number  \n'
                             'of copies. Used with --drop_ref_allele and    \n'
                             '--lp_with_original_ploidy [False]')

    parser.add_argument('--drop_ref_allele', dest='drop_ref_allele',
                        default='M',
                        help='Drop major ("M") or minor allele ("m") as the \n'
                             'normal copy. \n'
                             'format: \n'
                             '  > If more than one group is involved.   \n'
                             '      --drop_normal GROUP1:M,GROUP2:m,... \n'
                             '  > If only default group is involved.    \n'
                             '      --drop_normal M/m ',
                        metavar='')

    parser.add_argument('-e', dest='no_toleration',
                        default=False, action='store_true',
                        help='Upon failure at solving LP problem, do not \n'
                             'try to recover by deleting incompatible    \n'
                             'edges. [False]')

    parser.add_argument('--linear_virus', dest='linear_virus',
                        default=False, action='store_true',
                        help='If the virus concerned is linear in stead of \n'
                             'circular. [False]')

    parser.add_argument('-s', '--select_edge', dest='edge_selection_mode',
                        default=4, type=int,
                        help='Mode of edge selection [4]: \n'
                             '1.conservative \n'
                             '2.random \n'
                             '3.random_with_memory \n'
                             '4.try_all\n'
                             '',
                        metavar='')

    # === Output ===
    parser.add_argument('-v', dest='verbose',
                        default=False, action='store_true',
                        help='Verbose mode. [False]')

    parser.add_argument('-m', dest='max_loop', default=100, type=int,
                        help='Maximum number of loops while using "try_all"\n'
                             'option to select edges. [100]',
                        metavar='')

    parser.add_argument('-S', dest='choose_simple_contig',
                        default=False, action='store_true',
                        help='Only output simple contig composition. [False]')

    parser.add_argument('-L', dest='limit_simple_contig',
                        default=-1, type=int,
                        help='Limit number of simple contigs to be output. \n'
                             'Input a number less than 1 to disable. [-1]',
                        metavar='')
    parser.add_argument('--write_lp',
                        dest='write_lp', default=None,
                        help='Write the given LP problem to a .lp file.',
                        metavar='')

    # parser.add_argument('--print_result',
    #                     dest='print_result', default=False,
    #                     help='Print one of the solutions constructed with '
    #                          'unique cycles. By default, only unique cycles '
    #                          'and their copy number are displayed. [False]')

    '''
    infer normal junction happens before and during ensuring reachability.
    1. if infer_normal function is enabled, the program will infer normal junctions 
        before ensuring reachability.
    2. reachability is a must, certain normal junctions have to be created, regardless 
        if below argument is set True.
    '''

    params = parser.parse_args()
    selection_modes = {1: 'default',
                       2: 'random',
                       3: 'smart_random',
                       4: 'try_all'}
    SysSettings.VERBOSE = params.verbose
    SysSettings.EDGE_SELECTION = selection_modes[params.edge_selection_mode]
    SysSettings.DELETE_INCOMPATIBLE_EDGES = not params.no_toleration

    # === Weight ===
    SysSettings.HOST_SEG_LP_COE = params.host_seg_weight_factor
    SysSettings.VIRUS_SEG_LP_COE = params.virus_seg_weight_factor
    SysSettings.JUNCTION_LP_COE = params.junction_weight_factor

    SysSettings.INFERRED_JUNCTION_WEIGHT = params.inferred_junction_weight
    # SysSettings.INFER_COV_COE = params.inferred_seg_weight
    SysSettings.CIRCULAR_VIRUS = not params.linear_virus
    SysSettings.MAX_BRUTE_FORCE_LOOPS = params.max_loop
    SysSettings.CHOOSE_SIMPLE_CONTIG = params.choose_simple_contig
    SysSettings.LIMIT_SIMPLE_CONTIG = params.limit_simple_contig if params.limit_simple_contig > 0 else float('inf')
    SysSettings.OUTPUT_LP = params.write_lp
    SysSettings.INFER_NORMAL_JUNCTIONS = not params.no_infer_normal
    SysSettings.IF_DELETE_NORMAL_ALLELE = not params.keep_ref_allele
    SysSettings.LP_WITH_ORIGINAL_PLOIDY = params.lp_with_original_ploidy
    SysSettings.CAP_JUNCTION_WEIGHT = params.cap_jun_weight
    SysSettings.AUTO_WEIGHING = params.auto_weight
    SysSettings.SV_WEIGHT = params.cap_sv_weight
    SysSettings.INTEG_WEIGHT = params.cap_integ_weight
    if SysSettings.AUTO_WEIGHING:
        # format: [(0, 100, 0.1), (101, 1000, 0.5), ...] or None
        SysSettings.WEIGHT_BY_LENGTH = _read_range_weight_format(params.weight_by_length)
        # format: [(0, 100, 0.1), (101, 1000, 0.5), ...] or None
        SysSettings.WEIGHT_BY_CN = _read_range_weight_format(params.weight_by_cn)
    SysSettings.DELETE_NORMAL = _read_drop_normal_settings(params.drop_ref_allele)
    # SysSettings.OVERRIDE_WEIGHT = params.override_w_tag
    # if auto_weight is set, W tags are not effective anymore.
    # SysSettings.OVERRIDE_WEIGHT = params.auto_weight
    if not params.output_file:
        params.output_file = params.input_file.strip() + '.ucyc.txt'
        # if os.path.exists(params.output_file):
        #     pass
    # SysSettings.PRINT_RESULT = params.print_result

    if params.equally_assign_cycles_to_alleles:
        if params.keep_ref_allele or not params.lp_with_original_ploidy:
            print('Cannot use --equally_assign_cycles_to_alleles while keeping'
                  'reference alleles or not calculating LP with original '
                  'ploidy')
            exit(2)
        SysSettings.EQUALLY_ASSIGN_UCYC_TO_ALLELES = True
    return params


def _read_range_weight_format(param_str):
    '''For reading format like:
        -end1:value1,start2-end2:value2,.start3-:value3
    '''
    if not param_str:
        return None
    w_by_len = []
    for x in param_str.split(','):
        range_weight = x.split(':')
        st_ed = range_weight[0].split('-')
        if not st_ed[1]:  # end, 20-30,31-100,10-
            _range = int(st_ed[0]), float('inf')
        elif not st_ed[0]:  # start, -30,31-100
            _range = 0, int(st_ed[1])
        else:
            _range = map(int, st_ed)
        _weight = float(range_weight[1])
        w_by_len.append((_range[0], _range[1], _weight))
    return w_by_len


def _read_drop_normal_settings(param_str):
    d = {}
    if not param_str:
        return d
    if len(param_str.split(',')) > 1:
        for gp_al in param_str.split(','):
            gp, al = gp_al.split(':')
            if al not in ('M', 'm'):
                sys.exit('Wrong --drop_normal argument {}. '.format(al))
            d.update({gp: al})
    else:
        if param_str not in ('M', 'm'):
            sys.exit('Wrong --drop_normal argument {}. '.format(param_str))
        d.update({None: param_str})
    return d


if __name__ == '__main__':
    main()


from Components import *
from Util.Utilities import *
from fractions import gcd
import SysSettings
from collections import OrderedDict


class Graph:
    def __init__(self):
        self.mName = 'test'
        self.mLog = {'inferred_junctions': [],
                     'lp': None}
        self.mSegments = []
        self.mJunctions = []
        self.mCirclesFound = []
        self.mAvgCov = None
        self.mAvgCovRaw = None
        self.mSource = {}  # segment  {Group: sourceSegment}
        self.mSink = {}  # segment
        self.mMaxVIdx = 0
        self.mMinVIdx = 100000000
        self.mPurity = 1.0
        self.mAvgPloidy = 2.0
        self.mPloidy = {}  # two alleles, each with one copy {group: [allele1, allele2]}
        self.mPloidyRaw = {}
        self.mGroupPriority = [None]

        self.mCNDividedBy = 1
        self.mFreeVirusGenomeTemplate = {}  # {('V1+','V2+',...'V1+'): id}
        self.mDeletedNormalSegList = {}  # {group1: [segH1, segH2]}
        self.mDeleteNormalCNByGroup = {}  # {group: CN_to_del}
        # self.readGraph2(input_fd)

    def readGraph(self, lines_read):
        source_t, sink_t, source_id, sink_id = [], [], [], []
        for line_no, line in enumerate(lines_read):
            if not line.strip() or line.strip().startswith('#'):
                continue
            sp = line.strip().split()
            sp[0] = sp[0].upper()

            if sp[0] == 'SAMPLE':
                self.mName = sp[1]
            elif sp[0] in ('AVG_DP', 'AVG_DEPTH'):
                self.mAvgCov = float(sp[1])  # might be overridden
                self.mAvgCovRaw = self.mAvgCov
            elif sp[0] == 'AVG_PLOIDY':
                self.mAvgPloidy = float(sp[1])
            elif sp[0] == 'GROUP_PRIORITY':
                self.mGroupPriority = sp[1].split(';')
            elif sp[0] == 'SOURCE':
                node = sp[1].split(':')
                source_t.append(self.getTypeByName(node[0]))
                source_id.append(int(node[1]))
            elif sp[0] == 'SINK':
                node = sp[1].split(':')
                sink_t.append(self.getTypeByName(node[0]))
                sink_id.append(int(node[1]))
            elif sp[0] in ('TUMOR_PURITY', 'PURITY'):
                self.mPurity = float(sp[1])
            elif sp[0] == 'PLOIDY':
                if 'm' in sp[1]:
                    ploidy = map(int, sp[1].split('m'))
                    ploidy[0] -= ploidy[1]
                else:
                    ploidy = map(int, sp[1].split(':'))
                grouping = None
                if len(sp) > 2:
                    tags = self.readTag(sp[2])
                    if 'GROUP' in tags:
                        grouping = tags['GROUP']
                self.mPloidy.update({grouping: ploidy})
                self.mPloidyRaw.update({grouping: ploidy})

            elif sp[0] in ('SEGMENT', 'SEG', 'NODE', 'S'):
                # SEGMENT | HOST:2 | 32.12 | WEIGHT=1.0;POS=chr1:100000-200000
                node = sp[1].split(':')
                seg_type = self.getTypeByName(node[0])
                seg_id = int(node[1])
                seg_cov = max(float(sp[2]), 0.0)
                seg_cred = 1.0
                seg_pos = 'None'
                seg_length = -1
                no_lower_bound = False
                grouping = None
                _override_cred_setting = False
                if len(sp) > 3:
                    tags = self.readTag(sp[3])
                    if 'WEIGHT' in tags:
                        # user can override this
                        seg_cred = float(tags['WEIGHT'])
                        _override_cred_setting = True
                    if 'POS' in tags:
                        seg_pos = tags['POS']
                        seg_length = self._checkIntervalFormat(seg_pos)
                    if 'NO_LOWER_BOUND' in tags:
                        no_lower_bound = True
                    if 'GROUP' in tags:
                        grouping = tags['GROUP']
                # further calibration by depth and length
                if not _override_cred_setting:
                    seg_cred1 = self._getCredByRange(SysSettings.WEIGHT_BY_LENGTH,
                                                     seg_length)
                    seg_cred2 = self._getCredByRange(SysSettings.WEIGHT_BY_CN,
                                                     seg_cov)
                    seg_cred = min(seg_cred1, seg_cred2)
                s = self.addSegment(seg_type, seg_id, seg_cov, cred=seg_cred, position=seg_pos)
                s.mLowerBoundLimit = not no_lower_bound
                s.mGroup = grouping
                if seg_type == 'V':  # cuz virus may be circular
                    self.mMaxVIdx = max(seg_id, self.mMaxVIdx)
                    self.mMinVIdx = min(seg_id, self.mMinVIdx)

            elif sp[0] in ('JUNC', 'JUN', 'JUNCTION', 'J', 'E', 'EDGE'):
                # JUNCTION | HOST:2:+ | VIRUS:4:- | 23.2 | WEIGHT=1.0;PREFERRED=1;
                vtx1 = sp[1].split(':')
                vtx2 = sp[2].split(':')
                vtx1_type, vtx2_type = map(self.getTypeByName, [vtx1[0], vtx2[0]])
                vtx1_id, vtx2_id = map(int, [vtx1[1], vtx2[1]])
                vtx1_dir, vtx2_dir = vtx1[2], vtx2[2]
                jun_cov = float(sp[3])
                jun_cred, preferred, no_lower_bound = 1.0, False, False
                _override_cred_setting = False
                if len(sp) > 4:
                    tags = self.readTag(sp[4])
                    if 'WEIGHT' in tags:
                        jun_cred = float(tags['WEIGHT'])
                        _override_cred_setting = True
                    if 'PREFERRED' in tags:
                        preferred = int(tags['PREFERRED']) != 0
                    if 'NO_LOWER_BOUND' in tags:
                        no_lower_bound = True
                if not _override_cred_setting:
                    jun_cred = self._getCredByRange(SysSettings.WEIGHT_BY_CN,
                                                    jun_cred)
                j = self.addJunction(vtx1_type, vtx1_id, vtx1_dir,
                                     vtx2_type, vtx2_id, vtx2_dir,
                                     jun_cov, cred=jun_cred, preferred=preferred)
                j.mLowerBoundLimit = not no_lower_bound

            elif sp[0] == 'HOST_WEIGHT':
                SysSettings.HOST_SEG_LP_COE = float(sp[1])
            elif sp[0] == 'VIRUS_WEIGHT':
                SysSettings.VIRUS_SEG_LP_COE = float(sp[1])
            elif sp[0] == 'JUNC_WEIGHT':
                SysSettings.JUNCTION_LP_COE = float(sp[1])
            elif sp[0] == 'INFER_COV_COE':
                SysSettings.INFER_COV_COE = float(sp[1])
            elif sp[0] in ('FVGM', 'FREE_VIRUS_GENOME'):
                fvgm_id = int(sp[1])
                fvgm = tuple(self.sortVirusGenome(sp[2]))  # tuple of string.
                self.mFreeVirusGenomeTemplate.update({fvgm: fvgm_id})

        for i, t in enumerate(source_t):
            source = self.getSegment(t, source_id[i])
            self.mSource.update({source.mGroup: source})
        for i, t in enumerate(sink_t):
            sink = self.getSegment(t, sink_id[i])
            self.mSink.update({sink.mGroup: sink})

        # if not enough info is given.
        if not self.mAvgPloidy:
            stdout('[INFO] Using default average ploidy setting: 2.')
            self.mAvgPloidy = 2
        if not self.mPloidy:
            stdout('[INFO] Using default ploidy setting: 1 for major allele, 1 for minor allele.')
            self.mPloidy = {None: [1, 1]}
        if len(self.mGroupPriority) != len(self.mPloidy):
            sys.exit('Number of groups declared in GROUP_PRIORITY should\n'
                     'be the same as the ones declared in SEGMENT field.\n'
                     'Group Priority: {}\nPloidy: {}'.format(self.mGroupPriority,
                                                             self.mPloidy))
        if not self.mPurity:
            stdout('Using default purity setting: 1.')
            self.mPurity = 1.0
        if not self.mAvgCov:
            stdout('Using default average depth setting: 2.')
            self.mAvgCov = 2
            self.mAvgCovRaw = 2

        self._adjustCovByPurity()

    def readGraph2(self, lines_read, verbose=False):
        source_t, sink_t, source_id, sink_id = [], [], [], []
        for line_no, line in enumerate(lines_read):
            if not line.strip() or line.strip().startswith('#'):
                continue
            sp = line.strip().split()
            sp[0] = sp[0].upper()

            if sp[0] == 'SAMPLE':
                self.mName = sp[1]
                if verbose: stdout('[INFO] Processing sample: {}'.format(self.mName))
            elif sp[0] in ('AVG_DP', 'AVG_DEPTH'):
                '''
                Average depth of a normal diploid segment. Used to calculate 
                copy number of segments.
                Optional, default is 2, stands for diploid.
                '''
                self.mAvgCov = float(sp[1])  # might be overridden
                self.mAvgCovRaw = self.mAvgCov
            elif sp[0] == 'AVG_PLOIDY':
                '''
                Optional, combining with avg_dp and purity, we calculate the
                expected haploid depth and update mAvgCov attribute. 
                Default is 2, when purity is 100%, this number keeps mAvgCov
                unchanged.
                '''
                self.mAvgPloidy = float(sp[1])
            elif sp[0] == 'GROUP_PRIORITY':
                # GROUP_PRIORITY\tGROUP1;GROUP2
                self.mGroupPriority = sp[1].split(';')
            elif sp[0] == 'SOURCE':
                _source_t, _source_id = self._readSeg(sp[1])
                source_t.append(_source_t)
                source_id.append(_source_id)
            elif sp[0] == 'SINK':
                _sink_t, _sink_id = self._readSeg(sp[1])
                sink_t.append(_sink_t)
                sink_id.append(_sink_id)
            elif sp[0] in ('TUMOR_PURITY', 'PURITY'):
                self.mPurity = float(sp[1])
            elif sp[0] == 'PLOIDY':
                '''
                Ploidy at source and sink. Using the format taken from Patchwork, 
                "AmB" means that the total copy number is A and minor allele CN 
                is B. Which means the CN at major allele is (A-B)
                
                PLOIDY\t3m1;GROUP=gp1
                '''
                tags = self._readTags(sp[1])
                _pl = tags[None]
                if 'm' in _pl:
                    ploidy = map(int, _pl.split('m'))
                    ploidy[0] -= ploidy[1]
                else:
                    ploidy = map(int, _pl.split(':'))
                grouping = tags['GROUP'] if 'GROUP' in tags else None
                self.mPloidy.update({grouping: ploidy})
                self.mPloidyRaw.update({grouping: ploidy})
            elif sp[0] in ('SEGMENT', 'SEG', 'NODE', 'S'):
                tags = self._readTags(sp[1])
                seg_type, seg_id = self._readSeg(tags['ID'])
                seg_cov = max(float(tags['CN_ORIG']), 0.0)
                seg_pos = tags['INTERVAL']
                seg_length = self._checkIntervalFormat(seg_pos)
                no_lower_bound = 'NO_LOWER_BOUND' in tags
                grouping = tags['GROUP'] if 'GROUP' in tags else None
                seg_cred = 1.0
                if 'W' in tags:
                    # user can override this with parameters
                    seg_cred = float(tags['W'])
                # further calibration by depth and length or other parameters.
                seg_cred1 = self._getCredByRange(SysSettings.WEIGHT_BY_LENGTH,
                                                 seg_length)
                seg_cred2 = self._getCredByRange(SysSettings.WEIGHT_BY_CN,
                                                 seg_cov)
                seg_cred = max(0.0, min(seg_cred, seg_cred1, seg_cred2))
                # Leave it to input file for final forced cred setting.
                seg_cred = max(0.0, float(tags['FORCE_W'])) if 'FORCE_W' in tags else seg_cred
                s = self.addSegment(seg_type, seg_id, seg_cov, cred=seg_cred, position=seg_pos)
                s.mLowerBoundLimit = not no_lower_bound
                s.mGroup = grouping
                if seg_type == 'V':  # cuz virus may be circular
                    self.mMaxVIdx = max(seg_id, self.mMaxVIdx)
                    self.mMinVIdx = min(seg_id, self.mMinVIdx)
                for _t in ('ID', 'CN_ORIG', 'GROUP'):
                    if _t in tags:
                        del tags[_t]
                s.mTags = tags
            elif sp[0] in ('JUNC', 'JUN', 'JUNCTION', 'J', 'E', 'EDGE'):
                tags = self._readTags(sp[1])
                jun_id = tags['ID']
                seg_link = tags['SEGLINK']
                vtx1, vtx2 = seg_link.split(',')
                vtx1_type, vtx1_id = self._readSeg(vtx1[:-1])
                vtx1_dir = vtx1[-1]
                vtx2_type, vtx2_id = self._readSeg(vtx2[:-1])
                vtx2_dir = vtx2[-1]
                jun_cov = float(tags['CN_ORIG'])
                jun_cred = 1.0
                if 'W' in tags:
                    # user can override this
                    jun_cred = float(tags['W'])
                preferred = 'PREFERRED' in tags
                no_lower_bound = 'NO_LOWER_BOUND' in tags
                # Leave it to input file for final forced cred setting.
                jun_cred = max(0.0, float(tags['FORCE_W'])) if 'FORCE_W' in tags else jun_cred
                j = self.addJunction(vtx1_type, vtx1_id, vtx1_dir,
                                     vtx2_type, vtx2_id, vtx2_dir,
                                     jun_cov, cred=jun_cred, preferred=preferred,
                                     aJuncId=jun_id)
                if 'FORCE_W' not in tags:
                    # Further "smart" settings of credibility.
                    jun_cred1 = self._getCredByRange(SysSettings.WEIGHT_BY_CN,
                                                     jun_cov)
                    jun_cred_seg1 = jun_cred_seg2 = 1
                    if j.mSource.getLength() != -1:
                        jun_cred_seg1 = self._getCredByRange(SysSettings.WEIGHT_BY_LENGTH,
                                                             j.mSource.getLength())
                    if j.mTarget.getLength() != -1:
                        jun_cred_seg2 = self._getCredByRange(SysSettings.WEIGHT_BY_LENGTH,
                                                             j.mTarget.getLength())
                    jun_cred_cap = SysSettings.CAP_JUNCTION_WEIGHT
                    jun_type_cap = SysSettings.SV_WEIGHT \
                        if j.mSource.mType == j.mTarget.mType \
                        else SysSettings.INTEG_WEIGHT
                    jun_cred = max(0.0, min(jun_cred,
                                            jun_cred1,
                                            jun_cred_seg1 * jun_cred_seg2,
                                            jun_cred_cap,
                                            jun_type_cap))
                    j.setCred(jun_cred)
                # j = self.addJunction(vtx1_type, vtx1_id, vtx1_dir,
                #                      vtx2_type, vtx2_id, vtx2_dir,
                #                      jun_cov, cred=jun_cred, preferred=preferred,
                #                      aJuncId=jun_id)

                j.mLowerBoundLimit = not no_lower_bound
                for _t in ('ID', 'SEGLINK', 'CN_ORIG'):
                    if _t in tags:
                        del tags[_t]
                j.mTags = tags
            elif sp[0] in ('FVGM', 'FREE_VIRUS_GENOME'):
                tags = self._readTags(sp[1])
                fvgm_id = tags['ID']
                fvgm = tuple(self.sortVirusGenome(tags['SEGSTR']))  # tuple of strings.
                self.mFreeVirusGenomeTemplate.update({fvgm: fvgm_id})
            elif sp[0] == 'ARGS':
                tags = self._readTags(sp[1])
                if 'H_LP_WEIGHT' in tags: SysSettings.HOST_SEG_LP_COE = float(tags['H_LP_WEIGHT'])
                if 'V_LP_WEIGHT' in tags: SysSettings.VIRUS_SEG_LP_COE = float(tags['V_LP_WEIGHT'])
                if 'J_LP_WEIGHT' in tags: SysSettings.JUNCTION_LP_COE = float(tags['J_LP_WEIGHT'])
                if 'INFER_COV_COE' in tags: SysSettings.INFER_COV_COE = float(tags['INFER_COV_COE'])
        for i, t in enumerate(source_t):
            source = self.getSegment(t, source_id[i])
            self.mSource.update({source.mGroup: source})
        for i, t in enumerate(sink_t):
            sink = self.getSegment(t, sink_id[i])
            self.mSink.update({sink.mGroup: sink})

        # if not enough info is given.
        if not self.mAvgPloidy:
            if verbose: stdout('[INFO] Using default average ploidy setting: 2.')
            self.mAvgPloidy = 2
        if not self.mPloidy:
            if verbose: stdout('[INFO] Using default ploidy setting: 1 for major allele, 1 for minor allele.')
            self.mPloidy = {None: [1, 1]}
        if len(self.mGroupPriority) != len(self.mPloidy):
            sys.exit('Number of groups declared in GROUP_PRIORITY should\n'
                     'be the same as the ones declared in SEGMENT field.\n'
                     'Group Priority: {}\nPloidy: {}'.format(self.mGroupPriority,
                                                             self.mPloidy))
        if not self.mPurity:
            if verbose: stdout('[INFO] Using default purity setting: 1.')
            self.mPurity = 1.0
        if not self.mAvgCov:
            if verbose: stdout('[INFO] Using default average depth setting: 2.')
            self.mAvgCov = 2.0
            self.mAvgCovRaw = 2.0

        self._adjustCovByPurity()

    def _readSeg(self, seg_str):
        seg_type, seg_idx = seg_str[0], int(seg_str[1:])
        if seg_type not in ('H', 'V'):
            sys.exit('Wrong segment type: {}'.format(seg_type))
        return seg_type, seg_idx

    @staticmethod
    def _readTags(tag_str):
        tags = OrderedDict()
        tags_sp = tag_str.split(';')
        for tag in tags_sp:
            if not tag.strip():
                continue
            sp = tag.split('=')
            if len(sp) == 1:
                tags.update({None: sp[0]})
            else:
                tags.update({sp[0]: sp[1]})
        return tags

    def readTag(self, tag_str):
        output_tags = {}
        tags = {}
        [tags.update({kv.split('=')[0]: '='.join(kv.split('=')[1:])}) for kv in tag_str.strip().split(';')]
        if 'WEIGHT' in tags:
            output_tags.update({'WEIGHT': float(tags['WEIGHT'])})
        if 'W' in tags:
            output_tags.update({'WEIGHT': float(tags['W'])})
        if 'WT' in tags:
            output_tags.update({'WEIGHT': float(tags['WT'])})
        if 'PREFERRED' in tags:
            output_tags.update({'PREFERRED': float(tags['PREFERRED'])})
        if 'PREF' in tags:
            output_tags.update({'PREFERRED': float(tags['PREF'])})
        if 'POS' in tags:
            output_tags.update({'POS': tags['POS']})
        if 'INTERVAL' in tags:
            output_tags.update({'POS': tags['INTERVAL']})
        if 'NOLB' in tags:
            output_tags.update({'NO_LOWER_BOUND': None})
        if 'GP' in tags:
            output_tags.update({'GROUP': tags['GP']})
        if 'GROUP' in tags:
            output_tags.update({'GROUP': tags['GROUP']})
        return output_tags

    def _checkIntervalFormat(self, interval_str):
        try:
            chrom, interval = interval_str.split(':')
            st, ed = map(int, interval.split('-'))
            return ed - st + 1
        except:
            sys.exit('Wrong interval format: {}'.format(interval_str))

    def _getCredByRange(self, range_weight_list, value):
        '''

        format: [(0, 100, 0.1), (101, 1000, 0.5), ...]
        '''
        if not range_weight_list:
            return 1.0
        for st, ed, cred in range_weight_list:
            if st <= value <= ed:
                return cred
        # If the range is not declared in the list. use smallest one.
        return range_weight_list[-1][2]

    def sortVirusGenome(self, virus_genome_str):
        virus_vertices = virus_genome_str.split(',')
        # In V1+,V2+,V3+,V1+, last V1+ is not needed.
        if virus_vertices[0] != virus_vertices[-1]:
            virus_vertices.append(virus_vertices[0])
        return self._reorderCycle(virus_vertices, by_virus=True)

    def _reorderCycle(self, cycle_str_list, by_virus=True):
        """e.g. Turn 5,3,2,4,5 to 2,4,5,3,2 """
        min_list_idx = -1
        min_vertex_idx = float('inf')
        for i, v_str in enumerate(cycle_str_list[:-1]):
            if by_virus and v_str[0] == 'V':
                if int(v_str[1:-1]) < min_vertex_idx:
                    min_vertex_idx = int(v_str[1:-1])
                    min_list_idx = i
            if not by_virus and v_str[0] == 'H':
                if int(v_str[1:-1]) < min_vertex_idx:
                    min_vertex_idx = int(v_str[1:-1])
                    min_list_idx = i
        return cycle_str_list[min_list_idx:-1] + \
               cycle_str_list[:min_list_idx] + \
               [cycle_str_list[min_list_idx]]

    def _adjustCovByPurity(self):
        """Adjust coverage according to purity information before
        adding missing edges.

        > We assume that the Normal sample is diploid, Given tumor purity [Pu],
          average tumor ploidy [AvgTumorPloidy], whole genome average depth [AvgDepth]
          we want the average sequencing depth for one copy of genome [Dh] (i.e. depth for haploid).

          Obviously, Dh = AvgDepth / AvgPloidy

          The meanCN column in Patchwork results stands for the copy number at the mean coverage. It
          describes the mean copy number in tumor sample.

          Therefore, AvgPloidy = Pu * AvgTumorPloidy + (1 - Pu) * 2

        > Then the copy number of a certain segment from tumor genome given [SegDepth], is computed as:

          Pu * SegCopy * Dh + (1 - Pu) * 2 * Dh = SegDepth

          Which means that: SegCopy = (SegDepth - (1 - Pu) * 2 * Dh ) / (Pu * Dh)

        Note that we take AvgPloidy and Pu from Patchwork.
        If you have done this sort of normalizations manually, you should set AVG_DEPTH and AVG_PLOIDY
        in your input as 2.
        """
        haploid_depth = self.mAvgCovRaw / ((self.mAvgPloidy - 2) * self.mPurity + 2)
        self.mAvgCov = 2 * haploid_depth
        for s in self.mSegments:
            if s.mType == 'H':
                # might be a negative number.
                # s.mCov.mCovAdjusted = max((s.mCov.mCov - (1 - self.mPurity) * self.mAvgCov) / self.mPurity, 0.0)
                # s.mCov.mCov = s.mCov.mCovAdjusted
                cov_adjusted = max((s.mCov.mCov - (1 - self.mPurity) * self.mAvgCov) / self.mPurity, 0.0)
                s.setCov(cov_adjusted)

    def setLoop(self):
        """Set a connection from sink+ => source+.

            This must be run after reachability has been ensured given the
            way the code for reachability is written (if not, there will always
            be a way to reach sink and source by traversing in the same
            direction.)
        """
        for gp in self.mSource:
            source = self.mSource[gp]
            sink = self.mSink[gp]

            self.addJunction(sink.mPosV.mType, sink.mPosV.mIdx, '+',
                             source.mPosV.mType, source.mPosV.mIdx, '+',
                             (source.mCov.mCov + sink.mCov.mCov) / 2,
                             source_to_sink=True, cred=0)

    def calculateWeight(self, drop=False, drop_setting=None):
        self.backupCov()
        if not drop:
            for s in self.mSegments:
                s.setWeight(2.0 * s.getCov() / self.mAvgCov)
            for j in self.mJunctions:
                j.setWeight(2.0 * j.getCov() / self.mAvgCov)
            return 1
        else:
            remove_copy = {}
            for gp in self.mPloidyRaw:
                self.mDeletedNormalSegList.update({gp: []})
                if drop_setting[gp] == 'M':
                    self.mPloidy.update({gp: (self.mPloidyRaw[gp][1], 0)})
                    remove_copy.update({gp: self.mPloidyRaw[gp][0]})
                elif drop_setting[gp] in ('m', '-'):
                    self.mPloidy.update({gp: (self.mPloidyRaw[gp][0], 0)})
                    remove_copy.update({gp: self.mPloidyRaw[gp][1]})
            '''
            First check if all groups can be converted to 1m0 safely.
            Or, as a compromise, find the GCD among groups' copy number.
            '''
            if SysSettings.LP_WITH_ORIGINAL_PLOIDY:
                groups_cn_gcd = 1
            else:
                groups_cn_gcd = reduce(gcd,
                                       [_[0] for _ in self.mPloidy.values()])
            # Update mPloidy by dividing the gcd
            for gp in self.mPloidy:
                self.mPloidy[gp] = (self.mPloidy[gp][0] / groups_cn_gcd, 0)
            self.mCNDividedBy = groups_cn_gcd
            # Set seg and junction weight
            for s in self.mSegments:
                if s.mType == 'H':
                    s.setCov(max(0, s.getCov() - remove_copy[s.getGroup()] * self.mAvgCov / 2))
                    self.mDeletedNormalSegList[s.getGroup()].append(s)
                s.setWeight(2.0 * s.getCov() / (self.mAvgCov * groups_cn_gcd))
            for j in self.mJunctions:
                if self.isNormalHumanJunction(j):
                    j.setCov(max(0, j.getCov() - remove_copy[j.mSource.getGroup()] * self.mAvgCov / 2))
                j.setWeight(2.0 * j.getCov() / (self.mAvgCov * groups_cn_gcd))
            self.mDeleteNormalCNByGroup = remove_copy
            for gp in self.mDeletedNormalSegList:
                self.mDeletedNormalSegList[gp] = \
                    sorted(self.mDeletedNormalSegList[gp], key=lambda x: x.mIdx)
                # Ensure it is a cycle, ending with the last segment.
                self.mDeletedNormalSegList[gp].append(self.mDeletedNormalSegList[gp][0])
            return groups_cn_gcd

    @staticmethod
    def getTypeByName(name):
        if name.lower() in ('virus', 'vir', 'v'):
            return 'V'
        if name.lower() in ('host', 'human', 'h'):
            return 'H'
        sys.exit('invalid segment type ' + name)

    def addSegment(self, aType, aIdx, aCov, cred=1.0, position='None'):
        s = Segment(aType, aIdx, aCov, cred, position)
        self.mSegments.append(s)
        return s

    def getSegment(self, aType, aIdx):
        for s in self.mSegments:
            if s.mType == aType and s.mIdx == aIdx:
                return s

    def getStartVertex(self):
        """20170907 update."""
        for gp in self.mGroupPriority:
            source_v = self.mSource[gp].mPosV
            if source_v.hasWeight():
                return source_v
        return False

    def getGroupStartVertex(self, group):
        return self.mSource[group].mPosV

    def addJunction(self,
                    aSourceType, aSourceIdx, aSourceDir,
                    aTargetType, aTargetIdx, aTargetDir,
                    aCov, cred=1.0, inferred=False, preferred=False,
                    source_to_sink=False, aJuncId=None):
        source = self.getSegment(aSourceType, aSourceIdx)
        target = self.getSegment(aTargetType, aTargetIdx)
        if not source or not target:
            sys.exit('Error while adding junction {} {} => {} {}. '
                     'Make sure the nodes exist.'.format(aSourceType,
                                                         aSourceIdx,
                                                         aTargetType,
                                                         aTargetIdx))
        j = Junction(source, aSourceDir,
                     target, aTargetDir,
                     aCov, cred, inferred, preferred,
                     aSourceSinkJunction=source_to_sink,
                     aJuncId=aJuncId)
        if self.exists(j.mIdStr):
            raise DuplicateJunctionException('{0} exists!'.format(j))
        j.addJunctionToVertices()
        self.mJunctions.append(j)
        if inferred:
            self.mLog['inferred_junctions'].append(j)
        stderr('[INFO] Added junction: {}, cred: {}'.format(j, j.mCredibility))
        return j

    def delJunctions(self, *args):
        for j in args:
            stderr('[INFO] Removed junction: ', j)
            j.selfDestroy()
            self.mJunctions.remove(j)

    def isNormalHumanJunction(self, junction):
        e = junction.mEdgeA
        if e.mSourceV.mType == e.mTargetV.mType == 'H':
            if e.connectsSameDir():
                if e.mSourceV.mDir == '+' and e.mSourceV.mIdx == e.mTargetV.mIdx - 1:
                    return True
                if e.mSourceV.mDir == '-' and e.mSourceV.mIdx == e.mTargetV.mIdx + 1:
                    return True
                if e.isImaginaryJunction():
                    return True
        return False

    def isSource(self, seg):
        try:
            return seg == self.mSource[seg.getGroup()]
        except:
            # in case that virus segments have None as group and cause KeyError.
            return False

    def isSink(self, seg):
        try:
            return seg == self.mSink[seg.getGroup()]
        except:
            return False

    def getJunctionIdx(self, edge):
        for idx, j in enumerate(self.mJunctions):
            if edge == j.mEdgeA or edge == j.mEdgeB:
                return idx
        return None

    def exists(self, junc_idstr_list):
        for j in self.mJunctions:
            if junc_idstr_list[0] in j.mIdStr and junc_idstr_list[1] in j.mIdStr:
                return True
        return False

    def setOrphan(self):
        """Called right after the creation of graph.

        This only shows whether the segment is orphan at the very BEGINNING.
        """
        for s in self.mSegments:
            if s.isOrphan():
                s.setAsOrphan()

    def resetVisitedFlag(self):
        for j in self.mJunctions:
            j.mEdgeA.mVisited = False
            j.mEdgeB.mVisited = False

    def setInferredJunctionCred(self, heuristic=False):
        """Deprecated. Called right before balancing."""
        # if heuristic:
        #     # Disabled, makes no sense.
        #     max_cred = 1
        #     for j in self.mJunctions:
        #         if j.mInferred:
        #             max_cred = max(max_cred, j.mCredibility ** 0.5)
        #     for j in self.mJunctions:
        #         if j.mInferred:
        #             j.mCredibility = j.mCredibility ** 0.5 / max_cred
        #             stderr('Inferred junction {}\tCred: {}'
        #                    .format(j, j.mCredibility))
        #     stderr()
        # else:
        #     for j in self.mJunctions:
        #         if j.mInferred:
        #             print SysSettings.INFERRED_JUNC_COE
        #             j.mCredibility = j.mCredibility * SysSettings.INFERRED_JUNC_COE
        #             stderr('Inferred junction {}\tCred: {}'
        #                    .format(j, j.mCredibility))
        pass

    def backupCov(self):
        for j in self.mJunctions:
            j.mCov.backupCov()
        for s in self.mSegments:
            s.mCov.backupCov()

    def backupWeight(self):
        for j in self.mJunctions:
            j.mCov.backup()
        for s in self.mSegments:
            s.mCov.backup()

    def refreshGraph(self):
        for j in self.mJunctions:
            j.mCov.restore()
        for s in self.mSegments:
            s.mCov.restore()
            s.mPosV.mLastEdgeChoice = None
            s.mNegV.mLastEdgeChoice = None

    def backupWeight2(self):
        j_w, s_w, v_last_choice = [], [], []
        for j in self.mJunctions:
            j_w.append(j.getWeight())
        for s in self.mSegments:
            s_w.append(s.getWeight())
            v_last_choice.append((s.mPosV.mLastEdgeChoice, s.mNegV.mLastEdgeChoice))
        return j_w, s_w, v_last_choice

    def restoreWeight2(self, j_w, s_w, v_last_choice):
        for i, j in enumerate(self.mJunctions):
            j.setWeight(j_w[i])
        for i, s in enumerate(self.mSegments):
            s.setWeight(s_w[i])
            s.mPosV.mLastEdgeChoice = v_last_choice[i][0]
            s.mNegV.mLastEdgeChoice = v_last_choice[i][1]

    def forkingVertices(self):
        V = []
        N = []
        for s in self.mSegments:
            if len(s.mPosV.mNextEdges) > 1:
                V.append(s.mPosV)
                N.append(len(s.mPosV.mNextEdges))
            if len(s.mNegV.mNextEdges) > 1:
                V.append(s.mNegV)
                N.append(len(s.mNegV.mNextEdges))
        return V, N

    def inSameGroup(self, vtx1, vtx2):
        return vtx1.mSeg.mGroup == vtx2.mSeg.mGroup

    def edgeReachesImaginaryJunction(self, edge):
        """Tell if edge connects to H1-"""
        if edge.mTargetV.mType == 'H':
            return edge.mTargetV == self.mSource[edge.mTargetV.mSeg.getGroup()].mNegV
        return False

    def isNormalAllele(self, cyc):
        # Do not pass H1+;H2+;H3+;H1+, pass H1+;H2+;H3+ instead
        group = cyc[0].mSeg.getGroup()
        if cyc[0] != self.mSource[group].mPosV or cyc[-1] != self.mSink[group].mPosV:
            return False
        if cyc[0].mType != 'H':
            return False
        start_idx = cyc[0].mIdx
        current_idx = start_idx
        for i, v in enumerate(cyc[1:]):
            if current_idx + 1 != v.mIdx or v.mType != 'H' or v.mDir != '+':
                return False
            current_idx += 1
        return True

    def getRawPloidy(self):
        return self.mPloidyRaw

    def printGraphInfo(self):
        stdout('------------- GRAPH INFO -------------')
        stdout('AVG: {0:.2f}'.format(self.mAvgCov))
        for s in self.mSegments:
            stdout('{0}{1}, cov: {2:.2f}, cov(original): {4:.2f} weight: {3:.3f}' \
                   .format(s.mType, s.mIdx, s.mCov.mCov, s.mCov.mWeight, s.mCov.mCovOriginal))
            stdout('\tPosV ({0}) nextEdge:\t'.format(s.mPosV.getAbbr()), end='')
            for e in s.mPosV.mNextEdges:
                stdout(e.getAbbr() + ':' + str(e.getWeight()) + '\t', end='')
            stdout('\n\tPosV ({0}) prevEdge:\t'.format(s.mPosV.getAbbr()), end='')
            for e in s.mPosV.mPrevEdges:
                stdout(e.getAbbr() + ':' + str(e.getWeight()) + '\t', end='')
            stdout('\n\tNegV ({0}) nextEdge:\t'.format(s.mNegV.getAbbr()), end='')
            for e in s.mNegV.mNextEdges:
                stdout(e.getAbbr() + ':' + str(e.getWeight()) + '\t', end='')
            stdout('\n\tNegV ({0}) prevEdge:\t'.format(s.mNegV.getAbbr()), end='')
            for e in s.mNegV.mPrevEdges:
                stdout(e.getAbbr() + ':' + str(e.getWeight()) + '\t', end='')
            stdout()
        stdout('--------------------------------------')

from collections import namedtuple
# from Util.Utilities import TraversalException


class Vertex:
    def __init__(self, aType, aIdx, aDir, aCovObj, aCred):
        self.mType = aType
        self.mIdx = aIdx
        self.mDir = aDir
        self.mCov = aCovObj
        self.mCredibility = aCred
        self.mNextEdges = []
        self.mPrevEdges = []
        self.mIsOrphan = False
        self.mLastEdgeChoice = None  # previous choice of next edge
        self.mSeg = None
        self.mSafeNextEdges = []

    def hasWeight(self):
        return self.mCov.mWeight >= 1

    def getWeight(self):
        return self.mCov.mWeight

    def setWeight(self, w):
        self.mCov.mWeight = w

    def getCov(self):
        return self.mCov.mCov

    def getOutCov(self):
        return sum([e.getCov() if not e.isPalindrome() else e.getCov() * 2
                    for e in self.mNextEdges])

    def getInCov(self):
        return sum([e.getCov() if not e.isPalindrome() else e.getCov() * 2
                    for e in self.mPrevEdges])

    def getAbbr(self):
        return '{0}{1}{2}'.format(self.mType, self.mIdx, self.mDir)

    def getComplementAbbr(self):
        _d = '-' if self.mDir == '+' else '+'
        return '{0}{1}{2}'.format(self.mType, self.mIdx, _d)

    def getGroup(self):
        return self.mSeg.getGroup()

    def getSeg(self):
        return self.mSeg

    def getSiblingVtx(self):
        return self.mSeg.mPosV if self.mDir == '-' else self.mSeg.mNegV

    def __str__(self):
        return '({0}{1}{2}, {3:.2f}, {4:.3f})'.format(self.mType,
                                                      self.mIdx,
                                                      self.mDir,
                                                      self.getCov(),
                                                      self.getWeight())

    def __repr__(self):
        return self.getAbbr()


class Segment:
    def __init__(self, aType, aIdx, aCov, aCred, position=None):
        self.mType = aType
        self.mIdx = aIdx
        self.mCov = Coverage(aCov)
        self.mPosV = Vertex(aType, aIdx, '+', self.mCov, aCred)
        self.mNegV = Vertex(aType, aIdx, '-', self.mCov, aCred)
        self.mPosV.mSeg = self
        self.mNegV.mSeg = self
        self.mCredibility = aCred
        self.mIsOrphan = False  # if originally is orphan
        self.mLowerBoundLimit = True
        self.mPosition = position
        self.mGroup = None
        self.mTags = {}
        self.lp_stat_template = namedtuple('LpStat',
                                           ['orig_cn', 'lp_cn',
                                            'lp_cred', 'diff_pct'])
        # namedtuple, for reporting Lp result in the output.
        self.mLpStat = self.lp_stat_template(self.getWeight(),
                                             'N/A',
                                             'N/A',
                                             'N/A')

    def setCov(self, c):
        self.mCov.mCov = c

    def getCov(self):
        return self.mCov.mCov

    def getGroup(self):
        return self.mGroup

    def getWeight(self):
        return self.mCov.mWeight

    def setWeight(self, w):
        self.mCov.mWeight = w

    def getAbbr(self):
        return '{0}{1}'.format(self.mType, self.mIdx)

    def isDeadEnd(self):
        if self.mPosV.mNextEdges and self.mPosV.mPrevEdges:
            return False
        return True

    def isOrphan(self):
        if self.mPosV.mNextEdges or self.mPosV.mPrevEdges:
            return False
        return True

    def setAsOrphan(self):
        self.mIsOrphan = True
        self.mPosV.mIsOrphan = True
        self.mNegV.mIsOrphan = True

    def setCred(self, cred):
        self.mCredibility = cred
        self.mPosV.mCredibility = cred
        self.mNegV.mCredibility = cred

    def setLpStat(self, newCn, cred, diff_pct):
        self.mLpStat = self.lp_stat_template(self.getWeight(),
                                             newCn,
                                             cred,
                                             diff_pct)
        # self.mLpStat = 'LP={:.2f}:{:.2f}:{:.2f}'.format(self.getWeight(), newCn, cred)

    def getLpStat(self):
        return self.mLpStat

    def getLength(self):
        if self.mPosition:
            st, ed = map(int, self.mPosition.split(':')[1].split('-'))
            return ed - st + 1
        else:
            return -1

    def __str__(self):
        return '({0}{1}, {2:.2f}, {3:.3f})'.format(self.mType,
                                                   self.mIdx,
                                                   self.getCov(),
                                                   self.getWeight())

    def __repr__(self):
        return self.getAbbr()


class Edge:
    def __init__(self, aSourceV, aTargetV, aCovObj, aCred):
        self.mSourceV = aSourceV
        self.mTargetV = aTargetV
        self.mCov = aCovObj
        self.mVisited = False
        # self.mCredibility = aCred 20171020 Not actually used.
        self.mPreferred = False
        self.mIsPalindrome = False

    def hasWeight(self):
        return self.mCov.mWeight > 0

    def getWeight(self):
        return self.mCov.mWeight

    def setWeight(self, w):
        self.mCov.mWeight = w

    def setCov(self, cov):
        self.mCov.mCov = cov

    def getCov(self):
        return self.mCov.mCov

    def traverse(self, cn=1):
        self.mCov.mWeight -= cn
        self.mSourceV.mCov.mWeight -= cn
        return self.mCov.mWeight

    def getAbbr(self):
        return '[{0} => {1}]'.format(self.mSourceV.getAbbr(),
                                     self.mTargetV.getAbbr())

    def connectsSameDir(self):
        return self.mSourceV.mDir == self.mTargetV.mDir

    def isImaginaryJunction(self):
        return self.mCov.mImaginary

    def isPalindrome(self):
        return self.mIsPalindrome

    def __str__(self):
        return '[{0} => {1}, {4:.2f}->{2:.2f}, {3:.2f}]' \
            .format(self.mSourceV, self.mTargetV, self.mCov.mCov,
                    self.mCov.mWeight, self.mCov.mCovOriginal)


class Junction:
    def __init__(self, aSource, aSourceDir, aTarget, aTargetDir, aCov, aCred,
                 aInferred=False, aPreferred=False, aSourceSinkJunction=False, aJuncId=None):
        self.JUNCID = aJuncId
        self.mSource = aSource
        self.mTarget = aTarget
        self.mSourceDir = aSourceDir
        self.mTargetDir = aTargetDir
        self.mCov = Coverage(aCov)
        self.mCov.mImaginary = aSourceSinkJunction
        self.mCredibility = aCred
        self.mLowerBoundLimit = True
        self.mTags = {}
        self.mIsPalindrome = False
        # self.mPreferred = aPreferred

        if aSourceDir == '+' and aTargetDir == '+':
            self.mEdgeA = Edge(aSource.mPosV, aTarget.mPosV, self.mCov, aCred)
            self.mEdgeB = Edge(aTarget.mNegV, aSource.mNegV, self.mCov, aCred)
        elif aSourceDir == '-' and aTargetDir == '-':
            self.mEdgeA = Edge(aSource.mNegV, aTarget.mNegV, self.mCov, aCred)
            self.mEdgeB = Edge(aTarget.mPosV, aSource.mPosV, self.mCov, aCred)
        elif aSourceDir == '+' and aTargetDir == '-':
            self.mEdgeA = Edge(aSource.mPosV, aTarget.mNegV, self.mCov, aCred)
            self.mEdgeB = Edge(aTarget.mPosV, aSource.mNegV, self.mCov, aCred)
        elif aSourceDir == '-' and aTargetDir == '+':
            self.mEdgeA = Edge(aSource.mNegV, aTarget.mPosV, self.mCov, aCred)
            self.mEdgeB = Edge(aTarget.mNegV, aSource.mPosV, self.mCov, aCred)

        if self.mEdgeA.getAbbr() == self.mEdgeB.getAbbr():
            # When junction is a "palindrome" e.g. H2+ -> H2- or H2- -> H2+
            self.setAsPalindrome()

        self.mIdStr = [self.mEdgeA.getAbbr(), self.mEdgeB.getAbbr()]
        self.mInferred = aInferred
        self.mEdgeA.mPreferred = aPreferred
        self.mEdgeB.mPreferred = aPreferred

        self.lp_stat_template = namedtuple('LpStat',
                                           ['orig_cn', 'lp_cn',
                                            'lp_cred', 'diff_pct'])
        # namedtuple, for reporting Lp result in the output.
        self.mLpStat = self.lp_stat_template(self.getWeight(),
                                             'N/A',
                                             'N/A',
                                             'N/A')


    def addJunctionToVertices(self):
        if self.mEdgeA == self.mEdgeB:
            # When junction is a "palindrome" e.g. H2+ -> H2- or H2- -> H2+
            if self.mSourceDir == '+' and self.mTargetDir == '-':
                self.mSource.mPosV.mNextEdges.append(self.mEdgeA)
                self.mSource.mNegV.mPrevEdges.append(self.mEdgeB)
            elif self.mSourceDir == '-' and self.mTargetDir == '+':
                self.mSource.mNegV.mNextEdges.append(self.mEdgeA)
                self.mSource.mPosV.mPrevEdges.append(self.mEdgeB)
        else:
            if self.mSourceDir == '+' and self.mTargetDir == '+':
                self.mSource.mPosV.mNextEdges.append(self.mEdgeA)
                self.mSource.mNegV.mPrevEdges.append(self.mEdgeB)
                self.mTarget.mNegV.mNextEdges.append(self.mEdgeB)
                self.mTarget.mPosV.mPrevEdges.append(self.mEdgeA)
            elif self.mSourceDir == '-' and self.mTargetDir == '-':
                self.mSource.mNegV.mNextEdges.append(self.mEdgeA)
                self.mSource.mPosV.mPrevEdges.append(self.mEdgeB)
                self.mTarget.mPosV.mNextEdges.append(self.mEdgeB)
                self.mTarget.mNegV.mPrevEdges.append(self.mEdgeA)
            elif self.mSourceDir == '+' and self.mTargetDir == '-':
                self.mSource.mPosV.mNextEdges.append(self.mEdgeA)
                self.mSource.mNegV.mPrevEdges.append(self.mEdgeB)
                self.mTarget.mPosV.mNextEdges.append(self.mEdgeB)
                self.mTarget.mNegV.mPrevEdges.append(self.mEdgeA)
            elif self.mSourceDir == '-' and self.mTargetDir == '+':
                self.mSource.mNegV.mNextEdges.append(self.mEdgeA)
                self.mSource.mPosV.mPrevEdges.append(self.mEdgeB)
                self.mTarget.mNegV.mNextEdges.append(self.mEdgeB)
                self.mTarget.mPosV.mPrevEdges.append(self.mEdgeA)

    def setAsPalindrome(self):
        self.mEdgeB = self.mEdgeA
        self.mEdgeA.mIsPalindrome = True
        self.mIsPalindrome = True

    def isPalindrome(self):
        return self.mIsPalindrome

    def setWeight(self, w):
        self.mCov.mWeight = w

    def getWeight(self):
        return self.mCov.mWeight

    def setCov(self, cov):
        self.mCov.mCov = cov

    def getCov(self):
        return self.mCov.mCov

    def setCred(self, cred):
        self.mCredibility = cred
        self.mEdgeA.mCredibility = cred
        self.mEdgeB.mCredibility = cred

    def getVariableName(self):
        sd = td = 'p'
        if self.mSourceDir == '-':
            sd = 'n'
        if self.mTargetDir == '-':
            td = 'n'
        return '{0}_{1}'.format(self.mSource.getAbbr() + sd,
                                self.mTarget.getAbbr() + td)

    def selfDestroy(self):
        if self.mSourceDir == '+' and self.mTargetDir == '+':
            self.mSource.mPosV.mNextEdges.remove(self.mEdgeA)
            self.mSource.mNegV.mPrevEdges.remove(self.mEdgeB)
            self.mTarget.mNegV.mNextEdges.remove(self.mEdgeB)
            self.mTarget.mPosV.mPrevEdges.remove(self.mEdgeA)
        elif self.mSourceDir == '-' and self.mTargetDir == '-':
            self.mSource.mNegV.mNextEdges.remove(self.mEdgeA)
            self.mSource.mPosV.mPrevEdges.remove(self.mEdgeB)
            self.mTarget.mPosV.mNextEdges.remove(self.mEdgeB)
            self.mTarget.mNegV.mPrevEdges.remove(self.mEdgeA)
        elif self.mSourceDir == '+' and self.mTargetDir == '-':
            self.mSource.mPosV.mNextEdges.remove(self.mEdgeA)
            self.mSource.mNegV.mPrevEdges.remove(self.mEdgeB)
            self.mTarget.mPosV.mNextEdges.remove(self.mEdgeB)
            self.mTarget.mNegV.mPrevEdges.remove(self.mEdgeA)
        elif self.mSourceDir == '-' and self.mTargetDir == '+':
            self.mSource.mNegV.mNextEdges.remove(self.mEdgeA)
            self.mSource.mPosV.mPrevEdges.remove(self.mEdgeB)
            self.mTarget.mNegV.mNextEdges.remove(self.mEdgeB)
            self.mTarget.mPosV.mPrevEdges.remove(self.mEdgeA)

    def setLpStat(self, newCn, cred, diff_pct):
        self.mLpStat = self.lp_stat_template(self.getWeight(),
                                             newCn,
                                             cred,
                                             diff_pct)

    def getLpStat(self):
        return self.mLpStat

    def isImaginaryJunction(self):
        return self.mCov.mImaginary

    def __str__(self):
        return str('{ ' + self.mIdStr[0] + ' / ' + self.mIdStr[1] + ' }')


class Coverage:
    def __init__(self, aCov):
        self.mCov = aCov  # Before calculating weight, not used after weight is calculated.
        self.mCovOriginal = aCov  # original input
        # self.mCovAdjusted = aCov  # as adjusted by purity
        self.mCovBackup = aCov  # cov before dropping normal copies.
        self.mWeight = 0
        self.mWeightBackup = 0
        self.mImaginary = False

    def backup(self):
        self.mWeightBackup = self.mWeight

    def restore(self):
        self.mWeight = self.mWeightBackup

    def backupCov(self):
        self.mCovBackup = self.mCov
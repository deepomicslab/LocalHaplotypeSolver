import sys
import ConfigParser
import os
import SysSettings


def stderr(*msgs, **kwargs):
    msgs = [str(i) for i in msgs]
    end = kwargs['end'] if 'end' in kwargs else '\n'
    if SysSettings.VERBOSE:
        if 'sep' in kwargs:
            sys.stderr.write(kwargs['sep'].join(msgs) + end)
        else:
            sys.stderr.write(' '.join(msgs) + end)


def stdout(*msgs, **kwargs):
    msgs = [str(i) for i in msgs]
    end = kwargs['end'] if 'end' in kwargs else '\n'
    color = kwargs.get('color', None)
    if color:
        sys.stdout.write(color)
    if 'sep' in kwargs:
        sys.stdout.write(kwargs['sep'].join(msgs) + end)
    else:
        sys.stdout.write(' '.join(msgs) + end)
    if color:
        sys.stdout.write(StrColors.ENDC)


def readConfig(path=os.path.join(os.path.dirname(__file__), '../config.ini')):
    conf = ConfigParser.ConfigParser()
    conf.read(path)
    SysSettings.VERBOSE = conf.get('logging', 'verbose').lower() == 'true'
    SysSettings.EDGE_SELECTION = conf.get('edge_selection',
                                          'selection_mode')
    SysSettings.DELETE_INCOMPATIBLE_EDGES = \
        conf.get('edge_selection',
                 'exclude_incompatible_edges').lower() == 'true'
    SysSettings.HOST_SEG_LP_COE = float(conf.get('linear_programming_weight',
                                                 'host_segment_overall'))
    SysSettings.VIRUS_SEG_LP_COE = float(conf.get('linear_programming_weight',
                                                  'virus_segment_overall'))
    SysSettings.JUNCTION_LP_COE = float(conf.get('linear_programming_weight',
                                                 'junction_overall'))


class DuplicateJunctionException(Exception):
    def __init__(self, message):
        super(DuplicateJunctionException, self).__init__(message)


class NotReachableException(Exception):
    pass


class ImaginaryJunctionReachedException(Exception):
    pass


class BreakPointNotEnumeratedException(Exception):
    breakPoint = None


class NotInSameReferenceGroupException(Exception):
    pass


class MixedGroupException(Exception):
    """When a cycle has human segments from more than one group."""
    pass


class TraversalException(Exception):
    pass


class TooManyTraversalException(Exception):
    pass


def printLGM(lgm, color=None):
    """This already handles print function.
    If use print printLGM(), there will be a "None" being printed.
    """
    for cyc in lgm:
        stdout('<', *[c.getAbbr() for c in cyc], end=' >\n', color=color)

def cycleToStr(cyc):
    return ','.join([c.getAbbr() for c in cyc])

def print_cnt_cyc_dict(d):
    # stdout('{')
    for i in sorted(d.keys()):
        if not d[i]:
            continue
        stdout('\t', i, ': ', end='')
        stdout('\t<', *[c.getAbbr() for c in d[i][0]], end=' >\n')
        for cyc in d[i][1:]:
            stdout('\t\t<', *[c.getAbbr() for c in cyc], end=' >\n')
    # stdout('}')


class StrColors:
    def __init__(self):
        pass
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
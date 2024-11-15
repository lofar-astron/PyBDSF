""" WARNING, ERROR, and CRITICAL are always output to screen and to log file.
INFO and USERINFO always go to the log file. DEBUG goes to log file if debug is
True. USERINFO goes to screen only if quiet is False.

Use as follows:

mylog = mylogger.logging.getLogger("name")

mylog.info('info') --> print to logfile, but not to screen
mylogger.userinfo(mylog, 'info') --> print to screen (if quiet==False)
                                     and to logfile
"""
import logging
from socket import gethostname
import copy

def init_logger(logfilename, quiet=False, debug=False):
    logging.USERINFO = logging.INFO + 1
    logging.addLevelName(logging.USERINFO, 'USERINFO')
    logger = logging.getLogger('PyBDSF')
    logger.setLevel(logging.DEBUG)

    # First remove any existing handlers (in case PyBDSF has been run
    # before in this session but the quiet or debug options have changed
    while len(logger.handlers) > 0:
        logger.removeHandler(logger.handlers[0])

    # File handlers
    fh = ColorStripperHandler(logfilename)
    if debug:
        # For log file and debug on, print name and levelname
        fh.setLevel(logging.DEBUG)
        fmt1 = MultiLineFormatter('%(asctime)s %(name)-20s:: %(levelname)-8s: '\
                                 '%(message)s',
                                 datefmt='%a %d-%m-%Y %H:%M:%S')
    else:
        # For log file and debug off, don't print name and levelname as
        # they have no meaning to the user.
        fh.setLevel(logging.INFO)
        fmt1 = MultiLineFormatter('%(asctime)s:: %(levelname)-8s: %(message)s',
                                 datefmt='%a %d-%m-%Y %H:%M:%S')
    fh.setFormatter(fmt1)
    logger.addHandler(fh)

    # Console handler for warning, error, and critical: format includes levelname
    # ANSI colors are used
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    fmt2 = logging.Formatter('\033[31;1m%(levelname)s\033[0m: %(message)s')
    ch.setFormatter(fmt2)
    logger.addHandler(ch)

    # Console handler for USERINFO only: format does not include levelname
    # (the user does not need to see the levelname, as it has no meaning to them)
    # ANSI colors are allowed
    chi = logging.StreamHandler()
    chi.addFilter(InfoFilter())
    if quiet:
        # prints nothing, since filter lets only USERINFO through
        chi.setLevel(logging.WARNING)
    else:
        # prints only USERINFO
        chi.setLevel(logging.USERINFO)
    fmt3 = logging.Formatter('%(message)s')
    chi.setFormatter(fmt3)
    logger.addHandler(chi)

class InfoFilter(logging.Filter):
    # Lets only USERINFO through
    def filter(self, rec):
        return rec.levelno == logging.USERINFO

class MultiLineFormatter(logging.Formatter):
    def format(self, record):
        str = logging.Formatter.format(self, record)
        header, footer = str.split(record.message)
        nocolor_header = strip_color(header)
        str = str.replace('\n', '\n' + ' '*len(nocolor_header))
        return str

def userinfo(mylog, desc_str, val_str=''):
    """Writes a nicely formatted string to the log file and console

    mylog = logger
    desc_str = description string / message
    val_str = value string

    Message is constructed as:
      'desc_str .... : val_str'
    """
    bc = '\033[1;34m' # Blue
    nc = '\033[0m'    # Normal text color

    if val_str == '':
        sep = ''
        if desc_str[:1] == '\n':
            bc += '\n'
            desc_str = desc_str[1:]
        desc_str = bc + '--> ' + desc_str + nc
    else:
        sep = ' : '
        if len(desc_str) < 40:
            desc_str += ' '
        if len(desc_str) < 40:
            while len(desc_str) < 41:
                desc_str += '.'
        else:
            while len(desc_str) < 41:
                desc_str += ' '
    mylog.log(logging.USERINFO, desc_str+sep+val_str)


class ColorStripperHandler(logging.FileHandler):
    def emit(self, record):
        """Strips ANSI color codes from file stream"""
        myrecord = copy.copy(record)
        nocolor_msg = strip_color(myrecord.msg)
        myrecord.msg = nocolor_msg
        logging.FileHandler.emit(self, myrecord)

def strip_color(msg):
    """Strips specific ANSI color codes from an input string

    The color codes are hard-coded to those used above
    in userinfo() and in WARNING, ERROR, and CRITICAL.
    """
    nocolor_msg = ''
    a = msg.split('\033[1;34m')
    for b in a:
        c = b.split('\033[0m')
        for d in c:
            e = d.split('\033[31;1m')
            for f in e:
                nocolor_msg += f
    return nocolor_msg

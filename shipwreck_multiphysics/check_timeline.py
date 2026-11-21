"""Match high-precision ledger timestamps to OpenFOAM's printed time names."""
import math
import re


TIME_NAME = re.compile(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?')


def time_directory(parent, time):
    """Find the unique directory written for this printed OpenFOAM time."""
    if not parent.is_dir():
        return None
    matches = (path for path in parent.iterdir()
               if TIME_NAME.fullmatch(path.name) and float(path.name) == time and path.is_dir())
    match = next(matches, None)
    return match if next(matches, None) is None else None


def clock_settings(case):
    control = (case/'system/controlDict').read_text()
    precision = re.search(r'\btimePrecision\s+(\d+)\s*;',control)
    mode = re.search(r'\btimeFormat\s+(\w+)\s*;',control)
    return int(precision.group(1)) if precision else 6, mode.group(1) if mode else 'general'


def same_clock(logged, exact, precision, mode):
    """Allow half a printed decimal unit and variance from float parsing."""
    if not all(math.isfinite(t) for t in (logged,exact)):
        return False
    if logged==exact:
        return True
    if mode=='fixed':
        exponent = -precision
    elif mode in ('general','scientific') and logged:
        exponent = math.floor(math.log10(abs(logged))) - precision
        if mode=='general':
            exponent += 1
    else:
        return False
    return abs(logged-exact) <= .5*10.**exponent + 2*max(math.ulp(logged),math.ulp(exact))


def matches_timeline(logged, records, precision, mode):
    if not logged or len(logged)!=len(records):
        return False
    previous_logged = previous_actual = -math.inf
    for printed, row in zip(logged, records):
        actual = row['time_s']
        if (printed <= previous_logged or actual <= previous_actual
                or not same_clock(printed, actual, precision, mode)):
            return False
        previous_logged, previous_actual = printed, actual
    return True

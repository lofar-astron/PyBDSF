#!/usr/bin/env python
"""Print the parts of a cProfile file used to analyse the rmsimage bottleneck

Usage: python analyze_profile.py [profile file, default ../filter_skymodel.prof]
"""
import os
import pstats
import sys

prof = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', 'filter_skymodel.prof')
p = pstats.Stats(prof)
p.sort_stats('cumulative').print_stats(35)
p.sort_stats('tottime').print_stats(25)
p.sort_stats('cumulative')
p.print_callees('rms_mean_map')
p.print_callees('rmsimage.py:135')
p.print_callers('parallel_map')
p.print_callers('posix.fork')
p.print_callers('label')
p.print_callers('numpy.array')

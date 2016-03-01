import tracelabeler
import gaussmix
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from optparse import OptionParser

parser = OptionParser("tracelabeler-windowTraceview")
parser.add_option("--trace", type="string", dest="trace")
parser.add_option("--dme")
parser.add_option("--unbam")
parser.add_option("--zmw")
parser.add_option("--targetBase")
parser.add_option("--outputprefix")
(options, args) = parser.parse_args()
if not options.trace:
    parser.print_help()
    sys.exit(1)

# tl = tracelabeler.tracelabeler( trace="/pbi/collections/315/3150005/r54004_20151201_015856/1_A01/m54004_151201_015904.trc.h5", \
# dme="/home/UNIXHOME/dalexander/Projects/SequelFrenzy/All4Mers-3150005-0026-bakeoff/T2B-mainline-dmeDump/m54004_151201_015904.dme-dump.h5", \
# unbam="/home/UNIXHOME/dalexander/Projects/SequelFrenzy/All4Mers-3150005-0026-bakeoff/T2B-mainline-dmeDump/m54004_151201_015904.subreads.bam", \

tl = tracelabeler.tracelabeler( trace=options.trace, dme=options.dme, unbam=options.unbam, albam=None, ref=None)

# set the zmw and targetbase to get data from subread
tl.setzmw(int(options.zmw), int(options.targetBase))

mysubreadid = tl.unbam.query_name

# step through every base in the subread specified by the zmw and targetBase
# and compute the mean green and red signal (not baseline subtracted as stored in dme....)
ofp = open("%s.meangr.dat" % options.outputprefix, "w")
ofp.write("query_name\tindex\tbase\tsf\tpw\tmeangreen\tmeanred\n")
for bb in range(len(tl.rseq)):
    base = tl.rseq[bb]
    sf = tl.rsf[bb]
    pw = tl.rpw[bb]
    tr = tl.traceraw(sf,(sf+pw))
    dd = tr[0,:]
    meangreen=(sum(dd)/float(len(dd)))
    dd = tr[1,:]
    meanred=(sum(dd)/float(len(dd)))
    ofp.write( "%s\t%d\t%s\t%d\t%d\t%f\t%f\n" % (mysubreadid, bb, base, sf, pw, meangreen, meanred))
ofp.close()

############################################################################
# Now do the same but with baseline subtracted.

# this is somewhat cumbersome as the dme windows aren't always nice
# consistent lengths. This is taken care of behind the scences

ofp = open("%s.meangr.blsub.dat" % options.outputprefix, "w")
ofp.write("query_name\tindex\tbase\tsf\tpw\tmeangreen\tmeanred\n")
for bb in range(len(tl.rseq)):
    base = tl.rseq[bb]
    sf = tl.rsf[bb]
    pw = tl.rpw[bb]

    tr = tl.tracebasesub(sf,(sf+pw))

    dd = tr[0,:]
    meangreen=(sum(dd)/float(len(dd)))
    dd = tr[1,:]
    meanred=(sum(dd)/float(len(dd)))
    ofp.write( "%s\t%d\t%s\t%d\t%d\t%f\t%f\n" % (mysubreadid, bb, base, sf, pw, meangreen, meanred))
ofp.close()



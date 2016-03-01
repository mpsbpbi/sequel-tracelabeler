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
parser.add_option("--albam")
parser.add_option("--ref")
parser.add_option("--zmw")
parser.add_option("--targetBase")
parser.add_option("--framewindow")
parser.add_option("--outputprefix")
(options, args) = parser.parse_args()
if not options.trace:
    parser.print_help()
    sys.exit(1)

framewindow = int(options.framewindow)
outprefix = options.outputprefix

# tl = tracelabeler.tracelabeler( trace="/pbi/collections/315/3150005/r54004_20151201_015856/1_A01/m54004_151201_015904.trc.h5", \
# dme="/home/UNIXHOME/dalexander/Projects/SequelFrenzy/All4Mers-3150005-0026-bakeoff/T2B-mainline-dmeDump/m54004_151201_015904.dme-dump.h5", \
# unbam="/home/UNIXHOME/dalexander/Projects/SequelFrenzy/All4Mers-3150005-0026-bakeoff/T2B-mainline-dmeDump/m54004_151201_015904.subreads.bam", \
# albam="/home/UNIXHOME/mbrown/mbrown/workspace2016Q1/sequel-squash/rt.subreads.all4mer44.bam", \
# ref="/mnt/secondary/iSmrtanalysis/install/smrtanalysis_2.4.0.140820/common/references/All4mer_V2_44_circular_72x_l50256/sequence/All4mer_V2_44_circular_72x_l50256.fasta")

tl = tracelabeler.tracelabeler( trace=options.trace, dme=options.dme, unbam=options.unbam, albam=options.albam, ref=options.ref)

# set the zmw and targetbase to get data
tl.setzmw(int(options.zmw), int(options.targetBase))
# compute alignment correspondances as we have them
tl.alignCorresp()

# what is min/max aligned read base?
print "For selected subread min/max aligned read base indices", tl.alignReadStart, tl.alignReadEnd

# what DME windows contain the aligned bases?
windowMin = tl.readBaseToWindowStartFrame(tl.alignReadStart)
windowMax = tl.readBaseToWindowStartFrame(tl.alignReadEnd-1)
print "windowMin/max", windowMin, windowMax

# check to make sure the given framewindow makes sense
if framewindow>(windowMax-windowMin):
    print "framewindow too big:", framewindow, ">" , (windowMax-windowMin)
    sys.exit(1)

thiswindow = windowMin + framewindow
(startframe, endframe) = tl.dmeWindowToStartEndFrame(thiswindow)

# For the frames in this window what are the start/end reads?
(readStart, readEnd) = tl.frameWindowToReadStartEnd(thiswindow)
# make sure we are looking at the aligned bases
readStart = max(readStart, tl.alignReadStart)
readEnd = min(readEnd, tl.alignReadEnd-1)
print "read start end", readStart, readEnd, "framewindow=", startframe, endframe

# now for read bases in consideration collect startframe pulsewidth base and alignment
windowdat = tl.subreadDat(readStart,readEnd)
print(windowdat["doc"])
for dd in windowdat["dat"]:
    print " ".join([str(xx) for xx in dd])

# look at the dme data for this window
dmedat = tl.windowToDMEDat(thiswindow)
print "dmedat", dmedat

gm = gaussmix.gaussmix(dmedat["MixtureFraction"], dmedat["Mean"], dmedat["Covariance"], dmedat["BaselineMean"])

# compute the qv of the dme base across all frames. Note we chose
# everything so it is in a single dme frame
qv = []
for ff in range(startframe,endframe):
    qv.append(gm.qv(gm.logcompprob(tl.trace[0,ff], tl.trace[1,ff])))

#### Now traceview it to make sure everything lines up.

# plot the trace data
plt.figure(figsize=(16,2))
for dd in windowdat["dat"]:
    plt.plot([dd[1],dd[1]],[0,450],'c-', linewidth=0.2) # pulse start
    plt.plot([dd[1]+dd[2],dd[1]+dd[2]],[0,450],'y-', linewidth=0.2) # pulse end

    plt.text(dd[1],400,dd[6],fontsize=5) # ref
    plt.text(dd[1],350,dd[5],fontsize=5) # read

# green and red trace data
plt.plot(range(startframe,endframe), tl.trace[0,startframe:endframe], 'g-', linewidth=0.2)
plt.plot(range(startframe,endframe), tl.trace[1,startframe:endframe], 'r-', linewidth=0.2)

# qvs
plt.plot(range(startframe,endframe), qv, 'k-', linewidth=0.2)

plt.savefig('%s-squashtracelabler' % outprefix,dpi=500)
plt.close()


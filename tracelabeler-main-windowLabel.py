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

################################
# given reference compute the hmm
# hmm model parameters
pwbasemean = {'A': 16.0, 'C': 16.0, 'T': 16.0, 'G': 16.0}
ipdmean = 32.0
hmmst = {} # to,from
hmmst["-A"] = 1.0/(pwbasemean["A"]+1)
hmmst["AA"] = 1.0-hmmst["-A"]
hmmst["-C"] = 1.0/(pwbasemean["C"]+1)
hmmst["CC"] = 1.0-hmmst["-C"]
hmmst["-G"] = 1.0/(pwbasemean["G"]+1)
hmmst["GG"] = 1.0-hmmst["-G"]
hmmst["-T"] = 1.0/(pwbasemean["T"]+1)
hmmst["TT"] = 1.0-hmmst["-T"]
hmmst["--"] = 1.0-1.0/(ipdmean+1)
hmmst["A-"] = (1.0-hmmst["--"])/4
hmmst["C-"] = (1.0-hmmst["--"])/4
hmmst["G-"] = (1.0-hmmst["--"])/4
hmmst["T-"] = (1.0-hmmst["--"])/4
print "hmmst", hmmst

baseIdx = {"-":0,"A":1,"C":2,"G":3,"T":4}

# Get the reference in the subread window and construct hmmstates from it. Pad then use local alignment.

#old:
#ref = reduce( lambda x,y: x+y, map( lambda ii: dat[ii][6], range(len(dat))))
#ref = ref.replace("-","")

dat = windowdat["dat"]
pad = 5
if tl.alignStrand == 1:
    reffirst = dat[0][3][0]
    reflast = dat[-1][3][0]
    refst = max(0, reffirst-pad)
    refend = min(len(tl.ref), reflast+pad)
    ref = tl.ref[refst:refend]
else:
    reffirst = dat[-1][3][0]
    reflast = dat[0][3][0]
    refst = max(0, reffirst-pad)
    refend = min(len(tl.ref), reflast+pad)
    ref = tl.revcomp(tl.ref[refst:refend])

print "ref", refst, refend, ref

hmmstates = ["-"]
for rr in ref:
    hmmstates.append(rr)
    hmmstates.append("-")

# the dme cluter probs for the frame pulse data
gmp = tl.dmeprobs(startframe, endframe)
(traceback, tracebackscore, newbasecalls) = tl.computehmm( hmmstates, hmmst, gmp, False, startframe) # local alignment

print "new basecalls"
for ii in range(len(newbasecalls)):
    print ii, newbasecalls[ii]

#### plot the new basecalls
# plot the trace data
plt.figure(figsize=(16,2))
for dd in newbasecalls:
    plt.plot([dd[1],dd[1]],[0,450],'c-', linewidth=0.2) # pulse start
    plt.plot([dd[1]+dd[2],dd[1]+dd[2]],[0,450],'y-', linewidth=0.2) # pulse end
    plt.text(dd[1],400,dd[0],fontsize=5) # ref
# green and red trace data
plt.plot(range(startframe,endframe), tl.trace[0,startframe:endframe], 'g-', linewidth=0.2)
plt.plot(range(startframe,endframe), tl.trace[1,startframe:endframe], 'r-', linewidth=0.2)

plt.savefig('%s-squashnewtracelabler'% outprefix,dpi=500)
plt.close()

import tracelabeler
import gaussmix
import matplotlib.pyplot as plt
import numpy as np
import math
import sys
from optparse import OptionParser

def sortindex(seq):
    return [x for x,y in sorted(enumerate(seq), key = lambda x: -x[1])]

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
found = tl.setzmw(int(options.zmw), int(options.targetBase))
if not found: sys.exit(1)

myqueryid = tl.unbam.query_name

############################################################################

# compute QVs for each of the called bases

# step through each new called base and look at the trace data between
# the start of the base to the start of the next base.

# compute 
# - mismatch: base=T alt=A,C,G
# - insert:   base=T alt=TA,TC,TG,TT
# - delete:   base=T alt=-

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

allqvsHeader = ["A","C","G","T", "IA", "IC","IG","IT","D"]
allqvsHeaderToIdx = dict([(kk,vv) for (vv,kk) in enumerate(allqvsHeader)])
allqvs = []

for bb in range(len(tl.rseq)-1):
    # TODO: doesn't score last base (-1) as I don't know where the "end" is
    thisbase = tl.rseq[bb]
    trstart = tl.rsf[bb]
    trend = tl.rsf[bb+1]

    # the dme cluter probs for the frame pulse data
    gmp = tl.dmeprobs(trstart, trend)

    newqvs = []
    #### mm
    (tr, trsc, bc) = tl.computehmm( ["A","-"], hmmst, gmp )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( ["C","-"], hmmst, gmp )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( ["G","-"], hmmst, gmp )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( ["T","-"], hmmst, gmp )
    newqvs.append(trsc[-1])

    #### ins
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "A","-"], hmmst, gmp )
    newqvs.append(trsc[-1]-math.log(8))
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "C","-"], hmmst, gmp )
    newqvs.append(trsc[-1]-math.log(8))
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "G","-"], hmmst, gmp )
    newqvs.append(trsc[-1]-math.log(8))
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "T","-"], hmmst, gmp )
    newqvs.append(trsc[-1]-math.log(8))

    #### del
    (tr, trsc, bc) = tl.computehmm( ["-"], hmmst, gmp )
    newqvs.append(trsc[-1]-math.log(0.25))

    allqvs.append(newqvs)

ofp = open("%s.labelqv.dat" % options.outputprefix, "w")
print >>ofp, "myqueryid index base sf alt qv"
for bb in range(len(tl.rseq)-1):
    thisbase = tl.rseq[bb]
    trstart = tl.rsf[bb]
    trend = tl.rsf[bb+1]

    print >>ofp, myqueryid, bb, thisbase, trstart,
#    for vv in allqvs[bb]:
#        print >>ofp, vv,
    thislp = allqvs[bb][allqvsHeaderToIdx[thisbase]]
    si = sortindex(allqvs[bb])
    if allqvsHeader[si[0]]==thisbase:
        nextlp = allqvs[bb][si[1]]
        nexttype = allqvsHeader[si[1]]
    else:
        # basecaller disagrees with this best
        nextlp = allqvs[bb][si[0]]
        nexttype = allqvsHeader[si[0]]
    qv = 10*(thislp-nextlp)/math.log(10.0)
    print >>ofp, nexttype, qv
ofp.close()


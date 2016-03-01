import tracelabeler
import gaussmix
import matplotlib.pyplot as plt
import numpy as np
import math

tl = tracelabeler.tracelabeler( trace="/pbi/dept/primary/traceSim/Sequel/RTvsOfflinePipeline/LVP1+releaseChemistryTestSet/TestSet_004_rs200.trc.h5", \
dme="/home/UNIXHOME/mbrown/mbrown/workspace2016Q1/sequel-tracelabeler/TestSet_004_rs200.dme-dump.h5", \
unbam="/home/UNIXHOME/mbrown/mbrown/workspace2016Q1/sequel-tracelabeler/TestSet_004_rs200.subreads.bam", \
albam="/home/UNIXHOME/mbrown/mbrown/workspace2016Q1/sequel-tracelabeler/current-TestSet_004_rs200_32x.pbalign.bam", \
ref="/home/UNIXHOME/mbrown/mbrown/workspace2016Q1/dmesim/TestSet_004_rs200_32x.fasta")

# the first zmw is 2097184 in this sim data
tl.setzmw(2097184)
tl.alignCorresp()

# what is min/max aligned read base?
print "min/maxAlignedReadBase", tl.alignReadStart, tl.alignReadEnd
# min/maxAlignedReadBase 122 3450

# what 1k window contains the first/last base?
windowMin = tl.readBaseToWindowStartFrame(tl.alignReadStart)
windowMax = tl.readBaseToWindowStartFrame(tl.alignReadEnd)
print "windowMin/max", windowMin, windowMax
# windowMin/max 8 98

####
thiswindow = windowMin + 89 # +89 
startframe = 1024*thiswindow
endframe = 1024*(thiswindow+1)

# For the frames in this window what are the start/end reads?
(readStart, readEnd) = tl.frameWindowToReadStartEnd(thiswindow)
print "read start end", readStart, readEnd, "framewindow=", 1024*thiswindow, 1024*(thiswindow+1)
# read start end 3404 3440 framewindow= 99328 100352

# now for read bases collect startframe pulsewidth base
windowdat = tl.subreadDat(readStart,readEnd)
print(windowdat["doc"])
for dd in windowdat["dat"]:
    print " ".join([str(xx) for xx in dd])

"""readbaseindex startframe pulsewidth refindex alignindex readbase alignrefbases
3404 99362 22 [4436] [3327] C C
3405 99390 14 [4435] [3328] C C
3406 99419 6 [4434] [3329] A A
3407 99434 12 [4433] [3330] A A
3408 99448 7 [4432] [3331] A A
3409 99491 16 [4431] [3332] T T
3410 99521 10 [4430] [3333] G G
3411 99545 2 [4429] [3334] A A
3412 99566 8 [4428] [3335] C C
3413 99586 8 [4427] [3336] G G
3414 99605 26 [4426] [3337] A A
3415 99640 9 [4425] [3338] A C
3416 99658 34 [4424] [3339] T T
3417 99718 4 [4423] [3340] T T
3418 99731 47 [4422] [3341] C C
3419 99784 14 [4421] [3342] T T
3420 99799 9 [4420] [3343] A A
3421 99811 18 [4419] [3344] C C
3422 99837 11 [4418] [3345] C C
3423 99855 14 [4417] [3346] A A
3424 99874 7 [4416] [3347] C C
3425 99888 5 [4415] [3348] A A
3426 99910 10 [4414] [3349] T T
3427 99966 8 [4413] [3350] C C
3428 100002 4 [4413] [3351] T -
3429 100016 7 [4412] [3352] T T
3430 100058 17 [4411] [3353] A A
3431 100117 8 [4410] [3354] T T
3432 100194 1 [4409] [3355] G T
3433 100198 5 [4408] [3356] G G
3434 100234 13 [4407] [3357] A A
3435 100263 12 [4406] [3358] C C
3436 100304 15 [4405] [3359] A A
3437 100326 3 [4404] [3360] G T
3438 100333 9 [4403] [3361] T T
3439 100350 25 [4402] [3362] A A
"""

dmedat = tl.windowToDMEDat(thiswindow)
"""{'Covariance':
array([[  355.7074585 ,   341.05419922,    67.3758316 ],
       [  406.21313477,   593.21142578,   110.36556244],
       [  483.96380615,  1116.84570312,   239.33477783],
       [  647.30444336,   441.50177002,   146.74969482],
       [  955.24627686,   531.82269287,   261.15966797]], dtype=float32),
'MixtureFraction': array([ 0.59466136,  0.06441299,  0.11632776,  0.13682321,  0.08777492], dtype=float32),
'Mean': array([[  -0.94825631,    2.05985498],
       [  36.88322449,  116.41873169],
       [  73.76644897,  232.83746338],
       [ 127.84963989,   61.96887589],
       [ 199.76504517,   96.82637024]], dtype=float32)}
"""
gm = gaussmix.gaussmix(dmedat["MixtureFraction"], dmedat["Mean"], dmedat["Covariance"], dmedat["BaselineMean"])

#### Now traceview it to make sure everything lines up.

# compute the qv of the dme base across all frames
qv = []
for ff in range(startframe,endframe):
    qv.append(gm.qv(gm.logcompprob(tl.trace[0,ff], tl.trace[1,ff])))

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

plt.savefig('tracelabler',dpi=500)
plt.close()

################################################################
################################################################

#### Now all I have to do is an super simple hmm... Because there is a
#### 1-1 between frame and hmm state with one output per transistion
#### everything is EASY!

#### ingredients:
# trace data: tl.trace[0,startframe:endframe]
# prob(base|traceFrame) from dme: gm.logcompprob(tl.trace[0,ff], tl.trace[1,ff])
# the true base sequence in alignment with basecalled read: windowdat[ii][6] # alignrefbases. can be "-"
# the basecalled read: windowdat[ii][5] #readbase along with pulse position and width

# the simplest hmm: star topology with null in center with simple
# geometric distributions: only 5 geometric parameters assuming all
# bases equally likely. Each state outputs frame with probabilities
# defined by DME

# hmmst["tostate,fromstate"] = P
# hmmout["state", x0, x1] = P

# The DP table has size (#frames, 2*#alignrefbases+1). Assumes the
# model retruns to baseline after each base. Filling out this small
# table and computing the Viterbi path should be fast.

# The only estimation are the geometrics in the simple hmm. We could
# use baum-welsch. Just use methods of moments.

# https://en.wikipedia.org/wiki/Geometric_distribution: P(k in
# 0,1,2,3... | p) = (1-p)^k*p. mean=(1-p)/p. So if an A-pulse had a
# mean length of 10 then: M=(1-p)/p. p=1/(M+1) = 1/11

# let's overfit to this window

dat=windowdat["dat"]
pwbasemean = {}
for mybase in ["A","C","G","T"]:
    correct = filter( lambda ii: ((dat[ii][6]==dat[ii][5]) and (dat[ii][5]==mybase)), range(len(dat)) )
    dd = [ dat[cc][2] for cc in correct ]
    pwbasemean[mybase]=sum(dd)/float(len(dd))
# pwbasemean
# {'A': 12.583333333333334, 'C': 16.333333333333332, 'T': 12.75, 'G': 7.666666666666667}

# now compute ipd given both bases correct (ii and ii+1)
# ipd = sf[ii+1]-(sf[ii]+pw[ii])
correct = filter( lambda ii: ((dat[ii][6]==dat[ii][5]) and (dat[ii+1][6]==dat[ii+1][5])), range(len(dat)-1) )
ipds = [ dat[ii+1][1]-(dat[ii][1]+dat[ii][2]) for ii in correct]
ipdmean = sum(ipds)/float(len(ipds))
# 16.074074074074073

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

# Get the reference in the subread window and construct hmmstates from it
ref = reduce( lambda x,y: x+y, map( lambda ii: dat[ii][6], range(len(dat))))
ref = ref.replace("-","")
hmmstates = ["-"]
for rr in ref:
    hmmstates.append(rr)
    hmmstates.append("-")

################################
# compute the hmm
(traceback, tracebackscore, newbasecalls) = tl.computehmm( hmmstates, hmmst, gm, startframe, endframe )

# plot dme to make sure it make sense
dme = []
for ff in range(startframe,endframe):
    dme.append(gm.logcompprob(tl.trace[0,ff], tl.trace[1,ff]))
plt.figure()
plt.plot(range(0,0+endframe-startframe), [(dme[ii][0]) for ii in range(endframe-startframe)], 'k-', linewidth=0.2)
plt.plot(range(0,0+endframe-startframe), [(dme[ii][1]) for ii in range(endframe-startframe)], 'r-', linewidth=0.2)
plt.plot(range(0,0+endframe-startframe), [(dme[ii][2]) for ii in range(endframe-startframe)], 'g-', linewidth=0.2)
plt.plot(range(0,0+endframe-startframe), [(dme[ii][3]) for ii in range(endframe-startframe)], 'b-', linewidth=0.2)
plt.plot(range(0,0+endframe-startframe), [(dme[ii][4]) for ii in range(endframe-startframe)], 'c-', linewidth=0.2)
plt.savefig('dmeprobs')

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
# qvs
plt.plot(range(startframe,endframe), qv, 'k-', linewidth=0.2)

plt.savefig('newtracelabler',dpi=500)
plt.close()

################################

# compute QVs for each of the called bases

# step through each new called base and look at the trace data between
# the start of the base to the start of the next base.

# compute 
# - mismatch: base=T alt=A,C,G
# - insert:   base=T alt=TA,TC,TG,TT
# - delete:   base=T alt=-

allqvsHeader = ["MMA","MMC","MMG","MMT", "INSA", "INSC","INSG","INST","DEL"]
allqvsHeaderToIdx = dict([(kk,vv) for (vv,kk) in enumerate(allqvsHeader)])
allqvs = []
for ii in range(len(newbasecalls)-1):
    newqvs = []
    thisbase = newbasecalls[ii][0]
    trstart = newbasecalls[ii][1]
    trend = newbasecalls[ii+1][1]-1

    #### mm
    (tr, trsc, bc) = tl.computehmm( ["A","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( ["C","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( ["G","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( ["T","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])

    #### ins
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "A","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "C","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "G","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])
    (tr, trsc, bc) = tl.computehmm( [thisbase, "-", "T","-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])

    #### del
    (tr, trsc, bc) = tl.computehmm( ["-"], hmmst, gm, trstart, trend )
    newqvs.append(trsc[-1])

    allqvs.append(newqvs)

for ii in range(len(newbasecalls)-1):
    thisbase = newbasecalls[ii][0]
    trstart = newbasecalls[ii][1]
    trend = newbasecalls[ii+1][1]-1
    print "calling tr", thisbase, trstart, trend,
    for vv in allqvs[ii]:
        print vv,
    myargmax = filter( lambda jj: allqvs[ii][jj]==max(allqvs[ii]), range(len(allqvs[ii])))
    print "best", allqvsHeader[myargmax[0]],
    print allqvsHeader[myargmax[0]]=="MM%s" % thisbase,
    ss = sorted(allqvs[ii], key= lambda x: -x)
    qv = 10*(ss[0]-ss[1])/math.log(10.0)
    print "qv", qv

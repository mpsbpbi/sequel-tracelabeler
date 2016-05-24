import h5py
import numpy as np
import bisect
import sys
import pysam
import pbcore.io
import math
import gaussmix
import os.path

################################
def AContainedInB(asta, aend, bsta, bend):
    return( (asta >= bsta) and (asta<=bend) and (aend >= bsta) and (aend<=bend) )

################################
class tracelabeler:
    """goal: label every frame of a trace correctly

- basecaller take trace, model and output basecalls, start, width,
confidence

- tracelabeler takes trace, model, reference, approxAlign and output
refBases, start, width, confidence.

the tracelabeler is the basecaller HMM "unrolled" using the reference
with the approxAlign used to limit the computation in the dynamic
program.

Note you can do baum-welsch like parameter estimates by running em. Or
simply tabulate statistics from which the model parameters are set viz
(P(base=A, pw=2).
"""

    rMap = dict(zip("ACGTacgt-","TGCAtgca-"))

    ################################
    def revcomp(self, instr):
        return("".join([tracelabeler.rMap[c] for c in instr[::-1]]))

    ################################
    def __init__( self, trace="tracefile", dme=None, unbam=None, albam=None, ref=None):
        self.tracef=trace
        self.dmef=dme
        self.unbamf=unbam
        self.albamf=albam
        self.reff=ref

        # do timestamp to make sure trace < dme < unbam < albam
        timesOK=True
        t0 = os.path.getmtime( self.tracef)
        if self.dmef is not None:
            t1 = os.path.getmtime( self.dmef)
            if not t1>t0:
                timesOK=False
                print >>sys.stderr, "ERROR: dme file is before trace!", t1, t0, self.dmef, self.tracef
        if self.unbamf is not None:
            t2 = os.path.getmtime( self.unbamf)
            if not t2>t0:
                timesOK=False
                print >>sys.stderr, "ERROR: unbamf file is before trace!", t2, t0, self.unbamf, self.tracef
            if not t2>t1:
                timesOK=False
                print >>sys.stderr, "ERROR: unbamf file is before dme!", t2, t1, self.unbamf, self.dmef
        if self.albamf is not None:
            t3 = os.path.getmtime( self.albamf)
            if not t3>t0:
                timesOK=False
                print >>sys.stderr, "ERROR: albamf file is before trace!", t2, t0, self.albamf, self.tracef
            if not t3>t1:
                timesOK=False
                print >>sys.stderr, "ERROR: albamf file is before dme!", t2, t1, self.albamf, self.dmef
            if not t3>t2:
                timesOK=False
                print >>sys.stderr, "ERROR: albamf file is before unbam!", t2, t1, self.albamf, self.unbamf
        if timesOK:
            print >>sys.stderr, "tracelabeler. file timestamps in correct order"

    ################################
    def setzmw( self, zmw, targetbase=0 ):
        "Set the zmw subread that contains targetbase (bam subreads) and pull all the data"

        self.zmw = zmw
        self.targetbase = targetbase

        #### trace data
        self.traceff=h5py.File(self.tracef,"r")
        self.tridx = bisect.bisect(self.traceff['TraceData']['HoleNumber'],self.zmw)-1
        if (self.zmw!=self.traceff['TraceData']['HoleNumber'][self.tridx]):
            print "ERROR: zmw not found. assert tracedata zmw=", self.zmw, "at index", self.tridx, "not equal", self.traceff['TraceData']['HoleNumber'][self.tridx]
            return(False)
        self.trace = self.traceff['TraceData']['Traces'][self.tridx]
        print "trace got data", zmw, targetbase, self.traceff['TraceData']['HoleNumber'][self.tridx]

        #### dme
        self.dmeStartFrame=None
        self.dmeEndFrame=None
        self.dmeBlockSize=None
        self.dmeIndex=None
        if self.dmef is not None:
            self.dmeff=h5py.File(self.dmef,"r")
            self.dmeidx = bisect.bisect(self.dmeff['HoleNumber'],self.zmw)-1
            assert(self.zmw==self.dmeff['HoleNumber'][self.dmeidx])
            # some dme blocks have BlockSize==0 and aren't used. (TODO: ???)
            self.dmeStartFrame = []
            self.dmeEndFrame = []
            self.dmeBlockSize = []
            self.dmeIndex = []
            for ii in range(len(self.dmeff['BlockSize'])):
                if self.dmeff['BlockSize'][ii]>0:
                    self.dmeBlockSize.append(self.dmeff['BlockSize'][ii])
                    self.dmeStartFrame.append(self.dmeff['StartFrame'][ii])
                    self.dmeEndFrame.append(self.dmeff['EndFrame'][ii])
                    # make sure dme blocks are ordered increasing
                    if len(self.dmeStartFrame)>1:
                        assert( self.dmeStartFrame[-1] == (self.dmeEndFrame[-2]))
                    self.dmeIndex.append(ii)
            # dmeWindow is contiguous 0:len and Index is into the DME data that has holes
            self.dmeWindowToIndex = dict(enumerate(self.dmeIndex))
            print "assert dme lengths:", self.dmeEndFrame[-1], len(self.trace[0])
            assert(self.dmeEndFrame[-1] == len(self.trace[0]))
            print "dme got data"

        #### unbam
        self.unbam = None
        self.rsf = None
        self.rpw = None
        self.rseq = None
        if self.unbamf is not None:
            unbamff = pysam.AlignmentFile(self.unbamf,"rb", check_sq=False)
            for dd in unbamff:
                # with subread bam there might be multiple subreads with zmw take one containing targetbase
                if ("/%d/" % self.zmw) in dd.query_name:
                    # TODO: do I really have to parse the subread location from the query_name???
                    (substart, subend) = [int(xx) for xx in dd.query_name.split("/")[-1].split("_")]
                    print "unbam found zmw:", dd.query_name,dd.qstart, dd.qend, dd.query_alignment_start, dd.query_alignment_end, substart,subend
                    if AContainedInB(targetbase,targetbase, substart,subend):
                        self.unbam = dd
                        print "unbam got %s" % dd.query_name
                        self.rsf = self.unbam.get_tag("sf") # start frame
                        self.rpw = self.unbam.get_tag("pw") # pulse width
                        self.rseq = self.unbam.seq
                        self.rsubstart = substart
                        self.rsubend = subend
                        break
            if not self.unbam:
                print "ERROR: unbam not found for %s!!!" % self.zmw
                return(False)

        #### albam
        self.albam = []
        if self.albamf is not None:
            # bam doesn't have the reference
            reader = pbcore.io.openIndexedAlignmentFile( self.albamf, self.reff )
            for h in reader:
                zmw=h.HoleNumber
                if zmw == self.zmw and AContainedInB(h.queryStart, h.queryEnd, self.rsubstart, self.rsubend):
                    self.albam.append(h)
                    numerr=(h.nMM + h.nIns + h.nDel)
                    rlen = (h.readEnd-h.readStart)
                    if numerr>0:
                      err = float(numerr)/float(rlen)
                    else:
                      err = 1.0/float(rlen+1)
                    if h.isReverseStrand:
                      tid = h.referenceName+"rc"
                    else:
                      tid = h.referenceName
                    moviename= h.movieName
                    rstart = h.readStart
                    rend = h.readEnd
                    qid = "%s/%s/%d_%d" % (moviename,zmw,rstart,rend)
                    read = h.read()
                    ref = h.reference()
                    print "albam got hit rlen=%d err=%f tid=%s qid=%s" % (rlen,err,tid,qid)
            if len(self.albam)==0:
                print "ERROR: albam not found for %s!!!" % self.unbam.query_name
                return(False)

        #### ref
        if self.reff is not None:
            myref = []
            dat = open(self.reff).read().splitlines()
            for ii in range(1,len(dat)):
                myref.append(dat[ii])
            self.ref = "".join(myref)
        
        return(True)

    ################################
    def alignCorresp( self ):
        "compute aligned correspondance between read and ref bases"

        # because we specify zmw and targetbase, this should be a single alignment even with multiple subreads in the subreads bam
        h = self.albam[0]

        strand = +1 if h.isForwardStrand else -1
        ts = h.tStart
        te = h.tEnd
        qs = h.aStart
        qe = h.aEnd

        # print for debug
        #tSeq = [tracelabeler.rMap[c] for c in h.reference()[::-1]]
        print "alignCorresp:", qs, qe, ts, te, strand
        #print h.reference()
        #print h.read()

        # swap strands if necessary to get ref coord
        tM =  [ 1 if x != '-' else 0 for x in h.reference()] 
        qM =  [ 1 if x != '-' else 0 for x in h.read()] 
        qOffset = qs-1 + np.cumsum(qM)
        if strand == +1:
            tOffset = ts-1 + np.cumsum(tM)
        else:
            tOffset = te - np.cumsum(tM)
        index = np.array(range(len(qOffset)))

        # map-ping: todo probably not too efficient
        self.alignReadStart = qs - self.rsubstart # TODO h.queryStart # subtract off start of subread
        self.alignReadEnd = qe - self.rsubstart # TODO h.queryStart # subtract off start of subread
        self.alignRefStart = ts
        self.alignRefEnd = te
        self.alignStrand = strand
        self.alignReadToRef = map(lambda ii: tOffset[qOffset==ii], range(qs,qe) )
        self.alignRefToRead = map(lambda ii: qOffset[tOffset==ii], range(ts,te) )
        self.alignReadToAlign = map(lambda ii: index[qOffset==ii], range(qs,qe) )
        self.alignRead = h.read()
        self.alignRef = h.reference()

    def readToRef( self, rr ):
        "correspondance between read base and ref base on already computed mapping"
        if rr<self.alignReadStart:
            return([])
        if rr>=self.alignReadEnd:
            return([])
            
        return(self.alignReadToRef[rr-self.alignReadStart])

    def readToAlign( self, rr ):
        "correspondance between read base and alignment position on already computed mapping"
        if rr<self.alignReadStart:
            return([0])
        if rr>=self.alignReadEnd:
            return([self.alignReadEnd])

        return(self.alignReadToAlign[rr-self.alignReadStart])
            
    ################################
    def frameWindowToReadStartEnd(self, window):
        "Given a DME window number return the readStart and readEnd of read bases in that start in that window of frames"
        readStart = bisect.bisect( self.rsf, self.dmeStartFrame[window])
        readEnd =   bisect.bisect( self.rsf, self.dmeEndFrame[window])
        return(readStart,readEnd)

    ################################
    def readBaseToWindowStartFrame(self, basenum):
        "Given read base, return window containing start frame"
        return( self.frameToDMEWindow(self.rsf[basenum]) )

    ################################
    def rawAlignRefReadFromReadSE(self, start=None, end=None):
        "raw alignment containing read start to end"
        if start is not None:
            als = self.readToAlign(start)[0]
            ale = self.readToAlign(end)[0]
        else:
            als = 0
            ale = len( self.albam[0].reference() )
        return( (self.albam[0].reference()[als:ale ], self.albam[0].read()[als:ale ]))

    ################################
    def subreadDat(self, start, end):
        "for subread start to end return the relevant data"
        result = {}
        result["dat"] = []
        result["doc"] = "readbaseindex startframe pulsewidth refindex alignindex readbase alignrefbases"
        for readindex in range(start,end):
            # TODO!!!!!!!!!: I can't get the mapping right to raw read!!!!!! 
            #resrseq = self.rseq[readindex] # +self.albam[0].aStart]
            resrseq = "".join( [self.alignRead[xx] for xx in self.readToAlign(readindex)])
            ##resreadToRef = "".join( [self.ref[xx] if self.alignStrand==1 else self.rMap[self.ref[xx]] for xx in self.readToRef(readindex)])
            resreadToAlign_ref = "".join( [self.alignRef[xx] for xx in self.readToAlign(readindex)])
            result["dat"].append([ readindex, self.rsf[readindex], self.rpw[readindex], self.readToRef(readindex), self.readToAlign(readindex), resrseq, resreadToAlign_ref ])
        return(result)
        
    ################################
    def frameToDMEWindow( self, frame):
        "Return the correct dme window containing the given frame"
        # ideally this should be floor(frame/BlockSize) but I guess anything can happen
        # the frame needs to occur between start and end
        right = bisect.bisect_right(self.dmeStartFrame,frame)-1
        #left  = bisect.bisect_left(self.dmeEndFrame,frame)
        #assert(left==right)
        return(right)

    ################################
    def dmeWindowToStartEndFrame( self, window):
        return( (self.dmeStartFrame[window], self.dmeEndFrame[window]) )

    ################################
    def windowToDMEDat( self, inwindow):
        # We_order_them_as:_[-,_A,_C,_G,_T]

        # map from contiguous window to DME index with holes
        window = self.dmeWindowToIndex[inwindow]

        MixtureFraction = self.dmeff["/SmoothedEstimates/MixtureFraction"][window][self.dmeidx]
        Mean = self.dmeff["/SmoothedEstimates/Mean"][window][self.dmeidx]
        Covariance = self.dmeff["/SmoothedEstimates/Covariance"][window][self.dmeidx]
        BaselineMean = self.dmeff["/BaselineMean"][window][self.dmeidx]
        return({"MixtureFraction": MixtureFraction, "Mean": Mean, "Covariance": Covariance, "BaselineMean": BaselineMean})

    ################################
    def traceraw( self, start, end):
        "raw trace data start to end"
        return(self.trace[:,start:end])

    ################################
    def tracebasesub( self, start, end):
        "trace data with baseline subtracted from dme. start to end might have multiple DME windows, deals with it"
        if ((end-start)<1): return([])
        firstdmewin = self.frameToDMEWindow(start)
        lastdmewin = self.frameToDMEWindow(end-1)
        currentdmewin = firstdmewin
        currentdmedat = self.windowToDMEDat(currentdmewin)
        if (firstdmewin==lastdmewin):
            # all data within one dme
            return((self.trace[:,start:end].transpose() - currentdmedat["BaselineMean"]).transpose())
        else:
            # data across multiple dme frames
            newstart = start
            newend = min(self.dmeEndFrame[firstdmewin], end)
            print >>sys.stderr, "tracesub", start, end, "not in one dme so computing", newstart, newend, "and then", newend, end
            dat = (self.trace[:,newstart:newend].transpose() - currentdmedat["BaselineMean"]).transpose()
            # recurse TODO remove recursion
            return(np.hstack((dat, self.tracebasesub(newend, end))))

    ################################
    def dmeprobs( self, start, end):
        "dme probabilities. start to end might have multiple DME windows, deals with it"
        if ((end-start)<1): return([])
        firstdmewin = self.frameToDMEWindow(start)
        lastdmewin = self.frameToDMEWindow(end-1)
        currentdmewin = firstdmewin
        currentdmedat = self.windowToDMEDat(currentdmewin)
        currentgm = gaussmix.gaussmix(currentdmedat["MixtureFraction"], currentdmedat["Mean"], currentdmedat["Covariance"], currentdmedat["BaselineMean"])
        if (firstdmewin==lastdmewin):
            # all data within one dme
            res = []
            for ii in range(start,end):
                res.append( currentgm.logcompprob( self.trace[0,ii],self.trace[1,ii]))
            return(res)
        else:
            # data across multiple dme frames
            newstart = start
            newend = min(self.dmeEndFrame[firstdmewin], end)
            # print >>sys.stderr, "dmeprobs", start, end, "not in one dme so computing", newstart, newend, "and then", newend, end
            res = []
            for ii in range(newstart,newend):
                res.append( currentgm.logcompprob( self.trace[0,ii],self.trace[1,ii]))
            # recurse: TODO remove recursion
            rest = self.dmeprobs(newend, end)
            res.extend(rest)
            return( res )
                
    ################################
    def computehmm( self, hmmstates, hmmStateTrans, dmeprobs, doGlobal=True, startframe=0 ):

        # dmeprobs are the probs over bases from startframe to endframe computed before coming in: 
        #  outprob = dme.logcompprob(self.trace[0,frame+startframe], self.trace[1,frame+startframe])

        baseIdx = {"-":0,"A":1,"C":2,"G":3,"T":4}

        # fill out simple hmm dp table
        hmmdp = np.zeros( ( len(dmeprobs), len(hmmstates) ) )-9.9E99 # log(0)
        hmmdt = np.zeros( ( len(dmeprobs), len(hmmstates) ), dtype=np.int_ ) -1
        #(1024, 71)

        # start in the first.. todo somewhat arbitrary
        if doGlobal:
            hmmdp[0,0] = dmeprobs[0][baseIdx[hmmstates[0]]]
        else:
            # local: start in any state with prob deletion=5%
            pp = math.log(0.95)
            qq = math.log(1-pp)
            for st in range(len(hmmstates)):
                # 0.95*0.05^skips*output(0th from state)
                hmmdp[0,st] = pp + st*qq + dmeprobs[0][baseIdx[hmmstates[st]]]
                hmmdt[0,st] = -1 # boundary

        # simply fill out table viterbi max
        for frame in range(1,len(dmeprobs)):
            outprob = dmeprobs[frame]
            for st in range(len(hmmstates)):
                if st>0:
                    prev = st-1
                else:
                    prev = 0
                ide = hmmstates[st]
                idprev = hmmstates[prev]
                probsame = ( hmmdp[frame-1, st]   + math.log(hmmStateTrans[ide+ide])  + outprob[baseIdx[ide]] ) # same state
                probprev = ( hmmdp[frame-1, prev] + math.log(hmmStateTrans[ide+idprev]) + outprob[baseIdx[ide]] ) # from previous
                if probsame > probprev:
                    hmmdp[frame,st] = probsame
                    hmmdt[frame,st] = st
                else:
                    hmmdp[frame,st] = probprev
                    hmmdt[frame,st] = prev

        # traceback
        traceback = []
        tracebackscore = []
        here = len(dmeprobs)-1
        if doGlobal:
            best = len(hmmstates)-1 # global must end at the last state
        else:
            # local find the best score out of all
            best = hmmdp[here,:].argmax(0)
        prev = hmmdt[here,best]
        traceback.append(int(best))
        tracebackscore.append(hmmdp[here,best])
        while here>0:
            here -= 1
            best = prev
            prev = hmmdt[here,best]
            traceback.append(int(best))
            tracebackscore.append(hmmdp[here,best])
        traceback.reverse()
        tracebackscore.reverse()

        # print "traceback"
        # for (k,v) in enumerate(zip(traceback,tracebackscore)):
        #   print k, v

        ## reduce the traceback to startframe and pw
        newbasecalls = []
        sf = startframe
        currentpos = traceback[0]
        currentbase = hmmstates[currentpos]
        ii = 0
        while ii<len(dmeprobs):
            if traceback[ii] != currentpos:
                if hmmstates[traceback[ii]] == "-": # if current is "-" then previous is base so store
                    newbasecalls.append( [currentbase, sf, startframe + ii - sf] ) # base, startframe, pulsewidth
                sf = startframe+ii
                currentpos+=1
                currentbase = hmmstates[currentpos]
            ii+=1
        # end case
        if hmmstates[currentpos] != "-":
            newbasecalls.append( [currentbase, sf, startframe + ii - sf] ) # base, startframe, pulsewidth

        # base, sf, pw
        # ('G', 49154, 15)
        # print "newbasecalls"
        # for kk in newbasecalls:
        #     print kk

        return( traceback,tracebackscore, newbasecalls )

    ################################
    def computehmmfullconnect( self, hmmstates, hmmStateTrans, dmeprobs, startframe=0 ):

        def sortindex(seq):
            return [x for x,y in sorted(enumerate(seq), key = lambda x: -x[1])]

        # dmeprobs are the probs over bases from startframe to endframe computed before coming in: 
        #  outprob = dme.logcompprob(self.trace[0,frame+startframe], self.trace[1,frame+startframe])

        baseIdx = {"-":0,"A":1,"C":2,"G":3,"T":4, "=":0,"B":1,"D":2,"H":3,"U":4}

        # fill out simple hmm dp table
        hmmdp = np.zeros( ( len(dmeprobs), len(hmmstates) ) )-9.9E99 # log(0)
        hmmdt = np.zeros( ( len(dmeprobs), len(hmmstates) ), dtype=np.int_ ) -1
        #(1024, 71)

        # start 0th frame in each of the states
        for st in range(len(hmmstates)):
            hmmdp[0,st] = dmeprobs[0][baseIdx[hmmstates[st]]]
            hmmdt[0,st] = st

        # simply fill out table viterbi max with len(hmmstates) (5) possibilities
        for frame in range(1,len(dmeprobs)):
            outprob = dmeprobs[frame]

            for st in range(len(hmmstates)):
                allprobs =[0]*len(hmmstates)
                ide = hmmstates[st]

                for prevst in range(len(hmmstates)):
                    idprev = hmmstates[prevst]
                    if ide+idprev not in hmmStateTrans:
                        allprobs[prevst] = -9.99E999
                    else:
                        allprobs[prevst] = hmmdp[frame-1, prevst] + (hmmStateTrans[ide+idprev]) + outprob[baseIdx[hmmstates[st]]]

                # find max
                si = sortindex(allprobs)
                hmmdp[frame,st] = allprobs[si[0]]
                hmmdt[frame,st] = si[0]

        # traceback
        traceback = []
        tracebackscore = []
        here = len(dmeprobs)-1
        # find the best score out of all
        best = hmmdp[here,:].argmax(0)

        prev = hmmdt[here,best]
        traceback.append(int(best))
        tracebackscore.append(hmmdp[here,best])
        while here>0:
            here -= 1
            best = prev
            prev = hmmdt[here,best]
            traceback.append(int(best))
            tracebackscore.append(hmmdp[here,best])
        traceback.reverse()
        tracebackscore.reverse()

        #### dump the data
        if False:
            np.savetxt("computehmmfullconnect.gmp.numpy", dmeprobs)
            np.savetxt("computehmmfullconnect.hmmdp.numpy", hmmdp)
            np.savetxt("computehmmfullconnect.hmmdt.numpy", hmmdt)
            ofp = open("computehmmfullconnect.traceback","w")
            for (k,v) in enumerate(zip(traceback,tracebackscore)):
                print >>ofp, k, v[0], v[1]
            ofp.close()

        ## reduce the traceback to startframe and pw
        newbasecalls = []
        sf = startframe
        currentpos = traceback[0]
        currentscore=  tracebackscore[0]
        currentbase = hmmstates[currentpos]
        ii = 0
        while ii<len(dmeprobs):
            if traceback[ii] != currentpos:
                newbasecalls.append( [currentbase, sf, startframe + ii - sf, currentscore] ) # base, startframe, pulsewidth, score
                sf = startframe+ii
                currentpos = traceback[ii]
                currentscore = tracebackscore[ii]
                currentbase = hmmstates[currentpos]
            ii+=1
        # end case
        newbasecalls.append( [currentbase, sf, startframe + ii - sf, currentscore] ) # base, startframe, pulsewidth, score


        # base, sf, pw
        # ('G', 49154, 15)
        # print "newbasecalls"
        # for kk in newbasecalls:
        #     print kk

        return( traceback,tracebackscore, newbasecalls )

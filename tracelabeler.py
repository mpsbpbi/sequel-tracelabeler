import h5py
import numpy as np
import bisect
import sys
import pysam
import pbcore.io
import math

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
    def __init__( self, trace="tracefile", dme="dmefile", unbam="unalignedbamfile", albam="alignedbamfile", ref="reffile"):
        self.tracef=trace
        self.dmef=dme
        self.unbamf=unbam
        self.albamf=albam
        self.reff=ref

    ################################
    def setzmw( self, zmw, targetbase=0 ):
        "Set the zmw subread that contains targetbase (bam subreads) and pull all the data"

        self.zmw = zmw
        self.targetbase = targetbase

        #### trace data
        self.traceff=h5py.File(self.tracef,"r")
        self.tridx = bisect.bisect(self.traceff['TraceData']['HoleNumber'],self.zmw)-1
        assert(self.zmw==self.traceff['TraceData']['HoleNumber'][self.tridx])
        self.trace = self.traceff['TraceData']['Traces'][self.tridx]
        print "trace got data"

        #### dme
        self.dmeff=h5py.File(self.dmef,"r")
        self.dmeidx = bisect.bisect(self.dmeff['HoleNumber'],self.zmw)-1
        assert(self.zmw==self.dmeff['HoleNumber'][self.dmeidx])
        print "dme got data"

        #### unbam
        self.unbam = None
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
                    break
        if not self.unbam: print "ERROR: unbam not found!!!"

        #### albam
        self.albam = []
        # bam doesn't have the reference
        reader = pbcore.io.openIndexedAlignmentFile( self.albamf, self.reff )
        for h in reader:
            zmw=h.HoleNumber
            # TODO: do I really have to parse the subread location from the query_name???
            (substart, subend) = [int(xx) for xx in self.unbam.query_name.split("/")[-1].split("_")]
            if zmw == self.zmw and AContainedInB(h.queryStart, h.queryEnd, substart, subend):
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

        #### ref
        myref = []
        dat = open(self.reff).read().splitlines()
        for ii in range(1,len(dat)):
            myref.append(dat[ii])
        self.ref = "".join(myref)
        
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
        self.alignReadStart = qs - h.queryStart # subtract off start of subread
        self.alignReadEnd = qe - h.queryStart # subtract off start of subread
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
    def baseToTrace(self):
        "for each base in read set (start,end) frames in trace"
        self.rsf = self.unbam.get_tag("sf") # start frame
        self.rpw = self.unbam.get_tag("pw") # pulse width
        self.rseq = self.unbam.seq
        
    ################################
    def frameWindowToReadStartEnd(self, window):
        "Given a 1kb window number return the readStart and readEnd of read bases in that 1kb window of frames"
        readStart = bisect.bisect( self.rsf, 1024*window)
        readEnd =   bisect.bisect( self.rsf, 1024*(window+1))
        return(readStart,readEnd)

    ################################
    def readBaseToWindowStartFrame(self, basenum):
        "Given read base, return window containing start frame"
        return(int(self.rsf[basenum] / 1024.0 + 0.5))

    ################################
    def rawAlignRefReadFromReadSE(self, start, end):
        "raw alignment containing read start to end"
        als = tl.readToAlign(readStart)[0]
        ale = tl.readToAlign(readEnd)[0]
        return( (tl.albam[0].reference()[als:ale ], tl.albam[0].read()[als:ale ]))

    ################################
    def subreadDat(self, start, end):
        "for subread start to end return the relevant data"
        result = {}
        result["dat"] = []
        result["doc"] = "readbaseindex startframe pulsewidth refindex alignindex readbase alignrefbases"
        for readindex in range(start,end):
            resrseq = self.rseq[readindex]
            #resreadToRef = "".join( [self.ref[xx] if self.alignStrand==1 else self.rMap[self.ref[xx]] for xx in self.readToRef(readindex)])
            resreadToAlign_ref = "".join( [self.alignRef[xx] for xx in self.readToAlign(readindex)])
            #resreadToAlign_read = "".join( [self.alignRead[xx] for xx in self.readToAlign(readindex)])
            result["dat"].append([ readindex, self.rsf[readindex], self.rpw[readindex], self.readToRef(readindex), self.readToAlign(readindex), resrseq, resreadToAlign_ref ])
        return(result)
        
    ################################
    def windowToDMEDat( self, window):
        # We_order_them_as:_[-,_A,_C,_G,_T]

        # TODO: hack! 0th window is always 0????
        if window==0: window =1

        MixtureFraction = self.dmeff["/SmoothedEstimates/MixtureFraction"][window][self.dmeidx]
        Mean = self.dmeff["/SmoothedEstimates/Mean"][window][self.dmeidx]
        Covariance = self.dmeff["/SmoothedEstimates/Covariance"][window][self.dmeidx]
        BaselineMean = self.dmeff["/BaselineMean"][window][self.dmeidx]
        return({"MixtureFraction": MixtureFraction, "Mean": Mean, "Covariance": Covariance, "BaselineMean": BaselineMean})

    ################################
    def computehmm( self, hmmstates, hmmStateTrans, dme, startframe, endframe ):

        baseIdx = {"-":0,"A":1,"C":2,"G":3,"T":4}

        # fill out simple hmm dp table
        hmmdp = np.zeros( ( (endframe-startframe), len(hmmstates) ) )-9.9E99 # log(0)
        hmmdt = np.zeros( ( (endframe-startframe), len(hmmstates) ), dtype=np.int_ ) -1
        #(1024, 71)

        # start in the first.. todo somewhat arbitrary
        frame = 0
        outprob = dme.logcompprob(self.trace[0,frame+startframe], self.trace[1,frame+startframe])
        hmmdp[0,0] = outprob[baseIdx[hmmstates[0]]]

        # simply fill out table viterbi max
        for frame in range(1,(endframe-startframe)):
            outprob = dme.logcompprob(self.trace[0,frame+startframe], self.trace[1,frame+startframe])
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
        here = (endframe-startframe)-1
        best = len(hmmstates)-1 # global must end at the last state
        # "local" find the best score out of all: best = hmmdp[here,:].argmax(0)
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
        currentpos = 0
        currentbase = hmmstates[currentpos]
        ii = 0
        while ii<(endframe-startframe):
            if traceback[ii] != currentpos:
                if hmmstates[traceback[ii]] == "-": # if current is "-" then previous is base so store
                    newbasecalls.append( (currentbase, sf, startframe + ii - sf) ) # base, startframe, pulsewidth
                sf = startframe+ii
                currentpos+=1
                currentbase = hmmstates[currentpos]
            ii+=1
        # TODO: is last case handled correcly. need to add below?

        # base, sf, pw
        # ('G', 49154, 15)
        # print "newbasecalls"
        # for kk in newbasecalls:
        #     print kk

        return( traceback,tracebackscore, newbasecalls )

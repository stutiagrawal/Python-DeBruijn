import numpy as np
import matplotlib.pyplot as plt
from val import *
import sys
class graph(object):
    #Initialize the graph class
    def __init__(self, k, readsfile):
        self.kMerLen = k
        self.kMers = dict()
        self.reads = readsfile
        self.graph = dict()
    
    #Get the k-mers also output a sorted list of kmers
    def createKmers(self):
        if(is_valid_file(self.reads)):
            file = open(self.reads, "r")
            for line in file:
                if(not(line.startswith(">"))):
                    for i in xrange(0, len(line)-self.kMerLen):    
                        kMer = self.str2kmer(line[i:i+self.kMerLen])
                        if(kMer not in self.kMers):
                            self.kMers[kMer]=1
                        else:
                            self.kMers[kMer]+=1
                        
            file.close()
        return None

    #Get the prefix
    def getPrefix(self, kMer):
        return kMer>>2
    
    def str2kmer(self,s):
        result=0
        for c in s:
            if c=='A':
                result+=0
            elif c=='T':
                result+=1
            elif c=='G':
                result+=2
            elif c=='C':
                result+=3
            result=result<<2
        return result>>2
    
    def kmer2str(self,kmer):
        result=''
        while(len(kmer)>0):
            c=kmer&3
            if c==0:
                result='A'+result
            elif c==1:
                result='T'+result
            elif c==2:
                result='G'+result
            elif c==3:
                result='C'+result
            kmer=kmer>>2
        return result

    def kmer2str(kmer,count):
        result=''
        while(count>0):
            count-=1
            c=kmer&3
            if c==0:
                result='A'+result
            elif c==1:
                result='T'+result
            elif c==2:
                result='G'+result
            elif c==3:
                result='C'+result
            if kmer>0:
                kmer=kmer>>2
            
        return result

    def getSuffix(self,kMer):
        return kMer&(2**(self.K*2-3)-1)

    def createInterGraph(self):
        kMerList=sorted(self.kMers.keys())
        interGraph=dict()
        for kMer in kMerList:
            prefix = self.getPrefix(kMer)

        #may not need to check###############################################################
            if prefix not in self.graph:
                interGraph[prefix] = list()
            interGraph[prefix].append(kMer[-1])
        return interGraph
        #kMerList = self.createKmers(self)
  
    def createGraph(self):
        self.createKmers()
        print 'kmers complete'
        interGraph=self.createInterGraph()
        print 'intermediate graph complete'
        self.graph=dict()
        kMerList=self.kMers.keys()
        for kMer in kMerList:
            # print kMer,interGraph.get(self.getSuffix(kMer))
            self.graph[kMer]=interGraph.get(self.getSuffix(kMer),[])
        return None

    def plotKmerCoverage(self,filename):
        if len(self.kMers)==0:
            self.createkMers()
        else :
            values=self.kMers.values()
            plt.plot(values)
            plt.savefig(filename)

    def kMerStats(self):
        if len(self.kMers)==0:
            self.createKmers()
        values=self.kMers.values()
        return (np.mean(values),np.var(values))

    def connectivity(self):
        if len(self.graph)==0:
            self.createGraph()
        degrees=[]
        for value in self.graph.values():
            degrees.append(len(value))

        return np.mean(degrees)

    def size(self):
        return len(self.graph.keys()), len(self.graph.values())

    # def selfLoop(self):


###########################################################################
# a = graph(25, sys.argv[1])
# print a.reads
# print a.kMerLen
# #print a.createKmers()
# print a.getPrefix("ATTA")
# print a.createGraph(["ATTA", "ATTT"])

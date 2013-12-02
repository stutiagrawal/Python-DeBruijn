from val import *
import sys
class graph(object):
    #Initialize the graph class
    def __init__(self, k, readsfile):
        self.kMerLen = k
        self.kMers = set()
        self.reads = readsfile
        self.graph = dict()
    
    #Get the k-mers
    def createKmers(self):
        if(is_valid_file(self.reads)):
            file = open(self.reads, "r")
            for line in file:
                if(not(line.startswith(">"))):
                    for i in xrange(0, len(line)-self.kMerLen):    
                        kMer = line[i:i+self.kMerLen]
                        if(kMer not in self.kMers):
                            self.kMers.add(kMer)
            file.close()
        return sorted(list(self.kMers))

    #Get the prefix
    def getPrefix(self, kMer):
        return kMer[:-1]

    def createGraph(self, kMerList):
        for kMer in kMerList:
            prefix = self.getPrefix(kMer)
            if prefix not in self.graph:
                self.graph[prefix] = list()
            self.graph[prefix].append(kMer[-1])
        return self.graph
        #kMerList = self.createKmers(self)

###########################################################################
a = graph(25, sys.argv[1])
print a.reads
print a.kMerLen
#print a.createKmers()
print a.getPrefix("ATTA")
print a.createGraph(["ATTA", "ATTT"])


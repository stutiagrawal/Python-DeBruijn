def str2kmer(s):
	return s



def reads2kmers(reads,k):
	kmers=Set([])
	for read in reads:
		for i in xrange(len(read)):
			kmer=str2kmer(reads[i,i+k])
			kmers.add(kemer)
	return sorted(list(kmers))

def 
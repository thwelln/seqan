#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>	

using namespace seqan;

int main(int argc, char *argv[])
{
	// IMPORTANT VARIABLES
	
		unsigned readStartPos = atoi(argv[1]);
		unsigned readLength = atoi(argv[2]);
		double readErrorRate = atof(argv[3]);
	
	typedef GappedShape<HardwiredShape<2,2,3,2> > TInsideShape;	
	
	unsigned klen = 10; // length of k-mere devision in pattern
	
	//READ FILES IN
	
	CharString seqFileName = getAbsolutePath("/../Sequences/sequence.fasta");
	Dna5String seq;
	Dna5String read;
	CharString id;
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seq, seqFileIn);
    
    char outpath [256];
	sprintf(outpath, "/../Sequences/Reads/read_%d_%d_%.2f.fasta",readStartPos,readLength,readErrorRate);
	CharString readFileName = getAbsolutePath(outpath);  
    SeqFileIn readFileIn(toCString(readFileName));
    readRecord(id, read, readFileIn);    
    
	
	std::cout << "READYYYY!" << std::endl;
	
	// BUILDING INDEX
	
	Index<Dna5String, IndexQGram<TInsideShape> > qgindex(seq);
    //stringToShape(indexShape(qgindex), "111");
	
	//SEARCHING
	unsigned tp = 0;
	unsigned fn = 0;
	unsigned fp = 0;
	
	for (unsigned ki=0; ki<(length(read)/klen);++ki)
	{
		unsigned compareStartpos = ki*klen;
		Dna5String qgram = infixWithLength(read, (compareStartpos), klen);
		hash(indexShape(qgindex), begin(qgram));
		bool found = 0;
		for (unsigned i=0; i<length(getOccurrences(qgindex, indexShape(qgindex))); ++i)
		{
			unsigned findPos = getOccurrences(qgindex, indexShape(qgindex))[i];
			
			//std::cout << ki << " : " << findPos << "\t";
			if (findPos == readStartPos+ki*klen)
			{
				found = 1;
			}
			else
			{
				++fp;
			}
		}
		if (found)
		{
			++tp;
		}
		else
		{
			++fn;
		}
		//std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "TP:	" << tp << std::endl;
	std::cout << "FN:	" << fn << std::endl;
	std::cout << "FP:	" << fp << std::endl;	
		
	return 0;
}

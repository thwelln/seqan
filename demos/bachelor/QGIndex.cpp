#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>

#include <time.h>

using namespace seqan;

typedef GappedShape<HardwiredShape<1> > TInsideShape;	


int main(int argc, char *argv[])
{
	// IMPORTANT VARIABLES
	
		unsigned readStartPos = atoi(argv[1]);
		unsigned readLength = atoi(argv[2]);
		double readErrorRate = atof(argv[3]);
	
	
	unsigned klen = 3; // length of k-mere devision in pattern
	
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
    

	
	// BUILDING INDEX
	std::cout << 0 << std::endl;
	double tim = sysTime();
	Index<Dna5String, IndexQGram<TInsideShape> > qgindex(seq);
    //stringToShape(indexShape(qgindex), "111");
    hash(indexShape(qgindex),begin("AAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    std::cout  << sysTime() - tim << std::endl;
	
	//SEARCHING
	unsigned tp = 0;
	unsigned fn = 0;
	unsigned long fp = 0;
	tim = sysTime();
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
    std::cout  << sysTime() - tim << std::endl;	
	std::cout << tp << std::endl;
	std::cout << fn << std::endl;
	std::cout << fp << std::endl;	
		
	return 0;
}

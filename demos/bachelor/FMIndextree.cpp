#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>	

using namespace seqan;

int main(int argc, char *argv[])
{
		unsigned readStartPos = atoi(argv[1]);
		unsigned readLength = atoi(argv[2]);
		double readErrorRate = atof(argv[3]);
	
	unsigned klen = 10; // length of k-mere devision in pattern
	
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
    

	reverse(read); //ONLY FOR FM
	
	std::cout << "READYYYY!" << std::endl;

// BUILDING INDEX

	
	typedef Index<Dna5String, FMIndex<> > TFMIndex;
	TFMIndex fmindex(seq);
	
	Iterator<TFMIndex, TopDown<> >::Type fmit(fmindex);
	
// SEARCHING
	
	unsigned tp = 0;
	unsigned fn = 0;
	unsigned fp = 0;
	
	for (unsigned ki=0; ki<(length(read)/klen);++ki)
	{
		unsigned compareStartpos = ki*klen;
		if (!goDown(fmit, infixWithLength(read, (compareStartpos), klen))) // compare full k-mere
		{
			std::cout << "FN!" << std::endl;
			goRoot(fmit);
			++fn;
			continue;
		}
		unsigned compareLength = klen;
		while (countOccurrences(fmit)>1 && compareStartpos+compareLength<length(read))
		{
			if (!goDown(fmit, read[compareStartpos+compareLength])) // compare next letter
			{
				//std::cout << "ERROR!" <<std::endl;
				break;
			}
		++compareLength;
		}
		bool found = 0;
		for (unsigned i=0; i<countOccurrences(fmit); ++i)
		{
			unsigned findPos = getOccurrences(fmit)[i]+compareLength;
			
			std::cout << ki << " : " << findPos << "\t";
			if (findPos == readStartPos+length(read)-ki*klen)
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
		goRoot(fmit);
		std::cout << std::endl;
	}
	
// OUTPUT
	std::cout << std::endl;
	std::cout << "TP:	" << tp << std::endl;
	std::cout << "FN:	" << fn << std::endl;
	std::cout << "FP:	" << fp << std::endl;	
		
	return 0;
}

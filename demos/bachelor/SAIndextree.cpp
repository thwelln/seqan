#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>	

using namespace seqan;

int main()
{
	unsigned readStartPos=270000;
	unsigned readLength=150;
	unsigned readErrorRate=20;
	
	unsigned klen = 15; // length of k-mere devision in pattern
	
	CharString seqFileName = getAbsolutePath("/../Sequences/sequence.fasta");
	Dna5String seqin;
	Dna5String read;
	CharString id;
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seqin, seqFileIn);
    
    char outpath [256];
	sprintf(outpath, "/../Sequences/Reads/read_%d_%d_%d.fasta",readStartPos,readLength,readErrorRate);
	CharString readFileName = getAbsolutePath(outpath);  
    SeqFileIn readFileIn(toCString(readFileName));
    readRecord(id, read, readFileIn);    
    
	//DnaString text = "AGTATCTCATTGACTTAACG";
	//DnaString pattern = "TTACGTC";
	unsigned partialLength = 200000;
	Dna5String seq = infixWithLength(seqin, readStartPos-(partialLength/2), partialLength);
	//reverse(read);
	
	std::cout << "READYYYY!" << std::endl;
	
	
	typedef Index<Dna5String, IndexSa<> > TSAIndex;
	TSAIndex saindex(seq);
	
	Iterator<TSAIndex, TopDown<> >::Type sait(saindex);
	
	unsigned tp = 0;
	unsigned fn = 0;
	unsigned fp = 0;
	
	for (unsigned ki=0; ki<(length(read)/klen);++ki)
	{
		unsigned compareStartpos = ki*klen;
		if (!goDown(sait, infixWithLength(read, (compareStartpos), klen))) // compare full k-mere
		{
			std::cout << "FN!" << std::endl;
			goRoot(sait);
			++fn;
			continue;
		}
		unsigned compareLength = klen;
		while (countOccurrences(sait)>1 && compareStartpos+compareLength<length(read))
		{
			if (!goDown(sait, read[compareStartpos+compareLength])) // compare next letter
			{
				//std::cout << "ERROR!" <<std::endl;
				break;
			}
		++compareLength;
		}
		bool found = 0;
		for (unsigned i=0; i<countOccurrences(sait); ++i)
		{
			unsigned findPos = getOccurrences(sait)[i];
			
			std::cout << ki << " : " << findPos << "\t";
			if (findPos == (partialLength/2)+ki*klen)
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
		goRoot(sait);
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout << "TP:	" << tp << std::endl;
	std::cout << "FN:	" << fn << std::endl;
	std::cout << "FP:	" << fp << std::endl;	
		
	return 0;
}

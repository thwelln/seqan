#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>	

using namespace seqan;


template <typename TIterator>
int backtrack(const Dna5String & read, TIterator & it, unsigned errors, const unsigned & windowStart, unsigned & maxLevel, String<unsigned> & occ)
{
	unsigned compareLevel = nodeDepth(it);
	//std::cout << "CALL Level" << compareLevel << "ERRORS " << errors << std::endl;
	if(windowStart+compareLevel>=length(read))
	{
		return -1;
	}
	if (errors == 0)
	{
		if (!goDown(it, read[windowStart+compareLevel]))
		{
			if (maxLevel < compareLevel)
			{
				clear(occ);
				appendValue(occ, 5);
				//append(occ,getOccurrences(it));
				maxLevel = compareLevel;
				//std::cout << "NEW MAX " << maxLevel << " ";
				//printUnsignedString(occ);
			}
			else if (maxLevel == compareLevel)
			{
				//append(occ,getOccurrences(it));
				//std::cout << "MAX REACHED " << maxLevel << " ";
				//printUnsignedString(occ);		
			}
			return 0;
		}
		else
		{
			backtrack(read, it, 0, windowStart, maxLevel, occ);
			goUp(it);
			return 0;
		}
	}
	else
	{
		Dna5String all = "ACGTN";
		for (unsigned i = 0; i < length(all); i++)
		{
			//std::cout << "!!!" << std::endl;
			TIterator it2 = it;
			unsigned errors2 = errors;
			if (read[windowStart+compareLevel] != all[i])
			{
				--errors2;
			}
			if (!goDown(it2, all[i]))
			{
				if (maxLevel < compareLevel)
				{
					clear(occ);
					//append(occ,getOccurrences(it));
					maxLevel = compareLevel;
					//std::cout << "NEW MAX " << maxLevel << " ";
					//printUnsignedString(occ);
					
					
				}
				else if (maxLevel == compareLevel)
				{
					//append(occ,getOccurrences(it));
					//std::cout << "MAX REACHED " << maxLevel << " ";
					//printUnsignedString(occ);			
				}
			}
			else
			{
				backtrack(read, it2, errors2, windowStart, maxLevel, occ);
			}			
		}
		goUp(it);
	}
	return 0;
}

int main(int argc, char *argv[])
{
		unsigned readStartPos = atoi(argv[1]);
		unsigned readLength = atoi(argv[2]);
		double readErrorRate = atof(argv[3]);
		unsigned allowedErrors = atoi(argv[4]);
	
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
	
	Iterator<TFMIndex, TopDown<ParentLinks<> > >::Type fmit(fmindex);
	
// SEARCHING
	
	unsigned tp = 0;
	unsigned fn = 0;
	unsigned fp = 0;
	
	for (unsigned ki=0; ki<(length(read)/klen);++ki)
	{
		unsigned compareStartpos = ki*klen;

		
		String <unsigned> occs;
		unsigned maxL = 0;
	
		backtrack(read, fmit, allowedErrors, compareStartpos, maxL, occs);

		if (maxL < klen)
		{
			fn++;
		}
		else
		{
			bool found = 0;
			for (unsigned i=0; i<length(occs); ++i)
			{
				unsigned findPos = occs[i]+maxL;
				
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
		}
		goRoot(fmit);
		//std::cout << std::endl;
	}
	
	// OUTPUT
	std::cout << std::endl;
	std::cout << "TP:	" << tp << std::endl;
	std::cout << "FN:	" << fn << std::endl;
	std::cout << "FP:	" << fp << std::endl;	
		
	return 0;
}

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>	

using namespace seqan;

void printUnsignedString(String<unsigned> st)
{
	std::cout << "[";
	for (unsigned i=0; i<length(st); i++)
	{
		std::cout << st[i] << ", ";
	}
	std::cout << "]" << std::endl;
}

template <typename TIterator>
int backtrack(const Dna5String & read, TIterator & it, unsigned errors, const unsigned & windowStart, unsigned & maxLevel, String<unsigned> & occ)
{
	unsigned compareLevel = nodeDepth(it);
	std::cout << "CALL Level" << compareLevel << "ERRORS " << errors << std::endl;
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
				append(occ,getOccurrences(it));
				maxLevel = compareLevel;
				//std::cout << "NEW MAX " << maxLevel << " ";
				//printUnsignedString(occ);
			}
			else if (maxLevel == compareLevel)
			{
				append(occ,getOccurrences(it));
				//std::cout << "MAX REACHED " << maxLevel << " ";
				//printUnsignedString(occ);		
			}
			return 0;
		}
		else
		{
			backtrack(read, it, 0, windowStart, (compareLevel+1), maxLevel, occ);
			goUp(it);
			return 0;
		}
	}
	else
	{
		Dna5String all = "ACGTN";
		for (unsigned i = 0; i < length(all); i++)
		{
			std::cout << "!!!" << std::endl;
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
					append(occ,getOccurrences(it));
					maxLevel = compareLevel;
					//std::cout << "NEW MAX " << maxLevel << " ";
					//printUnsignedString(occ);
					
					
				}
				else if (maxLevel == compareLevel)
				{
					append(occ,getOccurrences(it));
					//std::cout << "MAX REACHED " << maxLevel << " ";
					printUnsignedString(occ);			
				}
			}
			else
			{
				backtrack(read, it2, errors2, windowStart, (compareLevel+1), maxLevel, occ);
			}			
		}
		goUp(it);
	}
	return 0;
}

int main(int argc, char *argv[])
{
	Dna5String seq = "ACGAACTATAG";
	Dna5String read = "ANNNNTAGTACG";

	typedef Index<Dna5String, IndexSa<> > TSAIndex;	
	//typedef Index<Dna5String, FMIndex<> > TSAIndex;
	TSAIndex saindex(seq);
	
	Iterator<TSAIndex, TopDown<ParentLinks<>> >::Type sait(saindex);
	
	String <unsigned> occs;
	unsigned maxL = 0;
	
	backtrack(read, sait, 2, 0, maxL, occs);
	printUnsignedString(occs);
	
	return 0;
}

// ==========================================================================
//                         benchmark_sa_construction
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Thimo Wellner <thimo.wellner@fu-berlin.de> (Modified from Sascha Meiers)
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/pipemyown.h>
#include <seqan/indexmyown.h>
#include <omp.h>


#include <time.h>

using namespace seqan;


template <typename TIterator>
int backtrack(const String<unsigned> & read, TIterator it, unsigned errors, const unsigned & windowStart, unsigned compareLevel, unsigned & maxLevel, String<unsigned> & occ, const unsigned & sigma)
{
	//std::cout << "CALL Level" << compareLevel << "ERRORS " << errors << std::endl;
	if(windowStart+compareLevel>=length(read))
	{
		return -1;
	}
	if (errors == 0)
	{
		if (!goDown(it, read[windowStart+compareLevel]))
		{		
			#pragma omp critical
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
			}
		return 0;
		}
		else
		{
			backtrack(read, it, 0, windowStart, (compareLevel+1), maxLevel, occ, sigma);
			return 0;
		}
	}
	else
	{
		#pragma omp parallel
		#pragma omp for
		for (unsigned i = 0; i < sigma; i++)
		{
			//std::cout << "!!!" << std::endl;
			TIterator it2 = it;
			unsigned errors2 = errors;
			if (read[windowStart+compareLevel] != i)
			{
				--errors2;
			}
			if (!goDown(it2, i))
			{
				#pragma omp critical
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
				}
			}
			else
			{
				backtrack(read, it2, errors2, windowStart, (compareLevel+1), maxLevel, occ, sigma);
			}			
		}
	}
	return 0;
}

void printUnsignedString(String<unsigned> st)
{
	std::cout << "[";
	for (unsigned i=0; i<length(st); i++)
	{
		std::cout << st[i] << ", ";
	}
	std::cout << "]" << std::endl;
}

// ==========================================================================
// Functions
// ==========================================================================



    template <
        typename TSA,
        typename TString,
        typename TSpec,
        typename TAlgSpec >
    void _makeDislexString(
        String <unsigned> &dis,
        unsigned & sig,
        TSA &suffixArray,
        StringSet<TString, TSpec> const &stringSet,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename Concatenator<StringSet<TString, TSpec> >::Type            TConcat;
        typedef typename MakeUnsigned_< typename Value<TConcat>::Type >::Type    TUValue;
        typedef typename StringSetLimits<StringSet<TString, TSpec> >::Type TLimits;
        typedef Multi<
            TAlgSpec,
            typename Value<TSA>::Type,
            TLimits >        				MultiConstrSpec;

        // specialization
        typedef Pipe< TConcat, Source<> >                src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
        typedef Pipe< unsigner_t, MultiConstrSpec >        creator_t;

        // instantiation and processing
        
        src_t        src(concat(stringSet));
        unsigner_t  unsigner(src);
        creator_t    creator(unsigner, stringSetLimits(stringSet));
		dis = creator.dislexString_ordered;
		sig = creator.sigma;
		
        //suffixArray << creator;

		
    }


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.


using namespace seqan;

int main(int argc, char *argv[])
{
		//IMPORTANT VARIABLES
		unsigned readStartPos = atoi(argv[1]);
		unsigned readLength = atoi(argv[2]);
		double readErrorRate = atof(argv[3]);
		unsigned allowedErrors = atoi(argv[4]);
	
	unsigned klen = 10; // length of k-mere devision in pattern
	
	typedef GappedShape<HardwiredShape<2,2,3,2> > TInsideShape;	
	typedef CyclicShape<FixedShape<0,TInsideShape, 0> > TShape;
	
	// 
	
	CharString seqFileName = getAbsolutePath("/../Sequences/sequence.fasta");
	Dna5String seqIn;
	Dna5String readIn;
	CharString id;

    
    char outpath [256];
	sprintf(outpath, "/../Sequences/Reads/read_%d_%d_%.2f.fasta",readStartPos,readLength,readErrorRate);
	CharString readFileName = getAbsolutePath(outpath);
    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seqIn, seqFileIn);
    SeqFileIn readFileIn(toCString(readFileName));
    readRecord(id, readIn, readFileIn);
    
    StringSet<seqan::Dna5String> seqs;
    appendValue(seqs, seqIn);
    appendValue(seqs, readIn);
	std::cout << "READ!" << std::endl <<std::endl;
	
		// TRANSLATING INTO DISLEX

            typedef Index<StringSet<Dna5String> > TIndex;
			typedef typename StringSetLimits<StringSet<Dna5String> >::Type TLimits;
			typedef Pair<unsigned, unsigned> TUPair;			
            TIndex index(seqs);
            String <unsigned> dislex;
			unsigned sigma;
			_makeDislexString(dislex, sigma, indexSA(index), indexText(index), DislexExternal<TShape, Skew7>());
			

			TLimits lim = stringSetLimits(seqs);
			//printUnsignedString(dislex);
			
			//std::cout << length(dislex) << std::endl;
			
			
			String <unsigned> seq = prefix(dislex, posGlobalize(TUPair(1,0),lim));
			String <unsigned> read = suffix(dislex, posGlobalize(TUPair(1,0),lim));


			//printUnsignedString(seq);
			//printUnsignedString(read);
			//std::cout << length(seq) << std::endl;
			//std::cout << length(read) << std::endl;
			
				// BUILDING INDEX
				typedef Index<String<unsigned>, IndexSa<> > TSAIndex;
				TSAIndex saindex(seq);
				//std::cout << "DONE!";
				Iterator<TSAIndex, TopDown<> >::Type sait(saindex);
				//std::cout << "DONE!";				
				//SEARCHING
				unsigned tp = 0;
				unsigned fn = 0;
				unsigned fp = 0;
				
	
				for (unsigned ki=0; ki<(length(read)/klen);++ki)
				{
					unsigned compareStartpos = ki*klen;
					
					String <unsigned> occs;
					unsigned maxL = 0;
				
					backtrack(read, sait, allowedErrors, compareStartpos, 0, maxL, occs,sigma);
			
					if (maxL == 0)
					{
						fn++;
					}
					else
					{
						bool found = 0;
						for (unsigned i=0; i<length(occs); ++i)
						{
							unsigned findPos = occs[i];
							
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
					}
					goRoot(sait);
					//std::cout << std::endl;
				}
				//OUTPUT
				std::cout << std::endl;
				std::cout << "TP:	" << tp << std::endl;
				std::cout << "FN:	" << fn << std::endl;
				std::cout << "FP:	" << fp << std::endl;	
			
				
    return 0;
}

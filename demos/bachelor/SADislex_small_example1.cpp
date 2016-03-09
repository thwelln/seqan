// ==========================================================================
//                         benchmark_sa_construction
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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


#include <time.h>

using namespace seqan;

typedef GappedShape<HardwiredShape<1> > TInsideShape;
typedef CyclicShape<FixedShape<0,TInsideShape, 0> > TShape;

template <typename TLimits>
unsigned getRealPos (TLimits lim, bool pattern,	unsigned dislexPos)
{
	typedef GappedShape<HardwiredShape<1> > TInsideShape;
	
	typedef CyclicShape<FixedShape<0,TInsideShape, 0> > TShape;
	typedef Pair<unsigned, unsigned> TUPair;
	typedef DislexReverseTransformMulti_<unsigned,
	TLimits, TUPair >                             TGetDislexReversePos;
	
		unsigned fullDisPos = posGlobalize(TUPair(pattern,dislexPos), lim);
		//std::cout << "test " << fullDisPos << std::endl;
		TUPair pair = TGetDislexReversePos(TShape::span,lim) (fullDisPos);
		//std::cout << pair << std::endl;
		return pair.i2;
}


template <typename TLimits, typename TUPair>
unsigned getDislexPos (TLimits lim, TUPair pair)
{
	typedef GappedShape<HardwiredShape<1> > TInsideShape;
	
	typedef CyclicShape<FixedShape<0,TInsideShape, 0> > TShape;
	typedef DislexTransformMulti_<TUPair,
	TLimits>                             TGetDislexPos;
	
		unsigned zeroPos = TGetDislexPos(TShape::span,lim) (pair);
		//std::cout << fullDisPos << pair << std::endl;
		return zeroPos;
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
		dis = creator.dislexString;
		
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
	
	unsigned klen = 1; // length of k-mere devision in pattern
	
	typedef GappedShape<HardwiredShape<1> > TInsideShape;	
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
    
    Dna5String text = "ACAGGTACTA";
    Dna5String pat = "CATGTA";
    
    StringSet<seqan::Dna5String> seqs;
    appendValue(seqs, text);
    appendValue(seqs,  pat);
	
		// TRANSLATING INTO DISLEX

            typedef Index<StringSet<Dna5String> > TIndex;
			typedef typename StringSetLimits<StringSet<Dna5String> >::Type TLimits;
			typedef Pair<unsigned, unsigned> TUPair;			
            TIndex index(seqs);
            String <unsigned> dislex;
			double tim = sysTime();
			_makeDislexString(dislex, indexSA(index), indexText(index), DislexExternal<TShape, Skew7>());
			

			TLimits lim = stringSetLimits(seqs);
			//printUnsignedString(dislex);
			
			std::cout << length(dislex) << std::endl;
			
			
			String <unsigned> seq = prefix(dislex, posGlobalize(TUPair(1,0),lim));
			String <unsigned> read = suffix(dislex, posGlobalize(TUPair(1,0),lim));
			std::cout  << sysTime() - tim << std::endl;
			printUnsignedString(seq);
			printUnsignedString(read);
			std::cout << length(seq) << std::endl;
			std::cout << length(read) << std::endl;
			for (unsigned i=0; i<10; i++)
			{
				std::cout << getRealPos(lim, 0, i)<< " ";
			}
			
			std::cout << std::endl;
			for (unsigned j=0; j<6; j++)
			{
				std::cout << getRealPos(lim, 1, j)<< " ";
			}
			std::cout << std::endl;
			std::cout << getDislexPos(lim, TUPair(1,0)) << std::endl;
			
				
    return 0;
}

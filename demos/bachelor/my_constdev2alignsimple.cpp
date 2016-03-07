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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/pipemyown.h>
#include <seqan/indexmyown.h>


#include <time.h>

using namespace seqan;

template <typename TLimits>
unsigned getRealPos (TLimits lim, bool pattern,	unsigned dislexPos)
{
	typedef GappedShape<HardwiredShape<1,1> > TInsideShape;
	
	//typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
	typedef CyclicShape<FixedShape<0,TInsideShape, 1> > TShape;
	typedef Pair<unsigned, unsigned> TUPair;
	typedef DislexReverseTransformMulti_<unsigned,
	TLimits, TUPair >                             TGetDislexReversePos;
	
		unsigned fullDisPos = posGlobalize(TUPair(pattern,dislexPos), lim);
		TUPair pair = TGetDislexReversePos(TShape::span,lim) (fullDisPos);
		//std::cout << fullDisPos << pair << std::endl;
		return pair.i2;
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
        //typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
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
		
        //suffixArray << creator;

		
    }


// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

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
	
	unsigned klen = 2; // length of k-mere devision in pattern
	
	CharString seqFileName = getAbsolutePath("/../Sequences/sequence.fasta");
	Dna5String seqIn;
	Dna5String readIn;
	CharString id;

    
    char outpath [256];
	sprintf(outpath, "/../Sequences/Reads/read_%d_%d_%.2f.fasta",readStartPos,readLength,readErrorRate);
	CharString readFileName = getAbsolutePath(outpath);
    SeqFileIn seqFileIn(toCString(readFileName));
    readRecord(id, seqIn, seqFileIn);
    //SeqFileIn readFileIn(toCString(readFileName));
    //readRecord(id, readIn, readFileIn);
    
    readIn = infixWithLength(seqIn, 60, 80);
    
    StringSet<seqan::Dna5String> seqs;
    appendValue(seqs, seqIn);
    appendValue(seqs, readIn);
	std::cout << "READ!" << std::endl;

			//typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
			typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,1> >, 1> > TShape;
            typedef Index<StringSet<Dna5String> > TIndex;
			typedef typename StringSetLimits<StringSet<Dna5String> >::Type TLimits;
			typedef Pair<unsigned, unsigned> TUPair;			
            TIndex index(seqs);
            String <unsigned> dislex;
			_makeDislexString(dislex, indexSA(index), indexText(index), DislexExternal<TShape, Skew7>());
			
			std::cout << "zrezz!" << std::endl;

			TLimits lim = stringSetLimits(seqs);
			//printUnsignedString(dislex);
			
			std::cout << length(dislex) << std::endl;
			
			String <unsigned> seq = prefix(dislex, 200);
			String <unsigned> read = suffix(dislex, 200);
			//printUnsignedString(seq);
			//printUnsignedString(read);
			std::cout << length(seq) << std::endl;
			std::cout << length(read) << std::endl;
			
			
			typedef Align<Dna5String, ArrayGaps> TAlignDna;      // align type
			TAlignDna align;
			resize(rows(align), 2);
			assignSource(row(align, 0), seqIn);
			assignSource(row(align, 1), readIn);
			int score = globalAlignment(align, Score<int, Simple>(1, 0, 0,-1));
			std::cout << "Score: " << score << std::endl;
			std::cout << align << std::endl;
				
			typedef Align<String<unsigned>, ArrayGaps> TAlignU;      // align type
			TAlignU align2;
			resize(rows(align2), 2);
			assignSource(row(align2, 0), seq);
			assignSource(row(align2, 1), read);
			int score2 = globalAlignment(align2, Score<int, Simple>(1, 0, 0,-1));
			std::cout << "Score: " << score2 << std::endl;
			std::cout << align2 << std::endl;
				
    return 0;
}
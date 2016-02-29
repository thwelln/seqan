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
#include <seqan/indexmy.h>


#include <time.h>

using namespace seqan;


int globalWrongMethods = 0;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class AppOptions
// --------------------------------------------------------------------------

struct AppOptions
{
    CharString infile;
    CharString mode;
    bool qsort;
    
    bool special1;
};


void printUnsignedString(String<unsigned> st)
{
	std::cout << "[";
	for (int i=0; i<length(st); i++)
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
        TSA &suffixArray,
        StringSet<TString, TSpec> const &stringSet,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
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
        typedef Pair<unsigned, unsigned> TUPair;
        typedef DislexReverseTransformMulti_<unsigned,
        TLimits, TUPair >                             TGetDislexReversePos;

        // instantiation and processing
        
        src_t        src(concat(stringSet));
        unsigner_t  unsigner(src);
        creator_t    creator(unsigner, stringSetLimits(stringSet));
        //String<unsigned> bla;
        suffixArray << creator;
        //std::cout << bla[0];        
        std::cout << creator.dislexString[43];
        printUnsignedString(creator.dislexString);
        printUnsignedString(creator.dislexString2);
        TUPair bub = TGetDislexReversePos(TShape::span, stringSetLimits(stringSet)) (43);
        std::cout << bub;
		
        suffixArray << creator;
        #ifdef SEQAN_TEST_INDEX
            //isSuffixArray(suffixArray, stringSet);
        #endif
    }


// --------------------------------------------------------------------------
// Function externalBenchmark()
// --------------------------------------------------------------------------

template <typename TIndex, typename TAlgo>
void externalBenchmark( TIndex & index,
                        TAlgo const & algo)
{
    double start = sysTime();
    _makeDislexString(indexSA(index), indexText(index), algo);
    std::cout << sysTime() - start << std::endl;
}



// --------------------------------------------------------------------------
// Function main()
// --------------------------------------------------------------------------

// Program entry point.

int main(int argc, char const ** argv)
{
 
    seqan::StringSet<seqan::CharString> ids;
    seqan::StringSet<seqan::Dna5String> seqs;

        
 
	CharString seqFileName = getAbsolutePath("~/seq.txt");
    CharString id;
    Dna5String seq;

    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seq, seqFileIn);
    		
			
            appendValue(seqs, seq);


	Dna5String seq2 = infixWithLength(seqStr, 200, 150);
    appendValue(seqs, seq2);


     typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<1,4,3,2,1,1> >, 1> > TShape;
        {
            typedef Index<StringSet<Dna5String> > TIndex;
            TIndex index(seqs);
            externalBenchmark(index, DislexExternal<TShape, Skew7>());
        }  
    return 0;
}

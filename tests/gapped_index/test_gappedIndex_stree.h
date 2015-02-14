// ==========================================================================
//                          test_gappedIndex_stree.h
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

#ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_STREE_H_
#define CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_STREE_H_

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


SEQAN_DEFINE_TEST(test_gappedIndex_top_down_traversal)
{
    typedef CyclicShape<FixedShape<0,GappedShape<HardwiredShape<3> >, 0> >      TShape; //1001
    typedef Index<Peptide, IndexSa<Gapped<ModCyclicShape<TShape> > > >  TIndex;        // Gapped<ModCyclicShape<TShape> >

    Peptide     s = "NHLDHCNMHLDHCQHLDHCCHLDHCKHLMLIHCVHLFKHCDHLEHC";
    TIndex index(s, TShape());            //, TShape()
    indexRequire(index, FibreSA());

    typedef Suffix<TIndex>::Type TModSuff;


    // Suffix array (ouput of:)
    //for (unsigned i=0; i<length(indexSA(index)); ++i)
    //    std::cout << "//  " << i << "\t" << indexSA(index)[i] << "\t" << suffix(index, indexSA(index)[i]) << std::endl;
    //
    //  0	45	C
    //  1	19	CDHHLIHHLHCL
    //  2	5	CHLCQDHHLCKMLCVFKDHH
    //  3	12	CLDCHHCLMHCLFCDE
    //  4	18	CLDKHLIVHKHHL
    //  5	39	CL
    //  6	32	CLFCDE
    //  7	24	CLMHCLFCDE
    //  8	16	DCHHCLMHCLFCDEH
    //  9	40	DE
    //  10	22	DKHLIVHKHHL
    //  11	3	DNMDHHLCCDHHLIHHLHCL
    //  12	10	DQHHCLDKHLIVHKHHL
    //  13	43	E
    //  14	36	FCDEH
    //  15	44	H                       *---
    //  16	41	HHC                         *---+---
    //  17	1	HHCHLCQDHHLCKMLCVFKDHHC
    //  18	8	HHCLDCHHCLMHCLFCDEH                 *---+---
    //  19	14	HHCLDKHLIVHKHHL
    //  20	20	HHCLMHCLFCDEH
    //  21	38	HHLC                            *---+---+---
    //  22	11	HHLCCDHHLIHHLHCL
    //  23	17	HHLCKMLCVFKDHH
    //  24	31	HHLHCL
    //  25	23	HHLIHHLHCL
    //  26	34	HKHHLC                      *---
    //  27	26	HLIVHKHHLC
    //  28	4	HMHHCLDCHHCLMHCLFCDEH
    //  29	30	IVHKHHLC                *---
    //  30	37	KDHHC
    //  31	25	KMLCVFKDHH
    //  32	42	LC
    //  33	15	LCCDHHLIHHLHCL
    //  34	21	LCKMLCVFKDHH
    //  35	2	LCNLDQHHCLDKHLIVHKHHL
    //  36	9	LCQDHHLCKMLCVFKDHH
    //  37	29	LCVFKDHH
    //  38	35	LHCL
    //  39	27	LIHHLHCL
    //  40	7	MDHHLCCDHHLIHHLHCLE
    //  41	28	MHCLFCDE
    //  42	0	NDHMHHCLDCHHCLMHCLFCDEH
    //  43	6	NLDQHHCLDKHLIVHKHHL
    //  44	13	QDHHLCKMLCVFKDHHC
    //  45	33	VFKDHHC


    Iterator<TIndex, TopDown<> >::Type treeIter(index);

    goDown(treeIter, 'H');
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(15,29));

    goDown(treeIter, 'H');
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(16,26));

    SEQAN_ASSERT_EQ(representative(treeIter), "HH");
    SEQAN_ASSERT_EQ(repLength(treeIter), 2u);

    goDown(treeIter, 'C');
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(16,21));

    goRight(treeIter);
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(21,26));

    goRoot(treeIter);
    goDown(treeIter, "HHLC");
    SEQAN_ASSERT_EQ(range(treeIter), Pair<unsigned>(21,24));
    SEQAN_ASSERT_EQ(representative(treeIter), "HHLC");
    SEQAN_ASSERT_EQ(repLength(treeIter), 4u);


    Iterator<TIndex, TopDown< ParentLinks<Preorder> > >::Type myIterator(index);
    for (unsigned i=0; i<13; ++i, ++myIterator)
    {}

    SEQAN_ASSERT_EQ(representative(myIterator), "CDHHLIHHLHCLE");
    SEQAN_ASSERT(_isLeaf(myIterator, HideEmptyEdges()));
    SEQAN_ASSERT_NOT(goDown(myIterator));

    ++myIterator;
    SEQAN_ASSERT_EQ(representative(myIterator), "CH");

    Iterator<TIndex, TopDown< ParentLinks<Preorder> > >::Type iter2 = myIterator;

    ++myIterator;
    for (; !atEnd(myIterator); ++myIterator, ++iter2)
    {}
    SEQAN_ASSERT_EQ(representative(iter2), "VFKDHHC");

}


#endif  // #ifndef CORE_TESTS_GAPPEDINDEX_TEST_GAPPEDINDEX_STREE_H_

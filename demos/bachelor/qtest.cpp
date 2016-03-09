#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>

#include <time.h>

using namespace seqan;

typedef GappedShape<HardwiredShape<1> > TInsideShape; // SHAPE: 11 SPAN: 2 

Dna5String seq = "ACGACTACGCTACGATCGTATATTACGCAATGCA"
Dna5String read = "ACTATACG"

Index<Dna5String, IndexQGram<TInsideShape> > qgindex(seq);

Dna5String qgram1 = infixWithLength(read, 0, 2);
Dna5String qgram2 = infixWithLength(read, 0, 4);

hash(indexShape(qgindex),begin(qgram1));
std::cout << countOccurences(qgindex, indexShape(qgindex)) << std::endl;

hash(indexShape(qgindex),begin(qgram2));
std::cout << countOccurences(qgindex, indexShape(qgindex)) << std::endl;

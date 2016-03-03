#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/index.h>
#include <seqan/modifier.h>
#include <random>	

using namespace seqan;

int main(int argc, char *argv[])
{
	typedef GappedShape<HardwiredShape<1,1> > TInsideShape;
	Dna5String seq = "ACGATCATAGCCATACG";
	Dna5String read = "CAT";
	typedef Index<Dna5String, IndexQGram<TInsideShape > > TIndex;
    TIndex qgindex(seq);
    
    hash(indexShape(qgindex), begin(read));    
    //hash(indexShape(qgindex), "CAT");
    for (unsigned i = 0; i < length(getOccurrences(qgindex, indexShape(qgindex))); ++i)
        std::cout << getOccurrences(qgindex, indexShape(qgindex))[i] << std::endl;

    return 0;
}

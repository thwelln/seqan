#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/arg_parse.h>
#include <seqan/pipe.h>
#include <seqan/index.h>


#include <time.h>

using namespace seqan;
int main()
{
String <unsigned> seq;
for (unsigned i=0; i<400; ++i)
{
	appendValue(seq, i);
}
	typedef Index<String <unsigned>, FMIndex<> > TFMIndex;
	TFMIndex fmindex(seq);
	std::cout << "BLA!" << std::endl;
	
	Iterator<TFMIndex, TopDown<> >::Type fmit(fmindex);
	std::cout << "BLUB!" << std::endl;

}

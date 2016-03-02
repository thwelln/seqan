#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <stdlib.h>
#include <random>	

using namespace seqan;

int changeBase(DnaString *seqin, unsigned position, unsigned changevalue)
{
	(*seqin)[position] = Dna((((*seqin)[position]).value+changevalue)%4);
	return 0;
}

int main(int argc, char *argv[])
{
    CharString seqFileName = getAbsolutePath("/../Sequences/sequence.fasta");
    char outpath [256];  
    CharString id;
    Dna5String seq;

    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seq, seqFileIn);
    
    DnaString seq2 = seq;

			
		unsigned startpos = atoi(argv[1]);
		unsigned readlength = atoi(argv[2]);
		double error_rate = atof(argv[3]);
		unsigned seed = atoi(argv[4]);		
			
	    DnaString read = infixWithLength(seq2, startpos, readlength);
	    
	    std::cout << read << std::endl;
	    
	    std::mt19937 gen(seed);
	    
	    std::uniform_int_distribution<unsigned> posdist(0,length(read)-1);
	    std::uniform_int_distribution<unsigned> basedist(1,3);
	    
	       
	    unsigned erroramount = (unsigned) (readlength*error_rate);
	    std::cout << erroramount << std::endl;
	    for (unsigned i=1;i<=erroramount;i++)
	    {
			unsigned pos = posdist(gen);
			unsigned change = basedist(gen);
			std::cout << i <<": " << pos << ", " << read[pos] << "->";
			changeBase(&read, pos, change);
			std::cout << read[pos] << std::endl;
	    }
	    
	    std::cout << read << std::endl;
	    sprintf(outpath, "/../Sequences/Reads/read_%d_%d_%.2f.fasta",startpos,readlength,error_rate);
	    CharString readFileName = getAbsolutePath(outpath);  
	    SeqFileOut seqFileOut(toCString(readFileName));
	    writeRecord(seqFileOut, id ,read);

    
    return 0;
}

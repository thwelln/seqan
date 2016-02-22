#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <random>	

using namespace seqan;

int changeBase(DnaString *seqin, unsigned position, unsigned changevalue)
{
	(*seqin)[position] = Dna((((*seqin)[position]).value+changevalue)%4);
	return 0;
}

int main()
{
	unsigned startpos_array [] = {270000, 500000};
	unsigned readlength_array [] = {50,100,150,200};
	unsigned error_rate_array [] = {5,10,20,50}; //one error per amount
    CharString seqFileName = getAbsolutePath("/../Sequences/sequence.fasta");
    char outpath [256];  
    CharString id;
    Dna5String seq;

    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seq, seqFileIn);
    
    DnaString seq2 = seq;
    for (int startmode = 0; startmode<2; startmode++)
    {
	for (int lengthmode = 0; lengthmode<4; lengthmode++)
	{
	for (int errormode = 0; errormode<4; errormode++)
	{
			
		unsigned startpos = startpos_array[startmode];
		unsigned readlength = readlength_array[lengthmode];
		unsigned error_rate = error_rate_array[errormode];	
			
	    DnaString read = infixWithLength(seq2, startpos, readlength);
	    
	    std::cout << read << std::endl;
	    
	    std::mt19937 gen(41);
	    
	    std::uniform_int_distribution<unsigned> posdist(0,length(read)-1);
	    std::uniform_int_distribution<unsigned> basedist(1,3);    
	    
	    for (unsigned i=0;i<=readlength/error_rate;i++)
	    {
			unsigned pos = posdist(gen);
			unsigned change = basedist(gen);
			std::cout << i <<": " << pos << ", " << read[pos] << "->";
			changeBase(&read, pos, change);
			std::cout << read[pos] << std::endl;
	    }
	    
	    std::cout << read << std::endl;
	    sprintf(outpath, "/../Sequences/Reads/read_%d_%d_%d.fasta",startpos,readlength,error_rate);
	    CharString readFileName = getAbsolutePath(outpath);  
	    SeqFileOut seqFileOut(toCString(readFileName));
	    writeRecord(seqFileOut, id ,read);
	}
	}
	}    
    
    return 0;
}

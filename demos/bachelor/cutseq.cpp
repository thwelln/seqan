#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <stdlib.h>
#include <random>	

using namespace seqan;

int main(int argc, char *argv[])
{
    CharString seqFileName = getAbsolutePath("/../Sequences/sequence.fasta");
    char outpath [256];  
    CharString id;
    Dna5String seq;
    
    unsigned startpos = atoi(argv[1]);
    unsigned readlength = atoi(argv[2]);   

    SeqFileIn seqFileIn(toCString(seqFileName));
    readRecord(id, seq, seqFileIn);
    
	    DnaString5 read = infixWithLength(seq, startpos, readlength);

	    sprintf(outpath, "/../Sequences/incoming.fasta");
	    CharString readFileName = getAbsolutePath(outpath);  
	    SeqFileOut seqFileOut(toCString(readFileName));
	    writeRecord(seqFileOut, id ,read);

    
    return 0;
}

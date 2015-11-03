#ifndef FASTA_READER_H
#define FASTA_READER_H

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include "unistd.h"
#include "fcntl.h"

#ifndef KSEQ
// KSEQ_INIT(gzFile, gzread)
KSEQ_INIT(int, read)
#endif

class FastaReader {
    FILE * fp;
    // int fp; // file handler
    kseq_t *seq;
    int l;

public:

    FastaReader(const char * fname) {
        fp = fopen(fname, "r"); // TODO: check params
        if (fp == NULL) {
            printf("Could not open file %s\n", fname);
            return;
        }
        seq = kseq_init( fileno(fp) );
    }

    ~FastaReader() {
        kseq_destroy(seq);
        fclose( fp );
    }

    kseq_t * nextSequence() {
        l = kseq_read(seq);
        if (l < 0) return NULL;
        else return seq;
    }
};

#endif /* FASTA_READER_H */
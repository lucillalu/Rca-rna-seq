#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "abpoa.h"

int seq_msa(int n_seqs, uint8_t** bseqs, int* seq_lens) {
    int i, j;
    // for (i = 0; i < n_seqs; ++i)
    // {
    //     printf("seq_lens[%d] = %d\n", i, seq_lens[i]);
    // }
    
    // for (i = 0; i < n_seqs; ++i)
    // {
    //     for (j = 0; j < seq_lens[i]; ++j)
    //     {
    //         printf("%c", "ACGTN-"[bseqs[i][j]]);
    //     }
    //     printf("\n");
    // }
    // variables to store result
    uint8_t **cons_seq; int *cons_l, cons_n=0;
    uint8_t **msa_seq; int msa_l=0;

   // initialize variables
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
     
    // output options
    abpt->out_msa = 1; // generate Row-Column multiple sequence alignment(RC-MSA), set 0 to disable
    // abpt->out_cons = 0; // generate consensus sequence, set 0 to disable

    abpoa_post_set_para(abpt);

    // output to stdout
    // fprintf(stdout, "=== output to stdout ===\n");

    // perform abpoa-msa
    // abpoa_msa(ab, abpt, n_seqs, seq_lens, bseqs, stdout, NULL, NULL, NULL, NULL, NULL);

    // abpoa_reset_graph(ab, abpt, seq_lens[0]); // reset graph before re-use

    // perform abpoa-msa
    abpoa_msa(ab, abpt, n_seqs, seq_lens, bseqs, NULL, &cons_seq, &cons_l, &cons_n, &msa_seq, &msa_l);

    // fprintf(stdout, "=== output to variables ===\n");
    // printf("cons_n = %d\n", cons_n);
    for (i = 0; i < cons_n; ++i) {
        fprintf(stdout, ">Consensus_sequence\n");
        for (j = 0; j < cons_l[i]; ++j)
            fprintf(stdout, "%c", "ACGTN"[cons_seq[i][j]]);
        fprintf(stdout, "\n");
    }
    fprintf(stdout, ">Multiple_sequence_alignment\n");
    for (i = 0; i < n_seqs; ++i) {
        for (j = 0; j < msa_l; ++j) {
            fprintf(stdout, "%c", "ACGTN-"[msa_seq[i][j]]);
        }
        fprintf(stdout, "\n");
    }

    if (cons_n) {
        for (i = 0; i < cons_n; ++i) free(cons_seq[i]); 
        free(cons_seq); free(cons_l);
    }
    if (msa_l) {
        for (i = 0; i < n_seqs; ++i) free(msa_seq[i]); free(msa_seq);
    }

    /* generate DOT partial order graph plot */
    // abpt->out_pog = strdup("example.png"); // dump parital order graph to file
    // if (abpt->out_pog != NULL) abpoa_dump_pog(ab, abpt);

    for (i = 0; i < n_seqs; ++i) free(bseqs[i]); free(bseqs); free(seq_lens);
    abpoa_free(ab, abpt); abpoa_free_para(abpt); 
    return 0;
}

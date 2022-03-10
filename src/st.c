#include <stdio.h>

#include "fasta.h"
#include "fastq.h"
#include "sam.h"
#include <cstr.h>

// We hardwire that the alphabet is acgt + $
// and have a global variable for the alphabet.
// It is initialised in main().
static cstr_alphabet ALPHA;

// Linked list of preprocessed trees
struct suf_tree
{
    const char *chr_name;  // chromosome name; data stored in fasta record
    cstr_sslice *x_buf;    // buffer for the remapped sequence. We own it.
    cstr_suffix_tree *st;  // the suffix tree. We also own it.
    struct suf_tree *next; // next suffix tree in the linked list.
};
static struct suf_tree *process_fasta_rec(struct fasta_record *fa_rec,
                                          struct suf_tree *next)
{
    struct suf_tree *suf_tree = cstr_malloc(sizeof *suf_tree);
    suf_tree->chr_name = fa_rec->name;
    suf_tree->x_buf = cstr_alloc_sslice(fa_rec->seq.len);
    bool map_ok = cstr_alphabet_map(*suf_tree->x_buf, fa_rec->seq, &ALPHA);
    if (!map_ok)
    {
        abort(); // A little drastic, but with the current choice it is all that I can do.
    }
    suf_tree->st = cstr_mccreight_suffix_tree(&ALPHA, CSTR_SLICE_CONST_CAST(*suf_tree->x_buf));
    suf_tree->next = next;
    return suf_tree;
}
static void free_suf_trees(struct suf_tree *suf_tree)
{
    struct suf_tree *next;
    for (; suf_tree; suf_tree = next)
    {
        next = suf_tree->next;
        // don't free chr_name, we don't own it.
        free(suf_tree->x_buf);
        cstr_free_suffix_tree(suf_tree->st);
        free(suf_tree);
    }
}

int main(int argc, char const *argv[])
{
    if (argc != 3)
    {
        fprintf(stderr, "%s genome reads", argv[0]);
        return 1;
    }
    const char *genome_fname = argv[1];
    const char *reads_fname = argv[2];

    // Initialise global alphabet.
    // NB! This assumes that we know the alphabet already and that it is "acgt" + $.
    // If we don't know it, we need to obtain it from the chromosomes.
    cstr_const_sslice letters = CSTR_SLICE_STRING0((const char *)"acgt");
    cstr_init_alphabet(&ALPHA, letters);

    struct fasta_records *genome = load_fasta_records(genome_fname);
    struct suf_tree *s_trees = 0;
    for (struct fasta_record *fa_rec = fasta_records(genome); fa_rec; fa_rec = fa_rec->next)
    {
        s_trees = process_fasta_rec(fa_rec, s_trees);
    }

    struct fastq_iter fqiter;
    struct fastq_record fqrec;
    static const size_t CIGAR_BUF_SIZE = 2048;
    char cigarbuf[CIGAR_BUF_SIZE];

    FILE *fq = fopen(reads_fname, "r");
    if (!fq)
    {
        abort(); // I can always implement better checking another day...
    }
    init_fastq_iter(&fqiter, fq);

    while (next_fastq_record(&fqiter, &fqrec))
    {
        snprintf(cigarbuf, CIGAR_BUF_SIZE, "%lldM", fqrec.seq.len);

        for (struct suf_tree *suf_tree = s_trees; suf_tree; suf_tree = suf_tree->next)
        {
            cstr_exact_matcher *m = cstr_st_exact_search_map(suf_tree->st, fqrec.seq);
            for (long long i = cstr_exact_next_match(m); i != -1; i = cstr_exact_next_match(m))
            {
                print_sam_line(stdout, (const char *)fqrec.name.buf,
                               suf_tree->chr_name, i, cigarbuf,
                               (const char *)fqrec.seq.buf);
            }
            cstr_free_exact_matcher(m);
        }
    }

    fclose(fq);
    free_suf_trees(s_trees);
    free_fasta_records(genome);

    return 0;
}

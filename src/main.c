/* insure strdup nad strndup are loaded correctly */
#ifdef __STDC_ALLOC_LIB__
#define __STDC_WANT_LIB_EXT2__ 1
#else
#define _POSIX_C_SOURCE 200809L
#endif

#include <stdlib.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <zlib.h>

#include "gd.h"

#include "args.h"
#include "palette.h"

#include "minimap2/minimap.h"
#include "minimap2/kseq.h"
#include "minimap2/kvec.h"
#include "minimap2/khash.h"

#define QRY 0
#define REF 1

KSEQ_INIT(gzFile, gzread);

typedef struct {
    struct { int id; long start, end;} qry, ref;
    int reversed, score, mapq;
    long matches, length;
} alignment_t;

typedef struct {
    int index, filter, order;
    long length, offset;
    struct { long total, *indv; } coverage;
    struct { int id; long start; } closest;
    char * name;
} sequence_t;

long* get_offsets(kvec_t(sequence_t)* seqs, size_t * map ){
    size_t i, size;
    size_t * index_map = map;

    size = kv_size(*seqs);

    /* If no mapping is given, setup 1-1 mapping */
    if(index_map == NULL){
        index_map = malloc(size * sizeof(size_t));
        for(i = 0; i < size; i++)
            index_map[i] = i;
    }

    /* order lengths, array is oversized so that index =size= will be the total
     * length of all included sequences. Also makes the code for dealing with
     * filtered sequences less cumbersome*/
    long * ordered_offsets = calloc(size + 1, sizeof(long));
    for( i = 0; i < size; i++){
        sequence_t * seq = &kv_A(*seqs, i);
        ordered_offsets[index_map[seq->index]] = seq->length;
    }

    /* convert lengths to offsets */
    for(i = 1; i <= size; i ++)
        ordered_offsets[i] += ordered_offsets[i-1];

    /* shift offsets */
    for(i = size-1; i >= 1; i --)
        ordered_offsets[i] = ordered_offsets[i-1];
    ordered_offsets[0] = 0;

    /* create original index ordered offset array  */
    long * offsets = malloc((size+1) * sizeof(long));
    for(i = 0; i < size; i ++)
        offsets[i] = ordered_offsets[index_map[i]];
    offsets[size] = ordered_offsets[size];

    free(ordered_offsets);
    if(map == NULL)
        free(index_map);

    return offsets;
}

int main(int argc, char *argv[]) {
    arguments_t args = parse_options(argc, argv);

    mm_idxopt_t iopt;
    mm_mapopt_t mopt;

    mm_verbose = 2; // disable message output to stderr
    mm_set_opt(NULL, &iopt, &mopt); // init options
    mm_set_opt("asm5", &iopt, &mopt); // set options to asm5 presets

    // open index reader
    mm_idx_reader_t *r = mm_idx_reader_open(args.ref, &iopt, NULL);

    gzFile f = gzopen(args.qry, "r");
    kseq_t *ks = kseq_init(f);

    kvec_t(alignment_t) alns;
    kv_init(alns);

    kvec_t(sequence_t) seqs[2];
    kv_init(seqs[REF]);
    kv_init(seqs[QRY]);
    while (kseq_read(ks) >= 0) {
        sequence_t tmp = { .index = kv_size(seqs[QRY]),
        .length = ks->seq.l,
        .name = strndup(ks->name.s, ks->name.l)};

        kv_push(sequence_t, NULL, seqs[QRY], tmp);
    }


    mm_idx_t *mi;
    while ((mi = mm_idx_reader_read(r, args.threads)) != 0) { // traverse each part of the index
        mm_mapopt_update(&mopt, mi); // this sets the maximum minimizer occurrence; TODO: set a better default in mm_mapopt_init()!

        /* Save ref names in index list */
        unsigned int i;
        for(i = 0; i < mi->n_seq; i++) {
            sequence_t tmp = { .index = i,
                               .length = mi->seq[i].len,
                               .name = strdup(mi->seq[i].name)};

            kv_push(sequence_t, NULL, seqs[REF], tmp);
        }

#pragma omp parallel           \
    num_threads(args.threads)                           \
    default(none)                                       \
    shared(mi, mopt, ks, f, args, alns)
        {

#pragma omp single nowait
            {
                gzrewind(f);
                kseq_rewind(ks);


                int qry_index = 0;
                /* Save qry names in index list */
                while (kseq_read(ks) >= 0) { // each kseq_read() call reads one query sequence

                    char * sequence = calloc(ks->seq.m, sizeof(char));
                    strncpy(sequence, ks->seq.s, ks->seq.m);
                    size_t length = ks->seq.l;


#pragma omp task                                \
    default(none)                               \
    firstprivate(sequence, length, qry_index)                 \
    shared(mi, mopt, args, alns)
                    {
            mm_tbuf_t *tbuf = mm_tbuf_init(); // thread buffer; for multi-threading, allocate one tbuf for each thread
                        mm_reg1_t *reg;
                        int j, n_reg;
                        reg = mm_map(mi, length, sequence, &n_reg, tbuf, &mopt, 0); // get all hits for the query
                        for (j = 0; j < n_reg; ++j) { // traverse hits and print them out
                            mm_reg1_t *r = &reg[j];

                            if(r->blen > args.minimum_align_length){
                                alignment_t aln = { .qry.id    = qry_index,
                                                    .qry.start = r->qs,
                                                    .qry.end   = r->qe,

                                                    .ref.id    = r->rid,
                                                    .ref.start = r->rs,
                                                    .ref.end   = r->re,

                                                    .reversed  = r->rev,
                                                    .matches   = r->mlen,
                                                    .length    = r->blen,
                                                    .score     = r->score,
                                                    .mapq      = r->mapq
                                                   };

#pragma omp critical
                                {
                                    kv_push(alignment_t, NULL, alns, aln);
                                }
                            }
                        }
                        free(reg);
                        mm_tbuf_destroy(tbuf);
                        free(sequence);
                    }
                    qry_index++;
                }
            }
        }
        mm_idx_destroy(mi);
    }
    mm_idx_reader_close(r); // close the index reader

    printf("Sizes %d %d\n", kv_size(seqs[QRY]), kv_size(seqs[REF]));

    /* alloc space for indv coverage */
    size_t i;
    for(i = 0; i < kv_size(seqs[QRY]); i++){
        sequence_t* seq = &kv_A(seqs[QRY], i);
        seq->coverage.total = 0;
        seq->coverage.indv = calloc(kv_size(seqs[REF]), sizeof(long));
    }

    /* Accumulate total coverage for each query */
    for(i = 0; i < kv_size(alns); i++){
        alignment_t* aln = &kv_A(alns, i);
        sequence_t* seq = &kv_A(seqs[QRY], aln->qry.id);
        seq->coverage.total += aln->length;
        seq->coverage.indv[aln->ref.id] += aln->length;
    }


    /* Find best match from ref for each qry sequence */
    size_t j;
    for(i = 0; i < kv_size(seqs[QRY]); i++){
        sequence_t* seq = &kv_A(seqs[QRY], i);

        long max = 0;
        size_t idx = kv_size(seqs[REF]);

        /* skip low coverage query contigs */
        if(seq->coverage.total > args.minimum_query_length){

            /* find ref contig index with largest coverage */
            for(j = 0; j <= kv_size(seqs[REF]); j++){
                if(max < seq->coverage.indv[j]){
                    max = seq->coverage.indv[j];
                    idx = j;
                }
            }
        }

        seq->closest.id = idx;
        seq->closest.start = LONG_MAX;
    }

    /* find earliest alignment for the closest reference to query seq */
    long score = 0, passed = 0;
    for(i = 0; i < kv_size(alns); i++){
        alignment_t* aln = &kv_A(alns, i);
        sequence_t* seq = &kv_A(seqs[QRY], aln->qry.id);

        if(aln->ref.id != seq->closest.id) continue;
        if(aln->ref.start > seq->closest.start) continue;

        seq->closest.start = aln->ref.start;
        score += aln->score;
        passed++;
    }
    score = 2 * score / passed;


    /* Create query index map, removing (setting to size_qry) index that don't
     * have enough coverage */
    struct seq_list_s {struct seq_list_s *next; sequence_t * seq;} *seq_list;
    seq_list = NULL;

    for( i = 0; i < kv_size(seqs[QRY]); i++ ){
        sequence_t * seq = &kv_A(seqs[QRY], i);
        if(seq->closest.id == kv_size(seqs[REF])) continue;

        struct seq_list_s * tmp = malloc(sizeof(struct seq_list_s));
        tmp->next = NULL;
        tmp->seq = seq;

        struct seq_list_s * cur = seq_list;
        struct seq_list_s * prev = NULL;

        /* Loop through linked list while
         * 1) there is an item
         * 2) the list item closest ref id is less then the seq closest ref id
         * 3) the list item starts before seq on the same ref */
        while(cur != NULL &&
              ( cur->seq->closest.id < tmp->seq->closest.id ||
                ( cur->seq->closest.id == tmp->seq->closest.id &&
                  cur->seq->closest.start < tmp->seq->closest.start))){
            prev = cur;
            cur = cur->next;
        }

        /* Insert seq into correct position */
        tmp->next = cur;
        if(prev == NULL){
            seq_list = tmp;
        }else{
            prev->next = tmp;
        }

        /* print out list (debug) */
        /* printf("Added %s to list [%x]\n", seq->name, seq_list); */
        /* cur = seq_list; */
        /* while(cur != NULL){ */
        /*     printf( "\t%x = { .next = %x; .seq = %x = {.name = %s; .index = %d}}\n", */
        /*             cur, cur->next, cur->seq, cur->seq->name, cur->seq->index ); */
        /*     cur = cur->next; */
        /* } */
    }

    /* Create index_map */
    size_t * index_map = malloc(kv_size(seqs[QRY]) * sizeof(size_t));
    for( i = 0; i < kv_size(seqs[QRY]); i++)
        index_map[i] = kv_size(seqs[QRY]);

    i = 0;
    while(seq_list){
        index_map[seq_list->seq->index] = i++;

        struct seq_list_s * tmp = seq_list->next;
        free(seq_list);
        seq_list = tmp;
    }

    long * offset[2] = { get_offsets(&(seqs[QRY]), index_map),
                         get_offsets(&(seqs[REF]), NULL     ) };

    int bin_size[2] = { (offset[QRY][kv_size(seqs[QRY])]/args.size) + 1,
                        (offset[REF][kv_size(seqs[REF])]/args.size) + 1 };



    /* Create image and palette */
    gdImagePtr img = gdImageCreate(args.size, args.size);
    int * colors = load_palette(img, args.pal);

    /* draw lines */
    for( i = 0; i < kv_size(seqs[QRY]); i++ ){
        if(offset[QRY][i] == offset[QRY][kv_size(seqs[QRY])]) continue;
        long coord[4] = { 0,         offset[QRY][i]/bin_size[QRY],
                          args.size, offset[QRY][i]/bin_size[QRY]};
        gdImageLine( img, coord[0], coord[1], coord[2], coord[3], colors[25] );
    }
    for( i = 0; i < kv_size(seqs[REF]); i++ ){
        if(offset[REF][i] == offset[REF][kv_size(seqs[REF])]) continue;
        long coord[4] = { offset[REF][i]/bin_size[REF], 0,
                          offset[REF][i]/bin_size[REF], args.size};
        gdImageLine( img, coord[0], coord[1], coord[2], coord[3], colors[25] );
    }

    /* draw alignments */
    for(i = 0; i < kv_size(alns); i++){
        alignment_t* aln = &kv_A(alns, i);

        struct {struct {long x,y;} str,end; } pnt;

        pnt.str.x = (aln->ref.start + offset[REF][aln->ref.id]) / bin_size[REF];
        pnt.end.x = (aln->ref.end   + offset[REF][aln->ref.id]) / bin_size[REF];
        pnt.str.y = (aln->qry.start + offset[QRY][aln->qry.id]) / bin_size[QRY];
        pnt.end.y = (aln->qry.end   + offset[QRY][aln->qry.id]) / bin_size[QRY];

        if(aln->reversed){
            long tmp = pnt.str.y;
            pnt.str.y = pnt.end.y;
            pnt.end.y = tmp;
        }

        unsigned long color_idx = 255 * aln->score / aln->matches;
        if( color_idx > 255 ) color_idx = 255;
        if(!args.color_similarity) color_idx = 200;

        printf("%d\t%d\t%d\t%f\t%d\n", aln->score, aln->length, aln->matches, 1.0*aln->score/aln->matches, color_idx);
        gdImageLine( img,
                     pnt.str.x, pnt.str.y,
                     pnt.end.x, pnt.end.y,
                     colors[color_idx] );
    }


    /* Write image */
    FILE* out = fopen(args.out, "wb");
    switch(args.type){
        case png:  gdImagePng(img, out);     break;
        case jpeg: gdImageJpeg(img, out, 95); break;
        case tiff: gdImageTiff(img, out);    break;
        case bmp:  gdImageBmp(img, out, 0);  break;
    }
    fclose(out);

    printf("Number of alignments ... %d\n"
           "Number of reference seqs ... %d\n"
           "Number of query seqs ... %d\n",
           alns.n, kv_size(seqs[REF]), kv_size(seqs[QRY]));
    gdImageDestroy(img);
    return 0;
}

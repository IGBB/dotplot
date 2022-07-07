#ifndef ARGS_H_
#define ARGS_H_

#include <stdio.h>

extern const char* const program_version;

typedef enum {
    png, bmp, jpeg, tiff
} outtype_t;

typedef struct {
    char * region, *ref, *qry, *out, *font;
    int size, threads, minimum_align_length, minimum_query_length, color_similarity;
    outtype_t type;
    const int *pal;
} arguments_t;

arguments_t parse_options(int argc, char **argv);


#endif // ARGS_H_

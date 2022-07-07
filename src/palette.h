#ifndef PALETTE_H_
#define PALETTE_H_

#include "gd.h"

extern const int magma[256];
extern const int inferno[256];
extern const int mako[256];
extern const int rocket[256];
extern const int grey[256];

int * load_palette(gdImagePtr, const int[] );

#endif // PALETTE_H_

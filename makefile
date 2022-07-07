src = $(wildcard src/*.c)
obj = $(src:.c=.o)

CFLAGS  += -Wall -std=c99
LDFLAGS += -Lsrc/minimap2 -lminimap2 -lgd -lpng -lfreetype -lm -lz -pthread -std=c99

# Optimizations
CFLAGS  += -O2 -ggdb -fgnu89-inline -std=c99 -march=native -mtune=native -fopenmp
LDFLAGS += -O2 -ggdb -std=c99 -march=native -mtune=native -fopenmp

dotplot: $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) polycat

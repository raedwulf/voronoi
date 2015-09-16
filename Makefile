CC = cc
CFLAGS = -std=gnu11 -g -I.

LIBS=-lm
ARCH=$(uname -m)

all: voronoi-example

%.o: %.c
	$(CC) -c -o $@ $< $(CFLAGS)

voronoi-example: voronoi-example.o voronoi.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f *.o voronoi-example

/*
 * The author of this software is Shane O'Sullivan.  
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 *
 * Slightly modified to fit C API by Tai Chi Minh Ralph Eastwood.
 */

#include <stdio.h>
#include <search.h>
#include <malloc.h>
#include "voronoi.h"

int main (int argc,char **argv) 
{	

	struct point values[4] = {{-22,-9}, {-17,31}, {4,13}, {22,-5}};
	int count = 4;

	struct graph_edge *e = voronoi(values, count, -100, 100, -100, 100, 3);

	float x1,y1,x2,y2;

	printf("\n-------------------------------\n");
	for (struct graph_edge *e1 = e; e1; e1 = e1->next)
		printf("GOT Line (%f,%f)->(%f,%f)\n",e1->x1,e1->y1,e1->x2,e1->y2);
	
	voronoi_free(e);

	return 0;
}

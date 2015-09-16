/*
 * Copyright (c) 2015 by Tai Chi Minh Ralph Eastwood.
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose without fee is hereby granted, provided that this entire notice
 * is included in all copies of any software which is or includes a copy
 * or modification of this software and in all copies of the supporting
 * documentation for such software.
 * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED
 * WARRANTY.  IN PARTICULAR, NEITHER THE AUTHORS NOR AT&T MAKE ANY
 * REPRESENTATION OR WARRANTY OF ANY KIND CONCERNING THE MERCHANTABILITY
 * OF THIS SOFTWARE OR ITS FITNESS FOR ANY PARTICULAR PURPOSE.
 */

#ifndef __VORONOI_H__
#define __VORONOI_H__

struct point {
	float x,y;
};

struct graph_edge {
	float x1,y1,x2,y2;
	struct graph_edge* next;
};

struct graph_edge *voronoi(struct point *xy, int nsites, float min_x, float max_x, float min_y, float max_y, float min_dist);
void voronoi_free(struct graph_edge *e);

#endif

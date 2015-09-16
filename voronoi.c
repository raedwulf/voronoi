/**
 * The author of this software is Steven Fortune.
 * Copyright (c) 1994 by AT&T Bell Laboratories.
 *
 * Later modifications: Shane O'Sullivan, Tai Chi Minh Ralph Eastwood.
 * Full License Headers reproduced in LICENSE.md.
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

#include "voronoi.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define DELETED -2

#define le 0
#define re 1

struct free_node {
	struct free_node *nextfree;
};

struct free_node_al {
	struct free_node* memory;
	struct free_node_al* next;
};

struct free_list {
	struct free_node	*head;
	int		nodesize;
};

/* structure used both for sites and for vertices */
struct site {
	struct point	coord;
	int		sitenbr;
	int		refcnt;
};

struct edge {
	float a,b,c;
	struct site 	*ep[2];
	struct site	*reg[2];
	int		edgenbr;

};

struct half_edge {
	struct half_edge	*el_left, *el_right;
	struct edge	*el_edge;
	int		el_refcnt;
	char	el_pm;
	struct site	*vertex;
	float	ystar;
	struct half_edge *pq_next;
};

struct beg {
	struct free_node_al **b;
	int sqrt_nsites;
};

struct el {
	struct beg b;
	struct free_list hfl;
	struct half_edge **hash;
	struct half_edge *leftend, *rightend;
	int hashsize;
	int xmin, deltax;

};

struct geom {
	struct beg b;
	struct free_list efl, sfl;
	struct graph_edge *edges;
	struct site *bottomsite;
	int nvertices, nedges;
	float borderMinX, borderMaxX, borderMinY, borderMaxY, min_dist;
};

struct pq {
	int pq_hashsize;
	struct half_edge *pq_hash;
	int pq_count;
	int pq_min;
	int ymin, deltay;
};

static void makefree(struct free_node *curr, struct free_list *fl)
{
	curr->nextfree = fl->head;
	fl->head = curr;
}

static void *getfree(struct beg *beg, struct free_list *fl)
{
	struct free_node *t;
	struct free_node_al **b = beg->b, *next;

	if (!fl->head) {
		if (!(t = malloc(beg->sqrt_nsites * fl->nodesize)))
			return NULL;
		next = *b;
		if (!(*b = malloc(sizeof(struct free_node_al)))) {
			free(t);
			return NULL;
		}
		(*b)->memory = t;
		(*b)->next = next;

		for (int i = 0; i < beg->sqrt_nsites; i++)
			makefree((struct free_node *)((char *)t + i * fl->nodesize), fl);
	}
	t = fl->head;
	fl->head = fl->head->nextfree;
	return t;
}

static struct half_edge* he_create(struct el *el, struct edge *e, int pm)
{
	struct half_edge *answer = getfree(&el->b, &el->hfl);
	answer->el_edge = e;
	answer->el_pm = pm;
	answer->pq_next = NULL;
	answer->vertex = NULL;
	answer->el_refcnt = 0;
	return answer;
}

int el_init(struct el *el, struct free_node_al **b, int sqrt_nsites, int xmin, int deltax)
{
	el->b.b = b;
	el->b.sqrt_nsites = sqrt_nsites;
	el->hfl.head = NULL;
	el->hfl.nodesize = sizeof **el->hash;
	el->hashsize = 2 * sqrt_nsites;
	if (!(el->hash = malloc(sizeof *el->hash * el->hashsize)))
		return 0;

	for (int i = 0; i < el->hashsize; i += 1)
		el->hash[i] = NULL;
	el->leftend = he_create(el, NULL, 0);
	el->rightend = he_create(el, NULL, 0);
	el->leftend->el_left = NULL;
	el->leftend->el_right = el->rightend;
	el->rightend->el_left = el->leftend;
	el->rightend->el_right = NULL;
	el->hash[0] = el->leftend;
	el->hash[el->hashsize - 1] = el->rightend;
	el->xmin = xmin;
	el->deltax = deltax;

	return 1;
}

void el_free(struct el *el)
{
	if (el->hash)
		free(el->hash);
}

static void el_insert(struct half_edge *lb, struct half_edge *newHe)
{
	newHe->el_left = lb;
	newHe->el_right = lb->el_right;
	(lb->el_right)->el_left = newHe;
	lb->el_right = newHe;
}

/* Get entry from hash table, pruning any deleted nodes */
struct half_edge *el_gethash(struct el *el, int b)
{
	if (b < 0 || b >= el->hashsize)
		return NULL;
	struct half_edge *he = el->hash[b];
	if (!he || he->el_edge != (struct edge *)DELETED)
		return (he);

	/* Hash table points to deleted half edge.  Patch as necessary. */
	el->hash[b] = NULL;
	if (!(he->el_refcnt -= 1))
		makefree((struct free_node *)he, &el->hfl);
	return NULL;
}

/* returns 1 if p is to right of halfedge e */
int right_of(struct half_edge *el, struct point *p)
{
	struct edge *e = el->el_edge;
	struct site *topsite = e->reg[1];
	int right_of_site = p->x > topsite->coord.x;

	if (right_of_site && el->el_pm == le) return (1);
	if (!right_of_site && el->el_pm == re) return (0);

	int above;
	if (e->a == 1.0) {
		float dyp = p->y - topsite->coord.y;
		float dxp = p->x - topsite->coord.x;
		int fast = 0;
		if ((!right_of_site & (e->b < 0.0)) | (right_of_site & (e->b >= 0.0))) {
			above = dyp >= e->b * dxp;
			fast = above;
		} else {
			above = p->x + p->y * e->b > e-> c;
			if (e->b < 0.0) above = !above;
			if (!above) fast = 1;
		}
		if (!fast) {
			float dxs = topsite->coord.x - (e->reg[0])->coord.x;
			above = e->b * (dxp * dxp - dyp * dyp) <
				dxs * dyp * (1.0 + 2.0 * dxp / dxs + e->b * e->b);
			if (e->b < 0.0) above = !above;
		}
	} else { /*e->b==1.0 */
		float yl = e->c - e->a * p->x;
		float t1 = p->y - yl;
		float t2 = p->x - topsite->coord.x;
		float t3 = yl - topsite->coord.y;
		above = t1 * t1 > t2 * t2 + t3 * t3;
	}
	return (el->el_pm == le ? above : !above);
}

struct half_edge * el_leftbnd(struct el *el, struct point *p)
{
	/* Use hash table to get close to desired halfedge */
	/* use the hash function to find the place in the hash map that this Halfedge should be */

	int bucket = (int)((p->x - el->xmin) / el->deltax * el->hashsize);
	/* make sure that the bucket position in within the range of the hash array */
	if (bucket < 0) bucket = 0;
	if (bucket >= el->hashsize) bucket = el->hashsize - 1;

	struct half_edge *he = el_gethash(el, bucket);
	/* if the HE isn't found, search backwards and forwards in the hash
	 * map for the first non-null entry */
	if (!he) {
		for (int i = 1; 1; i++) {
			if ((he = el_gethash(el, bucket - i)))
				break;
			if ((he = el_gethash(el, bucket + i)))
				break;
		}
	}
	/* Now search linear list of halfedges for the correct one */
	if (he == el->leftend || (he != el->rightend && right_of(he, p))) {
		/* keep going right on the list until either the end is reached,
		 * or you find the 1st edge which the point */
		do {
			he = he->el_right;
		} while (he != el->rightend && right_of(he, p));
		/* isn't to the right of */
		he = he->el_left;
	} else
		/* if the point is to the left of the Halfedge, then search
		 * left for the HE just to the left of the point */
		do {
			he = he->el_left;
		} while (he != el->leftend && !right_of(he, p));

	/* Update hash table and reference counts */
	if (bucket > 0 && bucket < el->hashsize - 1) {
		if (el->hash[bucket])
			el->hash[bucket]->el_refcnt--;
		el->hash[bucket] = he;
		el->hash[bucket]->el_refcnt++;
	};
	return he;
}

/* This delete routine can't reclaim node, since pointers from hash
table may be present.   */
void el_delete(struct half_edge *he)
{
	(he->el_left)->el_right = he->el_right;
	(he->el_right)->el_left = he->el_left;
	he->el_edge = (struct edge *)DELETED;
}

void geom_init(struct geom *geom, struct free_node_al **b, int sqrt_nsites, float minX, float minY, float maxX, float maxY, float min_dist)
{
	geom->b.b = b;
	geom->b.sqrt_nsites = sqrt_nsites;
	geom->edges = NULL;
	geom->efl.head = NULL;
	geom->efl.nodesize = sizeof(struct edge);
	geom->sfl.head = NULL;
	geom->sfl.nodesize = sizeof(struct site);
	geom->nvertices = 0;
	geom->nedges = 0;
	geom->borderMinX = minX;
	geom->borderMinY = minY;
	geom->borderMaxX = maxX;
	geom->borderMaxY = maxY;
	geom->min_dist = min_dist;
}

static void deref(struct geom *geom, struct site *v)
{
	if (!--v->refcnt) makefree((struct free_node *)v, &geom->sfl);
}

static void ref(struct site *v)
{
	v->refcnt++;
}

struct edge *bisect(struct geom *geom, struct site *s1, struct site *s2)
{
	struct edge *newedge = getfree(&geom->b, &geom->efl);

	/* store the sites that this edge is bisecting */
	newedge->reg[0] = s1;
	newedge->reg[1] = s2;
	ref(s1);
	ref(s2);
	/* to begin with, there are no endpoints on the bisector
	 * - it goes to infinity */
	newedge->ep[0] = NULL;
	newedge->ep[1] = NULL;
	
	/* get the difference in x dist between the sites */
	float dx = s2->coord.x - s1->coord.x;
	float dy = s2->coord.y - s1->coord.y;
	/* make sure that the difference in positive */
	float adx = dx > 0 ? dx : -dx;					
	float ady = dy > 0 ? dy : -dy;
	/* get the slope of the line */
	newedge->c = (float)(s1->coord.x * dx + s1->coord.y * dy +
			(dx * dx + dy * dy) * 0.5);

	if (adx > ady) {
		newedge->a = 1.0;
		newedge->b = dy / dx;
		newedge->c /= dx;
		/* set formula of line, with x fixed to 1 */
	} else {
		newedge->b = 1.0;
		newedge->a = dx / dy;
		newedge->c /= dy;
		/* set formula of line, with y fixed to 1 */
	}

	newedge->edgenbr = geom->nedges;

	//printf("\nbisect(%d) ((%f,%f) and (%f,%f)",nedges,s1->coord.x,s1->coord.y,s2->coord.x,s2->coord.y);

	geom->nedges += 1;
	return newedge;
}

/**
 * create a new site where the Halfedges el1 and el2 intersect - note that the
 * point in the argument list is not used, don't know why it's there
 */
struct site *intersect(struct geom *geom, struct half_edge *el1, struct half_edge *el2)
{
	struct edge *e1 = el1->el_edge, *e2 = el2->el_edge;
	if (!e1 || !e2)
		return NULL;

	//if the two edges bisect the same parent, return null
	if (e1->reg[1] == e2->reg[1])
		return NULL;

	float d = e1->a * e2->b - e1->b * e2->a;
	if (-1.0e-10 < d && d < 1.0e-10)
		return NULL;

	float xint = (e1->c * e2->b - e2->c * e1->b) / d;
	float yint = (e2->c * e1->a - e1->c * e2->a) / d;

	struct edge *e;
	struct half_edge *el;
	if ((e1->reg[1]->coord.y < e2->reg[1]->coord.y) ||
	    (e1->reg[1]->coord.y == e2->reg[1]->coord.y &&
	     e1->reg[1]->coord.x < e2->reg[1]->coord.x)) {
		el = el1;
		e = e1;
	} else {
		el = el2;
		e = e2;
	}

	int right_of_site = xint >= e->reg[1]->coord.x;
	if ((right_of_site && el->el_pm == le) || (!right_of_site && el->el_pm == re))
		return NULL;

	/* create a new site at the point of intersection - this is a new
	 * vector event waiting to happen */
	struct site *v = getfree(&geom->b, &geom->sfl);
	v->refcnt = 0;
	v->coord.x = xint;
	v->coord.y = yint;
	return v;
}

void clip_line(struct geom *geom, struct edge *e)
{
	float x1 = e->reg[0]->coord.x;
	float x2 = e->reg[1]->coord.x;
	float y1 = e->reg[0]->coord.y;
	float y2 = e->reg[1]->coord.y;

	/* if the distance between the two points this line was created from is less than
	 * the square root of 2, then ignore it */
	if (sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) < geom->min_dist)
		return;
	float pxmin = geom->borderMinX;
	float pxmax = geom->borderMaxX;
	float pymin = geom->borderMinY;
	float pymax = geom->borderMaxY;

	struct site *s1, *s2;
	if (e->a == 1.0 && e ->b >= 0.0) {
		s1 = e->ep[1];
		s2 = e->ep[0];
	} else {
		s1 = e->ep[0];
		s2 = e->ep[1];
	}

	if (e->a == 1.0) {
		y1 = pymin;
		if (s1 && s1->coord.y > pymin) y1 = s1->coord.y;
		if (y1 > pymax) y1 = pymax;
		x1 = e->c - e->b * y1;
		y2 = pymax;
		if (s2 && s2->coord.y < pymax) y2 = s2->coord.y;
		if (y2 < pymin) y2 = pymin;
		x2 = (e->c) - (e->b) * y2;
		if (((x1 > pxmax) & (x2 > pxmax)) | ((x1 < pxmin) & (x2 < pxmin))) return;
		if (x1 > pxmax) {
			x1 = pxmax;
			y1 = (e->c - x1) / e->b;
		}
		if (x1 < pxmin) {
			x1 = pxmin;
			y1 = (e->c - x1) / e->b;
		}
		if (x2 > pxmax) {
			x2 = pxmax;
			y2 = (e->c - x2) / e->b;
		}
		if (x2 < pxmin) {
			x2 = pxmin;
			y2 = (e->c - x2) / e->b;
		}
	} else {
		x1 = pxmin;
		if (s1 && s1->coord.x > pxmin) x1 = s1->coord.x;
		if (x1 > pxmax) x1 = pxmax;
		y1 = e->c - e->a * x1;
		x2 = pxmax;
		if (s2 && s2->coord.x < pxmax) x2 = s2->coord.x;
		if (x2 < pxmin) x2 = pxmin;
		y2 = e->c - e->a * x2;
		if (((y1 > pymax) & (y2 > pymax)) | ((y1 < pymin) & (y2 < pymin))) return;
		if (y1 > pymax) {
			y1 = pymax;
			x1 = (e->c - y1) / e->a;
		}
		if (y1 < pymin) {
			y1 = pymin;
			x1 = (e->c - y1) / e->a;
		}
		if (y2 > pymax) {
			y2 = pymax;
			x2 = (e->c - y2) / e->a;
		}
		if (y2 < pymin) {
			y2 = pymin;
			x2 = (e->c - y2) / e->a;
		}
	}

	//printf("\nPushing line (%f,%f,%f,%f)",x1,y1,x2,y2);
	/* new line */
	struct graph_edge* ge = malloc(sizeof(struct graph_edge));
	ge->next = geom->edges;
	geom->edges = ge;
	ge->x1 = x1;
	ge->y1 = y1;
	ge->x2 = x2;
	ge->y2 = y2;
}

void endpoint(struct geom *geom, struct edge *e, int lr, struct site * s)
{
	e->ep[lr] = s;
	ref(s);
	if (!e->ep[re - lr]) return;

	clip_line(geom, e);

	deref(geom, e->reg[le]);
	deref(geom, e->reg[re]);
	makefree((struct free_node*)e, &geom->efl);
}

float dist(struct site *s, struct site *t)
{
	float dx = s->coord.x - t->coord.x;
	float dy = s->coord.y - t->coord.y;
	return (float)(sqrt(dx * dx + dy * dy));
}

int pq_bucket(struct pq *pq, struct half_edge *he)
{
	int bucket = (int)((he->ystar - pq->ymin) / pq->deltay * pq->pq_hashsize);
	if (bucket < 0) bucket = 0;
	if (bucket >= pq->pq_hashsize) bucket = pq->pq_hashsize - 1 ;
	if (bucket < pq->pq_min) pq->pq_min = bucket;
	return bucket;
}

/* push the Halfedge into the ordered linked list of vertices */
void pq_insert(struct pq *pq, struct half_edge *he, struct site * v, float offset)
{
	struct half_edge *last, *next;

	he->vertex = v;
	ref(v);
	he->ystar = (float)(v->coord.y + offset);
	last = &pq->pq_hash[pq_bucket(pq, he)];
	while ((next = last->pq_next) &&
	       (he->ystar  > next->ystar  ||
		(he->ystar == next->ystar && v->coord.x > next->vertex->coord.x))) {
		last = next;
	}
	he->pq_next = last->pq_next;
	last->pq_next = he;
	pq->pq_count++;
}

/* remove the Halfedge from the list of vertices */
void pq_delete(struct pq *pq, struct geom *geom, struct half_edge *he)
{
	if (he->vertex) {
		struct half_edge *last = &pq->pq_hash[pq_bucket(pq, he)];
		while (last->pq_next != he)
			last = last->pq_next;

		last->pq_next = he->pq_next;
		pq->pq_count--;
		deref(geom, he->vertex);
		he->vertex = NULL;
	}
}

int pq_empty(struct pq *pq)
{
	return (pq->pq_count == 0);
}

struct point pq_min(struct pq *pq)
{
	while (!pq->pq_hash[pq->pq_min].pq_next)
		pq->pq_min++;
	struct point answer;
	answer.x = pq->pq_hash[pq->pq_min].pq_next->vertex->coord.x;
	answer.y = pq->pq_hash[pq->pq_min].pq_next->ystar;
	return answer;
}

struct half_edge *pq_extractmin(struct pq *pq)
{
	struct half_edge *curr = pq->pq_hash[pq->pq_min].pq_next;
	pq->pq_hash[pq->pq_min].pq_next = curr->pq_next;
	pq->pq_count--;
	return curr;
}

int pq_init(struct pq *pq, int sqrt_nsites, int ymin, int deltay)
{
	pq->pq_count = 0;
	pq->pq_min = 0;
	pq->pq_hashsize = 4 * sqrt_nsites;
	pq->pq_hash = malloc(pq->pq_hashsize * sizeof *pq->pq_hash);
	pq->ymin = ymin;
	pq->deltay = deltay;

	if (!pq->pq_hash)
		return 0;

	for (int i = 0; i < pq->pq_hashsize; i++)
		pq->pq_hash[i].pq_next = NULL;

	return 1;
}

void pq_free(struct pq *pq)
{
	if (pq->pq_hash)
		free(pq->pq_hash);
}

void voronoi_free(struct graph_edge *e)
{
	struct graph_edge* gec = e, *gep = NULL;
	while (gec) {
		gep = gec;
		gec = gec->next;
		free(gep);
	}
}

int scomp(const void *p1, const void *p2)
{
	struct point *s1 = (struct point *)p1, *s2 = (struct point *)p2;
	if (s1->y < s2->y) return -1;
	if (s1->y > s2->y) return 1;
	if (s1->x < s2->x) return -1;
	if (s1->x > s2->x) return 1;
	return 0;
}

static void out_site(struct site *s)
{
}

#define RIGHT_REG(he) (!he->el_edge ? bottomsite : (he->el_pm == le ? he->el_edge->reg[re] : he->el_edge->reg[le]));

/* implicit parameters: nsites, sqrt_nsites, xmin, xmax, ymin, ymax,
deltax, deltay (can all be estimates).
Performance suffers if they are wrong; better to make nsites,
deltax, and deltay too big than too small.  (?) */
struct graph_edge *voronoi(struct point *xy, int nsites, float min_x, float max_x, float min_y, float max_y, float min_dist)
{
	/* state reqeuired for Fortune's algorithm */
	struct pq pq;
	struct el el;
	struct geom geom;

	//struct free_node_al *mem = malloc(sizeof(struct free_node_al));
	struct free_node_al *mem = NULL;
	struct site *sites = malloc(nsites * sizeof(*sites));

	if (!sites) return NULL;

	float xmin = xy[0].x, xmax = xy[0].x;
	float ymin = xy[0].y, ymax = xy[0].y;

	for (int i = 0; i < nsites; i++) {
		sites[i].coord = xy[i];
		sites[i].sitenbr = i;
		sites[i].refcnt = 0;

		if (xy[i].x < xmin)
			xmin = xy[i].x;
		else if (xy[i].x > xmax)
			xmax = xy[i].x;

		if (xy[i].y < ymin)
			ymin = xy[i].y;
		else if (xy[i].y > ymax)
			ymax = xy[i].y;

		//printf("\n%f %f\n",xValues[i],yValues[i]);
	}
	float deltay = ymax - ymin;
	float deltax = xmax - xmin;

	qsort(sites, nsites, sizeof(*sites), scomp);

	float temp = 0;
	if (min_x > max_x) {
		temp = min_x;
		min_x = max_x;
		max_x = temp;
	}
	if (min_y > max_y) {
		temp = min_y;
		min_y = max_y;
		max_y = temp;
	}

	int sqrt_nsites = (int)sqrt((float)nsites + 4);
	geom_init(&geom, &mem, sqrt_nsites, min_x, min_y, max_x, max_y, min_dist);

	int siteidx = 0;

	pq_init(&pq, sqrt_nsites, ymin, deltay);
	struct site *bottomsite = sites + siteidx++;
	out_site(bottomsite);

	if (!el_init(&el, &mem, sqrt_nsites, xmin, deltax)) return NULL;

	struct site *newsite = (siteidx < nsites) ? sites + siteidx++ : NULL;
	while (1) {
		struct point newintstar;
		if (!pq_empty(&pq))
			newintstar = pq_min(&pq);

		/* if the lowest site has a smaller y value than the lowest vector
		 * intersection, process the site otherwise process the vector
		 * intersection */

		if (newsite && (pq_empty(&pq) || newsite->coord.y < newintstar.y
				|| (newsite->coord.y == newintstar.y && newsite->coord.x < newintstar.x))) {
			/* new site is smallest - this is a site event*/
			out_site(newsite); /* output the site */
			/* get the first half_edge to the LEFT & RIGHT of site */
			struct half_edge *lbnd = el_leftbnd(&el, &(newsite->coord));
			struct half_edge *rbnd = lbnd->el_right;
			/* if this halfedge has no edge, bot = bottom site (whatever that is) */
			struct site *bot = RIGHT_REG(lbnd);
			/* create a new edge that bisects */
			struct edge *e = bisect(&geom, bot, newsite);
			/* create new half edge */
			struct half_edge *bisector = he_create(&el, e, le);
			/* insert new bisector edge */
			el_insert(lbnd, bisector);
			/* if the new bisector intersects with the left edge,
			 * remove the left edge's vertex, and put in the new one */
			struct site *p;
			if ((p = intersect(&geom, lbnd, bisector))) {
				pq_delete(&pq, &geom, lbnd);
				pq_insert(&pq, lbnd, p, dist(p, newsite));
			}
			lbnd = bisector;
			/* create new half edge */
			bisector = he_create(&el, e, re);
			/* insert the new HE to the right of the original
			 * bisector earlier in the IF stmt */
			el_insert(lbnd, bisector);

			/* if this new bisector intersects with the */
			/* push the HE into the ordered linked list of vertices */
			if ((p = intersect(&geom, bisector, rbnd)))
				pq_insert(&pq, bisector, p, dist(p, newsite));
			/* allocate newsite if sites left */
			newsite = (siteidx < nsites) ? sites + siteidx++ : NULL;
		} else if (!pq_empty(&pq)) { /* intersection is smallest - this is a vector event */
			/* pop the Halfedge with the lowest vector off the ordered list of vectors */
			struct half_edge *lbnd = pq_extractmin(&pq); 
			/* get the Halfedge to the left of the above HE */
			struct half_edge *llbnd = lbnd->el_left; 
			/* get the Halfedge to the right of the above HE */
			struct half_edge *rbnd = lbnd->el_right; 
			/* get the Halfedge to the right of the HE to the
			 * right of the lowest HE */
			struct half_edge *rrbnd = rbnd->el_right; 
			/* get the site to the left of the left HE which it bisects */
			struct site *bot = !lbnd->el_edge ? bottomsite : (lbnd->el_pm == le ?
				lbnd->el_edge->reg[le] : lbnd->el_edge->reg[re]);

			/* get the site to the right of the right HE which it bisects */
			struct site *top = RIGHT_REG(rbnd);

			/* output the triple of sites, stating that a
			 * circle goes through them */
			//out_triple(bot, top, RIGHT_REG(lbnd));
			/* get the vertex that caused this event */
			struct site *v = lbnd->vertex;
			/* makevertex - set the vertex number - couldn't do
			 * this earlier since we didn't know when it would be processed */
			v->sitenbr = geom.nvertices++;
			/* set the endpoint of the left & right Halfedge to be this vector */
			endpoint(&geom, lbnd->el_edge, lbnd->el_pm, v);
			endpoint(&geom, rbnd->el_edge, rbnd->el_pm, v);
			/* mark the lowest HE for deletion - can't delete yet because there might be pointers to it in Hash Map */
			el_delete(lbnd);
			/* remove all vertex events to do with the  right HE */
			pq_delete(&pq, &geom, rbnd);
			/* mark the right HE for deletion - can't delete yet because there might be pointers to it in Hash Map */
			el_delete(rbnd);
			int pm = le; /* set the pm variable to zero */
			if (bot->coord.y > top->coord.y) {	//if the site to the left of the event is higher than the site
				//to the right of it, then swap them and set the 'pm' variable to 1
				struct site *temp = bot;
				bot = top;
				top = temp;
				pm = re;
			}
			/* create an edge (or line) that is between the two sites. This creates 
			 * the formula of the line, and assigns a line number to it */
			struct edge *e = bisect(&geom, bot, top);
			/* create a HE from the edge 'e', and make it point to that edge with its el_edge field */
			struct half_edge *bisector = he_create(&el, e, pm);
			/* insert the new bisector to the right of the left HE */
			el_insert(llbnd, bisector);
			/* set one endpoint to the new edge to be the vector point 'v'. */
			endpoint(&geom, e, re - pm, v);
			/* If the site to the left of this bisector is higher than the right
			 * site, then this endpoint is put in position 0; otherwise in pos 1 */
			deref(&geom, v); /* delete the vector 'v' */

			/* if left HE and the new bisector don't intersect, then delete the left HE, and reinsert it */
			struct site *p;
			if ((p = intersect(&geom, llbnd, bisector))) {
				pq_delete(&pq, &geom, llbnd);
				pq_insert(&pq, llbnd, p, dist(p, bot));
			}

			/* if right HE and the new bisector don't intersect, then reinsert it */
			if ((p = intersect(&geom, bisector, rrbnd)))
				pq_insert(&pq, bisector, p, dist(p, bot));
		} else
			break;
	}

	for (struct half_edge *lbnd = el.leftend->el_right; lbnd != el.rightend; lbnd = lbnd->el_right)
		clip_line(&geom, lbnd->el_edge);

	/* cleanup the remaining memory */
	free(sites);

	struct free_node_al *current = mem, *prev;
	while (current) {
		prev = current;
		current = current->next;
		free(prev->memory);
		free(prev);
	}

	el_free(&el);
	pq_free(&pq);

	return geom.edges;
}

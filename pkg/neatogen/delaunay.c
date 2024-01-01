/*************************************************************************
 * Copyright (c) 2011 AT&T Intellectual Property 
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * which accompanies this distribution, and is available at
 * http://www.eclipse.org/legal/epl-v10.html
 *
 * Contributors: Details at https://graphviz.org
 *************************************************************************/

#include "config.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cgraph/alloc.h>
#include <cgraph/cgraph.h>     /* for agerr() and friends */
#include <neatogen/delaunay.h>
#include <common/memory.h>

static char* err = "Graphviz built without any triangulation library\n";
int* get_triangles (double *x, int n, int* tris)
{
    agerr(AGERR, "get_triangles: %s\n", err);
    return 0;
}
static v_data *delaunay_triangulation(double *x, double *y, int n) {
    agerr(AGERR, "delaunay_triangulation: %s\n", err);
    return 0;
}
int *delaunay_tri(double *x, double *y, int n, int* nedges)
{
    agerr(AGERR, "delaunay_tri: %s\n", err);
    return 0;
}
surface_t*
mkSurface (double *x, double *y, int n, int* segs, int nsegs)
{
    agerr(AGERR, "mkSurface: %s\n", err);
    return 0;
}
void
freeSurface (surface_t* s)
{
    agerr (AGERR, "freeSurface: %s\n", err);
}

static void remove_edge(v_data * graph, int source, int dest)
{
    int i;
    for (i = 1; i < graph[source].nedges; i++) {
	if (graph[source].edges[i] == dest) {
	    graph[source].edges[i] = graph[source].edges[--graph[source].nedges];
	    break;
	}
    }
}

v_data *UG_graph(double *x, double *y, int n) {
    v_data *delaunay;
    int i;
    double dist_ij, dist_ik, dist_jk, x_i, y_i, x_j, y_j;
    int j, k, neighbor_j, neighbor_k;

    if (n == 2) {
	int *edges = gv_calloc(4, sizeof(int));
	delaunay = gv_calloc(n, sizeof(v_data));
	delaunay[0].ewgts = NULL;
	delaunay[0].edges = edges;
	delaunay[0].nedges = 2;
	delaunay[0].edges[0] = 0;
	delaunay[0].edges[1] = 1;
	delaunay[1].edges = edges + 2;
	delaunay[1].ewgts = NULL;
	delaunay[1].nedges = 2;
	delaunay[1].edges[0] = 1;
	delaunay[1].edges[1] = 0;
	return delaunay;
    } else if (n == 1) {
	int *edges = gv_calloc(1, sizeof(int));
	delaunay = gv_calloc(n, sizeof(v_data));
	delaunay[0].ewgts = NULL;
	delaunay[0].edges = edges;
	delaunay[0].nedges = 1;
	delaunay[0].edges[0] = 0;
	return delaunay;
    }

    delaunay = delaunay_triangulation(x, y, n);

    // remove all edges v-u if there is w, neighbor of u or v, that is closer to both u and v than dist(u,v)
    for (i = 0; i < n; i++) {
        x_i = x[i];
        y_i = y[i];
        for (j = 1; j < delaunay[i].nedges;) {
            neighbor_j = delaunay[i].edges[j];
            x_j = x[neighbor_j];
            y_j = y[neighbor_j];
            dist_ij = (x_j - x_i) * (x_j - x_i) + (y_j - y_i) * (y_j - y_i);
            // now look at i'th neighbors to see whether there is a node in the "forbidden region"
            // we will also go through neighbor_j's neighbors when we traverse the edge from its other side
            bool removed = false;
            for (k = 1; k < delaunay[i].nedges && !removed; k++) {
                neighbor_k = delaunay[i].edges[k];
                dist_ik = (x[neighbor_k] - x_i) * (x[neighbor_k] - x_i) +
                    (y[neighbor_k] - y_i) * (y[neighbor_k] - y_i);
                if (dist_ik < dist_ij) {
                    dist_jk = (x[neighbor_k] - x_j) * (x[neighbor_k] - x_j) +
                        (y[neighbor_k] - y_j) * (y[neighbor_k] - y_j);
                    if (dist_jk < dist_ij) {
                        // remove the edge beteween i and neighbor j
                        delaunay[i].edges[j] = delaunay[i].edges[--delaunay[i].nedges];
                        remove_edge(delaunay, neighbor_j, i);
                        removed = true;
                    }
                }
            }
            if (!removed) {
                j++;
            }
        }
    }
    return delaunay;
}

void freeGraph (v_data * graph)
{
    if (graph != NULL) {
	free(graph[0].edges);
	free(graph[0].ewgts);
	free(graph);
    }
}

void freeGraphData(vtx_data * graph)
{
    if (graph != NULL) {
	free(graph[0].edges);
	free(graph[0].ewgts);
#ifdef DIGCOLA
	free(graph[0].edists);
#endif
	free(graph);
    }
}


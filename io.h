#ifndef _io_h_
#define _io_h_

#include "graph.h"

/*
*	Reads graph from stdin specified in 
*   PACE Challenge 2019 format.
*/
void read_graph(Graph & graph);

/*
*	Writes vertex cover to stdout (temporarily) as specified
*	in PACE Challenge 2019 solution format.
*/
void write_vertex_cover(VertexCover const & vc, Number const num_vertices);

/*
 *  Removes duplicate entries from adjacency list.
 *  To be used after reading in graph.
 */
void remove_duplicates(AdjacencyList & vertices);

#endif
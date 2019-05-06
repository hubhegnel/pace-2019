#include <iostream>
#include <cstring>
#include <string>
#include "graph.h"
#include "io.h"
#include <algorithm>

using namespace std;

int main(int argc, char *argv[])
{
    Graph graph;
    read_graph(graph);

    vector<Vertex> branch_set;
    ChangeVector empty;
    Number zero = 0;
    
    bool changed = true;
    while (changed)
    {
        graph.apply_reduction_rules_and_lower_bounds(empty, zero);
        changed = graph.apply_dominance_rule(empty, zero);
    }

    graph.compute_modulator_to_bipartite(branch_set);

    graph.branch(0, branch_set);
    graph.undo_changes(empty);
    if (graph.validate_vertex_cover())
    {
        write_vertex_cover(graph._best_vc, graph._num_vertices);
    }
    else
    {
        while (true)
        {
            //Timeout
        }
    }
    
}

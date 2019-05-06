#include <iostream>
#include <sstream>
#include <fstream>
#include <cstring>
#include <string>
#include <iterator>
#include <algorithm>
#include "io.h"

using namespace std;

void read_graph(Graph & graph)
{
    Number const infty = numeric_limits<Number>::max();
    std::ios::sync_with_stdio(false);

    string temp;

    //skip first batch of comments
    while (getline(cin, temp) and temp.at(0) == 'c')
    {
    }

    istringstream iss(temp);
    vector<string> results(istream_iterator<string>{iss}, istream_iterator<string>());

    //problem description line should come now
    if (results[0] != "p" or (results[1] != "td" and results[1] != "edge" and results[1] != "tw"))
    {
        cerr << "Invalid format" << endl;
    }
    else
    {
        graph._num_vertices = stoul(results[2]);
        graph._num_orig_vertices = graph._num_vertices;
        graph._num_edges = stoul(results[3]);
    }

    //initialize vertex states
    graph._state.resize(graph._num_vertices, undecided);

    //initialize Adjacency Lists of Graph
    graph._degrees.clear();
    graph._degrees.assign(graph._num_vertices,infty);

    graph._adj.resize(graph._num_vertices);

    for (Number i = 0; i < graph._num_edges; ++i)
    {
        cin >> temp;
        if (temp == "c")
        {
            //skip comment lines
            --i;
            getline(cin, temp);
        }
        else
        {
            //read edges
            Vertex v;
            Vertex w;
            v = stoul(temp);

            cin >> w;

            //-1 since vertices are numbered 1..n in PACE format
            //but internally we use 0..n-1
            if (v != w)
            {
                graph._adj[v - 1].push_back(w - 1);
                graph._adj[w - 1].push_back(v - 1);
            }
            else
            {
                graph._adj[v - 1].push_back(v - 1);
                graph._state[v - 1] = taken;
                cerr << "loop " << v << endl;
            }
        }
    }

    //remove duplicate edges
    graph._num_edges = 0;
    for (Vertex v = 0; v < graph._num_vertices; ++v)
    {
        remove_duplicates(graph._adj[v]);
        graph._num_edges += graph._adj[v].size();
    }
    graph._num_edges /= 2;

    graph._best_size = graph._num_vertices;

    graph._skip_degree_reduction_once = false;

    graph._folds.clear();
}

void write_vertex_cover(VertexCover const & vc, Number const num_vertices)
{
    cout << "s vc " << num_vertices << " " << vc.size() << '\n';

    for (auto it = vc.begin(); it != vc.end(); ++it)
    {
        cout << *it + 1 << '\n';
    }
}
#include <string>
#include <iostream>
#include <algorithm>
#include <queue>
#include "graph.h"
#include "io.h"


using namespace std;

Number const INFTY = numeric_limits<Number>::max();


void remove_duplicates(vector<Vertex> & vertices)
{
    if (not vertices.empty())
    {
        sort(vertices.begin(), vertices.end());
        vector<Vertex> copy(vertices);
        vertices.clear();
        vertices.push_back(copy[0]);
        Vertex read_last = copy[0];
        for (Number i = 1; i < copy.size(); ++i)
        {
            if (copy[i] != read_last)
            {
                vertices.push_back(copy[i]);
                read_last = copy[i];
            }
        }
    }
}

void Graph::remove_duplicates_with_degree_update(vector<Vertex> & vertices, vector<Vertex>  & process_stack)
{
    if (not vertices.empty())
    {
        sort(vertices.begin(), vertices.end());
        vector<Vertex> copy(vertices);
        vertices.clear();
        vertices.push_back(copy[0]);
        Vertex read_last = copy[0];
        for (Number i = 1; i < copy.size(); ++i)
        {
            if (copy[i] != read_last)
            {
                vertices.push_back(copy[i]);
                read_last = copy[i];
            }
            else
            {
                if (_degrees[read_last] == INFTY)
                {
                    _degrees[read_last] = degree(read_last) - 1;
                }
                else
                {
                    _degrees[read_last] = _degrees[read_last] - 1;
                }

                if (_degrees[read_last] < 3)
                {
                    process_stack.push_back(read_last);
                }
            }
        }
    }
}



StatusVector copy_vec(StatusVector & states)
{
    StatusVector new_vec(states.size());
    copy(states.begin(), states.end(), new_vec.begin());
    return new_vec;
}



Change::Change(Vertex v)
{
    _type = state_change;
    _vertex = v;
}

Change::Change(ChangeType type, Vertex v)
{
    _type = type;
    _vertex = v;
}



FoldData::FoldData(Vertex new_vertex, Vertex when_discarded, vector<Vertex> & when_taken)
{
    _new_vertex = new_vertex;
    _when_discarded = when_discarded;
    _when_taken = when_taken;
}



void Graph::add_vertex(vector<Vertex> & neighbors, ChangeVector & changes)
{
    //new vertex will be the last, so has id num_vertices
    _num_edges += neighbors.size();
    _state.push_back(undecided);
    _adj.push_back(neighbors);

    //also add new vertex into adjacency lists of neighbors
    for (Number i = 0; i < neighbors.size(); ++i)
    {
        _adj[neighbors[i]].push_back(_num_vertices);
    }
    changes.push_back(Change(new_vertex_due_to_fold, _num_vertices));
    ++_num_vertices;
}

void Graph::delete_last_vertex()
{
    --_num_vertices; //now last vertex has id num_vertices
    if (_num_vertices < _num_orig_vertices)
    {
        throw "Deleted too many vertices!";
    }
    _num_edges -= _adj[_num_vertices].size();

    //delete vertex in adjacency lists of neighbors
    AdjacencyList & neighbors = _adj[_num_vertices];
    for (Number i = 0; i < neighbors.size(); ++i)
    {
        if (_adj[neighbors[i]].back() != _num_vertices)
        {
            throw "Deleted vertex not last in neighbors adjacency list!";
        }
        _adj[neighbors[i]].pop_back();
    }
    _state.pop_back();
    _adj.pop_back();
}

void Graph::change_state(Vertex v, Status new_state, ChangeVector & changes)
{
    _state[v] = new_state;
    changes.push_back(v);
}

bool Graph::edge_exists(Vertex v, Vertex w)
{
    if (_state[v] != undecided or _state[w] != undecided)
    {
        return false;
    }

    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        if (_adj[v][i] == w)
        {
            return true;
        }
    }
    return false;
}

void Graph::add_undecided_neighbors_to_vector(Vertex v, vector<Vertex> & vertices)
{
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        Vertex neighbor = _adj[v][i];
        if (_state[neighbor] == undecided)
        {
            vertices.push_back(neighbor);
        }
    }
}

vector<ConnectedComponent> Graph::compute_connected_components_BFS(vector<BipState> & bipartite_states,
                                                                   vector<bool> & is_bipartite)
{
    vector<ConnectedComponent> cc;

    Number i = 0;
    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        if (bipartite_states[v] == unvisited)
        {
            cc.resize(i + 1);
            cc[i].push_back(v);
            bipartite_states[v] = left_side;
            is_bipartite.push_back(BFS_and_bipartite_test(v,  bipartite_states, cc, i));
            ++i;
        }
    }

    return cc;

}
bool Graph::BFS_and_bipartite_test(Vertex v,
                                   vector<BipState> & bipartite_states,
                                   vector<ConnectedComponent> & cc,
                                   Number i)
{
    bool is_bipartite = true;

    Vertex current, neighbor;
    queue<Vertex> BFSqueue;
    BFSqueue.push(v);
    bool left = true;
    while (not BFSqueue.empty())
    {
        current = BFSqueue.front();
        BFSqueue.pop();

        for (Number j = 0; j < _adj[current].size(); ++j)
        {
            neighbor = _adj[current][j];

            switch(bipartite_states[neighbor])
            {
                case unvisited:
                    if (bipartite_states[current] == left_side)
                    {
                        bipartite_states[neighbor] = right_side;
                    }
                    else if (bipartite_states[current] == right_side)
                    {
                        bipartite_states[neighbor] = left_side;
                    }
                    else
                    {
                        if (left)
                        {
                            bipartite_states[neighbor] = left_side;
                            left = not left;
                        }
                        else
                        {
                            bipartite_states[neighbor] = right_side;
                            left = not left;
                        }
                    }
                    cc[i].push_back(neighbor);
                    BFSqueue.push(neighbor);
                    break;
                case left_side:
                    if (bipartite_states[current] == left_side)
                    {
                        //conflict, so current is destructive
                        bipartite_states[current] = destructive;
                        is_bipartite = false;
                    }
                    break;
                case right_side:
                    if (bipartite_states[current] == right_side)
                    {
                        //conflict, so current is destructive
                        bipartite_states[current] = destructive;
                        is_bipartite = false;
                    }
                    break;
                default:
                    break;
            }
        }
    }

    return is_bipartite;
}



Degree Graph::degree(Vertex v)
{
    Degree deg = 0;
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        if (_state[_adj[v][i]] == undecided)
        {
            ++deg;
        }
    }
    return deg;
}


void Graph::recompute_degrees()
{
    _degrees.assign(_num_vertices, INFTY);
    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        _degrees[v] = degree(v);
    }
}

double Graph::small_degree_neighbors(Vertex v, Degree small_degree)
{
    double weight = 1.25;
    double count = 0;
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        Vertex neighbor = _adj[v][i];
        if (_state[neighbor] == undecided)
        {
            if (_degrees[neighbor] <= small_degree)
            {
                count += weight;
            }
            else
            {
                ++count;
            }
        }
    }
    return count;
}

Number Graph::branch(Number lb, vector<Vertex> branch_set)
{
    bool branch_cc = true; //if set to true, cc are considered
    ChangeVector changes;

    //first apply reduction rules
    apply_reduction_rules_and_lower_bounds(changes, lb);

    //compare current solution size with upper bound
    Number cur_solution_size = vertex_cover_size();

    if (cur_solution_size + lb >= _best_size)
    {
        //prune branch
        undo_changes(changes);
        return _num_vertices;
    }

    if (branch_cc)
    {
        //check, if graph is still connected and compute modular to bipartite set

        vector<BipState> bipartite_states(_num_vertices, non_existent);
        vector<bool> is_bipartite;
        //set vertices in G[undecided] to unvisited
        for (Vertex v = 0; v < _num_vertices; ++v)
        {
            if (_state[v] == undecided)
            {
                bipartite_states[v] = unvisited;
            }
        }

        vector<ConnectedComponent> all_cc = compute_connected_components_BFS(bipartite_states, is_bipartite);

        Number count_undecided = 0;
        Number count_destructive = 0;
        for (Number i = 0 ; i < bipartite_states.size(); ++i)
        {
            if (bipartite_states[i] != non_existent)
            {
                count_undecided++;
            }
            if (bipartite_states[i] == destructive)
            {
                count_destructive++;
            }
        }

        if (static_cast<double>(count_destructive) / static_cast<double>(count_undecided) <= 0.2-0.02*count_destructive
            and (branch_set.empty() or count_destructive < branch_set.size()))
        {
            branch_set.clear();
            for (Vertex v = 0; v < _num_vertices; ++v)
            {
                if (bipartite_states[v] == destructive)
                {
                    branch_set.push_back(v);
                }
            }
            sort(branch_set.begin(), branch_set.end(), [this](Vertex v, Vertex w)
            {
                return _degrees[v] < _degrees[w];
            });
        }

        if ((is_bipartite.size() == 1) and (is_bipartite[0]))
        {
            Number cur_size = solve_bipartite_vertex_cover(bipartite_states, changes);

            undo_changes(changes);
            return cur_size;
        }

        if (all_cc.size() > 1)
        {
            Number cur_size = solve_on_connected_components(all_cc, changes);
            undo_changes(changes);
            return cur_size;
        }
    }
    
    //no connected component -> branch on vertex

    //preprocess branch_set, ensure that last vertex is undecided
    if (not branch_set.empty())
    {
        Vertex v = branch_set.back();
        while (not branch_set.empty() and _state[v] != undecided)
        {
            branch_set.pop_back();
            if (not branch_set.empty())
            {
                v = branch_set.back();
            }
        }
    }

    Vertex branch_vertex = _num_vertices;
    if (not branch_set.empty())
    {
        //branch_set not empty, last vertex is undecided
        branch_vertex = branch_set.back();
        branch_set.pop_back();
    }
    else
    {
        //select highest degree vertex for branching
        //go to first undecided vertex
        Vertex start_vertex = 0;
        for (; start_vertex < _num_vertices and _state[start_vertex] != undecided; ++start_vertex)
        {}

        if (start_vertex == _num_vertices)
        {
            Number solution_size = vertex_cover_size();
            if (solution_size < _best_size)
            {
                //update optimal solution
                _best_size = solution_size;
                StatusVector tmp_state = copy_vec(_state);
                update_best_vertex_cover(tmp_state);
            }
            undo_changes(changes);
            return solution_size;
        }

        Vertex deg3_vertex = _num_vertices;
        Vertex opt_vertex = start_vertex;
        double opt_degree = small_degree_neighbors(opt_vertex, 3);
        for (Vertex v = start_vertex + 1; v < _num_vertices; ++v)
        {
            if (_state[v] == undecided)
            {
                Number weighted_degree = small_degree_neighbors(v, 3);
                if (weighted_degree > opt_degree)
                {
                    opt_vertex = v;
                    opt_degree = weighted_degree;
                }
                if (weighted_degree < 4)
                {
                    deg3_vertex = v;
                }
            }
        }

        if (opt_degree == 5 and deg3_vertex != _num_vertices)
        {
            branch_vertex = deg3_vertex;     
        }
        else 
        {
            branch_vertex = opt_vertex;
        }
    }

    changes.push_back(branch_vertex);

    Number opt_size = _num_vertices;
    if (_degrees[branch_vertex] == 3)
    {
        opt_size = branch_on_degree_three(branch_vertex, lb);
    }
    else if (_degrees[branch_vertex] == 4)
    {

        opt_size = branch_on_small_neighborhood(branch_vertex, lb);
    }
    else
    {
        opt_size = branch_on_vertex(branch_vertex, changes, lb, branch_set);
    }

    undo_changes(changes);

    return opt_size;
}

Number Graph::branch_on_vertex(Vertex const v, ChangeVector & changes, Number lb, vector<Vertex> branch_set)
{
    //put v into vertex cover
    _state[v] = taken; //dont want change_state here
    update_neighbor_degrees(v);
    if (lb > 0)
    {
        --lb;
    }
   
    Number first_branch = branch(lb, branch_set);

    //do not put v into vertex cover
    _state[v] = discarded; //dont want change_state here
    take_neighbors(v, changes, lb); // no need to update degrees since they are all INFTY after undo changes

    Number second_branch = branch(lb, branch_set);

    return min(first_branch, second_branch);
}

Number Graph::branch_on_degree_three(Vertex v, Number lb)
{
    //we keep only the neighbors of the neighbors in changes

    ChangeVector changes;
    vector<Vertex> empty;
    //we know that v has degree 3
    vector<Vertex> neighbors;
    neighbors.reserve(3);
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        Vertex w = _adj[v][i];
        if (_state[w] == undecided)
        {
            neighbors.push_back(w);
        }
    }

    //negative branch
    _state[v] = discarded;
    update_neighbor_degrees(v);

    take_neighbors(v, changes, lb);
    for (Number i = 0; i < neighbors.size(); ++i)
    {
        update_neighbor_degrees(neighbors[i]);
    }
    Number best_branch = branch(lb, empty);


    undo_changes(changes);

    if (lb > 0)
    {
        --lb;
    }
    _state[v] = taken;
    bool valid_case = true;
    for (Number i = 0; i < 4; ++i)
    {
        valid_case = true;
        for (Number j = 0; j < neighbors.size(); ++j)
        {
            _state[neighbors[j]] = undecided;
        }
        for (Number j = 0; j < 3; ++j)
        {
            if (j == i)
            {
                _state[neighbors[j]] = taken;
            }
            else
            {
                if (_state[neighbors[j]] == undecided)
                {
                    _state[neighbors[j]] = discarded;
                    take_neighbors(neighbors[j], changes, lb);
                }
                else
                {
                    valid_case = false;
                    break;
                }
            }
        }
        if (valid_case)
        {

            best_branch = min(best_branch, branch(lb, empty));

        }
        undo_changes(changes);
    }

    _state[v] = undecided;
    for (Number i = 0; i < neighbors.size(); ++i)
    {
        _state[neighbors[i]] = undecided;
    }
    return best_branch;
}

Number Graph::branch_on_small_neighborhood(Vertex v, Number lb)
{
    //we keep only the neighbors of the neighbors in changes

    ChangeVector changes;
    vector<Vertex> empty;

    vector<Vertex> neighbors;
    neighbors.reserve(_degrees[v]);
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        Vertex w = _adj[v][i];
        if (_state[w] == undecided)
        {
            neighbors.push_back(w);
        }
    }

    //negative branch
    _state[v] = discarded;
    if (lb > 0)
    {
        --lb;
    }
    update_neighbor_degrees(v);

    take_neighbors(v, changes, lb);
    for (Number i = 0; i < neighbors.size(); ++i)
    {
        update_neighbor_degrees(neighbors[i]);
    }
    Number best_branch = branch(lb, empty);
    undo_changes(changes);

    //positive branch, at least 2 neighbors of v needs to be discarded
    _state[v] = taken;
    if (lb > 0)
    {
        lb--;
    }


    bool valid_case = true;
    for (Number i = 0; i < neighbors.size(); ++i)
    {
        for (Number j = i + 1; j < neighbors.size(); ++j)
        {
            valid_case = true;
            //in this iteration neighbor[i] and neighbor[j] should be discarded and all other neighbors should be undecided
            for (Number k = 0; k < neighbors.size(); ++k)
            {
                _state[neighbors[k]] = undecided;
            }
            _state[neighbors[i]] = discarded;
            take_neighbors(neighbors[i], changes, lb);
            if (_state[neighbors[j]] == undecided) // test if neighbor[i] and neighbor[j] are adjacent
            {
                _state[neighbors[j]] = discarded;
                take_neighbors(neighbors[j], changes, lb);
            }
            else
            {
                valid_case = false;
            }

            if (valid_case)
            {
                best_branch = min(best_branch, branch(lb, empty));
            }
            undo_changes(changes);
        }
    }

    _state[v] = undecided;
    for (Number i = 0; i < neighbors.size(); ++i)
    {
        _state[neighbors[i]] = undecided;
    }

    return best_branch;
}


Number const Graph::vertex_cover_size()
{
    Number vc_size = 0;
    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        if (_state[v] == taken)
        {
            ++vc_size;
        }
    }
    return vc_size + (_num_vertices - _num_orig_vertices);
}

void Graph::update_best_vertex_cover(StatusVector & tmp_state)
{
    //for undoing folds we make some computations on the state vector
    //since we do not want to incur any changes on the graph, we do this on a copy of the states
    _best_vc.clear();
    for (Vertex v = 0; v < _num_orig_vertices; ++v)
    {
        if (tmp_state[v] == taken)
        {
            _best_vc.push_back(v);
        }
    }
    
    //important that we iterate in reverse for handling successive folds
    for (long i = _folds.size() - 1; i >= 0; --i)
    {
        FoldData & fold = _folds[i];
        Vertex artificial = fold._new_vertex;

        if (tmp_state[artificial] == discarded)
        {
            if (fold._when_discarded < _num_orig_vertices)
            {
                _best_vc.push_back(fold._when_discarded);
            }
            //must update state to account for successive folds
            tmp_state[fold._when_discarded] = taken;

            for (Number j = 0; j < fold._when_taken.size(); ++j)
            {        
                tmp_state[fold._when_taken[j]] = discarded;
            }
        }
        else if (tmp_state[artificial] == taken)
        {
            for (Number j = 0; j < fold._when_taken.size(); ++j)
            {
                if (fold._when_taken[j] < _num_orig_vertices)
                {
                    _best_vc.push_back(fold._when_taken[j]);
                }
                tmp_state[fold._when_taken[j]] = taken;
            }
            
            tmp_state[fold._when_discarded] = discarded;
        }
    }
}

bool const Graph::validate_vertex_cover()
{
    VertexCover const & vc = _best_vc;
    vector<bool> taken(_num_vertices, false);
    for (Number i = 0; i < vc.size(); ++i)
    {
        taken[vc[i]] = true;
    }

    bool result = true;
    //loop through all edges and check whether
    //at least one of the endpoints is taken
    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        for (Number i = 0; i < _adj[v].size(); ++i)
        {
            Vertex w = _adj[v][i];
            if (not taken[v] and not taken[w] and v < _num_orig_vertices and w < _num_orig_vertices)
            {
                cerr << "{ " << v << ", " << w << "} is not covered" << endl;
                result = false;
            }
        }
    }
    return result;
}

void Graph::apply_reduction_rules_and_lower_bounds(ChangeVector & changes, Number & lb)
{
    Number lower_bound = 0;
    Number undecided_vertices = 0;

    if (not _skip_degree_reduction_once)
    {
        vector<Vertex> process_stack;
        process_stack.reserve(_num_vertices);

        //initialize processing stack with all undecided vertices
        for (Vertex v = 0; v < _num_vertices; ++v)
        {
            if (_state[v] == undecided)
            {
                process_stack.push_back(v);
                undecided_vertices++;
            }
        }

        while (not process_stack.empty() )
        {
            Vertex v = process_stack.back();
            process_stack.pop_back();

            if (_state[v] == undecided)
            {
                if (_degrees[v] == INFTY)
                {
                    _degrees[v] = degree(v);
                }

                switch (_degrees[v])
                {
                    case 0:
                        apply_red_degree_zero(changes, v);
                        if (lb > 0)
                        {
                            lb--;
                        }
                        undecided_vertices--;
                        break;
                    case 1:
                        apply_red_degree_one(changes, v, process_stack);
                        if (lb > 0)
                        {
                            lb--;
                        }
                        undecided_vertices = undecided_vertices - 2;
                        break;
                    case 2:
                        apply_red_degree_two(changes, v, process_stack, lb);
                        undecided_vertices = undecided_vertices - 2;
                        break;
                }
            }

        }
    }
    else
    {
        recompute_degrees();
        while (apply_dominance_rule(changes, lb))
        {
            //do nothing
        }
        _skip_degree_reduction_once = false;
        apply_reduction_rules_and_lower_bounds(changes, lb);
    }

    //recompute the lowerbound if the fraction is less than 25%
    if ((double) lb / (double) undecided_vertices < 0.25)
    {
        lb = max(lower_bound, greedy_maximal_matching());
    }
}


void Graph::apply_red_degree_zero(ChangeVector & changes, Vertex v)
{
    change_state(v, discarded, changes);
}

void Graph::apply_red_degree_one(ChangeVector & changes, Vertex v, vector<Vertex> & process_stack)
{
    change_state(v, discarded, changes);
    //search for the neighbor of v
    Vertex u = -1;
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        if (_state[_adj[v][i]] == undecided)
        {
            u = _adj[v][i];
            break;
        }
    }

    change_state(u, taken, changes);

    //update degrees of neighbor of u
    for (Number i = 0; i < _adj[u].size(); ++i)
    {
        if (_state[_adj[u][i]] == undecided)
        {
            Vertex w = _adj[u][i];

            if (_degrees[w] == INFTY)
            {
                _degrees[w] = degree(w);
            }
            else
            {
                _degrees[w] = _degrees[w] - 1;
            }

            if (_degrees[w] < 3)
            {
                process_stack.push_back(w);
            }
        }
    }
}

void Graph::apply_red_degree_two(ChangeVector & changes, Vertex v, vector<Vertex> & process_stack, Number & lb)
{
    vector<Vertex> neighbors_of_v;

    add_undecided_neighbors_to_vector(v, neighbors_of_v);

    Vertex u = neighbors_of_v[0];
    Vertex w = neighbors_of_v[1];

    // check, if one neighbor has degree 1
    if (_degrees[u] == INFTY)
    {
        _degrees[u] = degree(u);
    }

    if (_degrees[w] == INFTY)
    {
        _degrees[w] = degree(w);
    }

    if (_degrees[u] == 1)
    {
        apply_red_degree_one(changes, u, process_stack);
        return;
    }
    if (_degrees[w] == 1)
    {
        apply_red_degree_one(changes, w, process_stack);
        return;
    }

    //check, if there is the edge {u,w} in the graph:
    if (edge_exists(u, w))
    {
        //There is the edge {u,w}. Thus, Take u and w
        change_state(v, discarded, changes);
        change_state(u, taken, changes);
        change_state(w, taken, changes);

        // update degree of neighbors of u and w
        for (Number i = 0; i < _adj[u].size(); ++i)
        {
            if (_state[_adj[u][i]] == undecided)
            {
                Vertex x = _adj[u][i];
                if (_degrees[x] == INFTY)
                {

                }//you must not compute the degree here, since it could be decreased in next for loop (which would be wrong)
                else
                {
                    _degrees[x] = _degrees[x] - 1;
                }

                if (_degrees[x] < 3)
                {
                    process_stack.push_back(x);
                }

            }
        }

        for (Number i = 0; i < _adj[w].size(); ++i)
        {
            if (_state[_adj[w][i]] == undecided)
            {
                Vertex x = _adj[w][i];

                if (_degrees[x] == INFTY)
                {
                    _degrees[x] = degree(x);
                }
                else
                {
                    _degrees[x] = _degrees[x] - 1;
                }

                if (_degrees[x] < 3)
                {
                    process_stack.push_back(x);
                }

            }
        }

        if (lb > 2)
        {
            lb = lb - 2;
        }
        else
        {
            lb = 0;
        }
    }
    else
    {
        //edge {u,w} does not exist, so fold
        change_state(v, folded, changes);
        vector<Vertex> neighbors;
        add_undecided_neighbors_to_vector(u, neighbors);
        add_undecided_neighbors_to_vector(w, neighbors);

        remove_duplicates_with_degree_update(neighbors, process_stack);

        FoldData fold_data(_num_vertices, v, neighbors_of_v);
        _folds.push_back(fold_data);

        change_state(u, folded, changes);
        change_state(w, folded, changes);

        add_vertex(neighbors, changes);

        //add the new vertex to the degrees
        _degrees.push_back(neighbors.size());
        //add it to the process_stack
        process_stack.push_back(_num_vertices - 1);
        if (lb > 0)
        {
            lb--;
        }
    }
}

bool is_subset_of(vector<Vertex> const & subset, vector<Vertex> const & superset)
{
    Number j = 0;
    for (Number i = 0; i < subset.size(); ++i)
    {
        while (j < superset.size() and superset[j] != subset[i])
        {
            ++j;
        }
        if (j == superset.size())
        {
            return false;
        }
    }
    return true;
}

bool Graph::apply_dominance_rule(ChangeVector & changes, Number & lb)
{
    //degrees are computed beforehand
    vector<Vertex> undecided_vertices;
    undecided_vertices.reserve(_num_vertices);
    vector<Vertex> undecided_vertices_inverse(_num_vertices, _num_vertices);

    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        if (_state[v] == undecided)
        {
            undecided_vertices_inverse[v] = undecided_vertices.size();
            undecided_vertices.push_back(v);
        }
    }

    //construct sorted adjacency lists for undecided vertices (closed neighborhood!)
    vector<AdjacencyList> sorted_adj_list(undecided_vertices.size());
    for (Number i = 0; i < undecided_vertices.size(); ++i)
    {
        Vertex v = undecided_vertices[i];
        sorted_adj_list[i].push_back(v);
        for (Number j = 0; j < _adj[v].size(); ++j)
        {
            Vertex w = _adj[v][j];
            if (_state[w] == undecided)
            {
                sorted_adj_list[i].push_back(w);
            }
        }
        sort(sorted_adj_list[i].begin(), sorted_adj_list[i].end());
    }

    vector<bool> take(undecided_vertices.size(), false);

    //check for dominance
    for (Number i = 0; i < undecided_vertices.size(); ++i)
    {
        if (not take[i])
        {
            Vertex v = undecided_vertices[i];

            //go to first neighbor with a larger index than v
            Number j = 0;
            for (; sorted_adj_list[i][j] < v; ++j)
            {
                //do nothing
            }
            ++j;

            for (; j < sorted_adj_list[i].size(); ++j)
            {
                Vertex w = sorted_adj_list[i][j];
                Number k = undecided_vertices_inverse[w];

                if (_degrees[w] < _degrees[v] and not take[i] and is_subset_of(sorted_adj_list[k], sorted_adj_list[i]))
                {
                    take[i] = true;
                }
                else if (_degrees[w] >= _degrees[v] and is_subset_of(sorted_adj_list[i], sorted_adj_list[k]))
                {
                    take[k] = true;
                }
            }
        }
    }

    Number count = 0;
    bool changed = false;
    for (Number i = 0; i < undecided_vertices.size(); ++i)
    {
        if (take[i])
        {
            ++count;
            changed = true;
            Vertex v = undecided_vertices[i];
            change_state(v, taken, changes);

            //update degrees
            for (Number j = 0; j < sorted_adj_list[i].size(); ++j)
            {
                Vertex w = sorted_adj_list[i][j];
                if (v != w)
                {
                    --_degrees[w];
                }
            }

            //update lb
            if (lb > 0)
            {
                --lb;
            }
        }
    }
    return changed;
}

void Graph::undo_changes(ChangeVector & changes)
{
    //important that we iterate in reverse for handling successive folds
    for (long i = changes.size() - 1; i >= 0; --i)
    {
        Change current = changes[i];
        switch (current._type)
        {
            case state_change:
                _state[current._vertex] = undecided;
                break;
            case new_vertex_due_to_fold:
                delete_last_vertex(); //also decreases _num_vertices
                _folds.pop_back();
                break;
        }

    }
    changes.clear();
    _degrees.clear();
    _degrees.assign(_num_vertices, INFTY);
}

void Graph::take_neighbors(Vertex v, ChangeVector & changes, Number & lb)
{
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        Vertex neighbor = _adj[v][i];
        if (_state[neighbor] == undecided)
        {
            change_state(neighbor, taken, changes);
            if (lb > 0)
            {
                --lb;
            }
        }
    }
}

map<Vertex, Vertex> Graph::compute_reassignments_from_graph_to_cc(ConnectedComponent const & cc)
{
    map<Vertex, Vertex> reassignments;
    Vertex next = 0;
    for (Number i = 0; i < cc.size(); ++i)
    {
        reassignments[cc[i]] = next;
        ++next;
    }
    return reassignments;
}

map<Vertex, Vertex> Graph::compute_reassignments_from_cc_to_graph(ConnectedComponent const & cc)
{
    map<Vertex, Vertex> reassignments;
    Vertex next = 0;
    for (Number i = 0; i < cc.size(); ++i)
    {
        reassignments[next] = cc[i];
        ++next;
    }
    return reassignments;
}

void Graph::copy_connected_component(Graph & cc_graph,
                                     ConnectedComponent const & cc,
                                     map<Vertex, Vertex> & reassignments)
{
    cc_graph._num_vertices = cc.size();
    cc_graph._num_orig_vertices = cc.size();
    cc_graph._adj.resize(cc.size());
    cc_graph._num_edges = 0;
    cc_graph._skip_degree_reduction_once = true;
    cc_graph._folds.clear();

    for (Number i = 0; i < cc.size() ; ++i)
    {
        cc_graph._adj[i].resize(0);
    }

    Vertex neighbor;
    for (Number i = 0; i < cc.size(); ++i)
    {
        for (Number j = 0; j < _adj[cc[i]].size(); ++j)
        {
            neighbor = _adj[cc[i]][j];
            if (_state[neighbor] == undecided)
            {
                cc_graph._adj[i].push_back(reassignments[neighbor]);
                cc_graph._num_edges += 1;
            }
        }
    }
    cc_graph._num_edges = cc_graph._num_edges / 2;

    //initialize vertex states
    cc_graph._state.resize(cc_graph._num_vertices, undecided);

    cc_graph._best_size = cc_graph._num_vertices;

}

Number Graph::solve_on_connected_components(vector<ConnectedComponent> & all_cc, ChangeVector & changes)
{
    Number total_vc_size = 0; // count the number of vertices we need, to cover all cc
    VertexCover vc_on_all_cc; // Store the best vertex cover of the cc

    Graph cc_graph;
    map<Vertex, Vertex> reassignments_g_to_cc;
    map<Vertex, Vertex> reassignments_cc_to_g;

    //solve each connected component separately
    for (Number current_cc_number = 0; current_cc_number < all_cc.size() ; ++current_cc_number)
    {
        reassignments_g_to_cc = compute_reassignments_from_graph_to_cc(all_cc[current_cc_number]);
        reassignments_cc_to_g = compute_reassignments_from_cc_to_graph(all_cc[current_cc_number]);

        copy_connected_component(cc_graph, all_cc[current_cc_number], reassignments_g_to_cc);

        cc_graph._degrees.assign(cc_graph._num_vertices, INFTY);

        for (Vertex v = 0; v < cc_graph._num_vertices ; ++v)
        {
            cc_graph._degrees[v] = cc_graph.degree(v);
        }

        vector<Vertex> empty;

        Number part_solution = cc_graph.branch(0, empty);
        total_vc_size += part_solution;

        //Store the best vertex cover of this connected component
        for (Number i = 0 ; i < cc_graph._best_vc.size(); ++i)
        {
            vc_on_all_cc.push_back(reassignments_cc_to_g[cc_graph._best_vc[i]]);
        }
    }

    Number vertex_cover_on_this_branch = vertex_cover_size(); // to get all taken vertices up to the connected components
    vertex_cover_on_this_branch += total_vc_size; // add the vertex cover of the connected components

    if (vertex_cover_on_this_branch < _best_size)
    {
        _best_size = vertex_cover_on_this_branch;

        //modify state vector accordingly to transfer the solutions on the CCs to the original graph
        StatusVector tmp_state = copy_vec(_state);

        for (Number i = 0; i < vc_on_all_cc.size() ; ++i)
        {
            tmp_state[vc_on_all_cc[i]] = taken;
        }
        for (Vertex v = _num_orig_vertices; v < _num_vertices; ++v)
        {
            if (tmp_state[v] != taken)
            {
                tmp_state[v] = discarded;
            }
        }
        update_best_vertex_cover(tmp_state);
    }

    return vertex_cover_on_this_branch;
}


Number Graph::greedy_maximal_matching()
{
    vector<bool> used(_num_vertices, false);

    Number count = 0;

    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        if (_state[v] == undecided and not used[v])
        {
            for (Number i = 0; i < _adj[v].size(); ++i)
            {
                Vertex w = _adj[v][i];
                if (_state[w] == undecided and not used[w])
                {
                    used[v] = true;
                    used[w] = true;
                    ++count;
                    break;
                }
            }
        }
    }

    return count;
}

void Graph::alternating_distances(vector<Vertex> & dfs_startpoints,
                                  vector<Vertex> const & left_side,
                                  vector<MaybeVertex> const & matching_partner,
                                  vector<Number> & distance)
{
   
    distance.assign(_num_vertices, INFTY);

    vector<Vertex> cur_layer;

    //put all exposed vertices in the first copy into the cur_layer
    for (Number i = 0; i < left_side.size(); ++i)
    {
        Vertex v = left_side[i];
        if (matching_partner[v] == -1)
        {
            cur_layer.push_back(v);
            distance[v] = 0;
        }
    }

    BipState cur_side = BipState::left_side;
    Number cur_distance = 0;
    vector<Vertex> next_layer;
    do
    {
        while (not cur_layer.empty())
        {
            Vertex v = cur_layer.back();
            cur_layer.pop_back();
            if (cur_side == BipState::left_side)
            {
                for (Number i = 0; i < _adj[v].size(); ++i)
                {
                    Vertex w = _adj[v][i];
                    if (_state[w] == undecided and distance[w] == INFTY)
                    {
                        distance[w] = cur_distance + 1;
                        next_layer.push_back(w);
                    }
                }
            }
            else
            {
                if (matching_partner[v] != -1)
                {
                    //w cannot have an already assigned distance
                    //because there are no edges inside the left side and w is in the left side
                    //and there is only one incident matching edge (i.e. from the right side)
                    //which we have just found, namely {v,w}
                    Vertex w = static_cast<Vertex>(matching_partner[v]);
                    distance[w] = cur_distance + 1;
                    next_layer.push_back(w);
                }
                else
                {
                    dfs_startpoints.push_back(v);
                }
            }
        }
        ++cur_distance;
        cur_side = cur_side == BipState::left_side ? BipState::right_side : BipState::left_side;
        cur_layer = next_layer;
        next_layer.clear();
    } while (not cur_layer.empty());
}

void Graph::backwards_dfs(vector<Vertex> const & dfs_startpoints,
                          vector<MaybeVertex> & matching_partner,
                          vector<Number> const & distance)
{
    //vertices that are not undecided should have infinite distance
    vector<Vertex> stack;
    vector<bool> visited(_num_vertices, false);
    vector<Vertex> predecessor(_num_vertices, INFTY);

    for (Number j = 0; j < dfs_startpoints.size(); ++j)
    {
        stack.push_back(dfs_startpoints[j]);
        while (not stack.empty())
        {
            Vertex v = stack.back();
            stack.pop_back();
            visited[v] = true;
            if (distance[v] == 0)
            {
                augment(v, predecessor, matching_partner);
                break;
            }
            else
            {
                if (distance[v] % 2 == 1) //right side
                {
                    for (Number i = 0; i < _adj[v].size(); ++i)
                    {
                        Vertex w = _adj[v][i];
                        if (not visited[w] and (distance[w] == distance[v] - 1))
                        {
                            stack.push_back(w);
                            predecessor[w] = v;
                        }
                    }
                }
                else //left side and not exposed
                {
                    Vertex w = static_cast<Vertex>(matching_partner[v]);
                    stack.push_back(w);
                    predecessor[w] = v;
                }
            }
        }
        stack.clear();
    }
}

void Graph::augment(Vertex v, vector<Vertex> & predecessor, vector<MaybeVertex> & matching_partner)
{
    int distance = 0;

    while (predecessor[v] != INFTY)
    {
        if (distance % 2 == 0)
        {
            matching_partner[v] = predecessor[v];
            matching_partner[predecessor[v]] = v;
        }
        else
        {
            //do nothing
        }
        v = predecessor[v];
        ++distance;
    }
}

void Graph::matching_to_vertex_cover(vector<Vertex> & vc,
                                     vector<BipState> const & bipartite_states,
                                     vector<Number> const & distance)
{
    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        if (bipartite_states[v] == BipState::left_side)
        {
            if (distance[v] == INFTY)
            {
                vc.push_back(v);
            }
        }
        else if (bipartite_states[v] == BipState::right_side)
        {
            if (distance[v] != INFTY)
            {
                vc.push_back(v);
            }
        }
    }
}

Number Graph::solve_bipartite_vertex_cover(vector<BipState> const & bipartite_states, ChangeVector & changes)
{
    VertexCover vc;
    vector<MaybeVertex> matching_partner(_num_vertices, -1);
    vector<Number> distance;

    vector<Vertex> left_side;
    left_side.reserve(_num_vertices);
    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        if (bipartite_states[v] == BipState::left_side)
        {
            left_side.push_back(v);
        }
    }

    vector<Vertex> dfs_startpoints;
    bool done = false;
    int counter = 0;
    while (not done)
    {
        dfs_startpoints.clear();
        alternating_distances(dfs_startpoints, left_side, matching_partner, distance);
        done = dfs_startpoints.empty();

        backwards_dfs(dfs_startpoints, matching_partner, distance);

        ++counter;
    }

    matching_to_vertex_cover(vc, bipartite_states, distance);

    Number cur_size = vertex_cover_size() + vc.size();

    if (cur_size < _best_size)
    {
        _best_size = cur_size;
        for (Number i = 0; i < vc.size(); ++i)
        {
            change_state(vc[i], taken, changes);
        }

        //set all undecided vertices to discarded
        for (Vertex v = 0; v < _num_vertices; ++v)
        {
            if (_state[v] == undecided)
            {
                change_state(v, discarded, changes);
            }
        }

        StatusVector tmp_state = copy_vec(_state);
        update_best_vertex_cover(tmp_state);
    }

    return cur_size;
}

void Graph::update_neighbor_degrees(Vertex v)
{
    for (Number i = 0; i < _adj[v].size(); ++i)
    {
        if (_state[_adj[v][i]] == undecided)
        {
            Vertex w = _adj[v][i];

            if (_degrees[w] == INFTY)
            {
                _degrees[w] = degree(w);
            }
            else
            {
                _degrees[w] = _degrees[w] - 1;
            }
        }
    }
}

void Graph::compute_modulator_to_bipartite(vector<Vertex> & branch_set)
{
    vector<BipState> bipartite_states(_num_vertices, BipState::non_existent);

    //set vertices in G[undecided] to unvisited
    for (Vertex v = 0; v < _num_vertices; ++v)
    {
        if (_state[v] == Status::undecided)
        {
            bipartite_states[v] = BipState::unvisited;
        }
    }

    vector<bool> is_bipartite;

    vector<ConnectedComponent> all_cc = compute_connected_components_BFS(bipartite_states, is_bipartite);

    Number count_undecided = 0;
    Number count_destructive = 0;

    for (Number i = 0 ; i < bipartite_states.size(); ++i)
    {
        if (bipartite_states[i] != non_existent)
        {
            count_undecided++;
        }
        if (bipartite_states[i] == destructive)
        {
            count_destructive++;
        }
    }

    if (count_destructive < 50)
    {
        branch_set.clear();
        for (Vertex v = 0; v < _num_vertices; ++v)
        {
            if (bipartite_states[v] == destructive)
            {
                branch_set.push_back(v);
            }
        }
        sort(branch_set.begin(), branch_set.end(), [this](Vertex v, Vertex w)
        {
            return _degrees[v] < _degrees[w];
        });
    }
}
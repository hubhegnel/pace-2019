#ifndef _graph_h_
#define _graph_h_

#include <vector>
#include <map>
typedef unsigned Vertex;
typedef unsigned Degree;
typedef unsigned Number;
typedef int MaybeVertex;
typedef std::vector<Vertex> AdjacencyList;
typedef std::vector<Vertex> VertexCover;
typedef std::vector<Vertex> ConnectedComponent;
typedef std::vector<Number> DistanceVector;



/*
 *  Removes duplicate entries in vector vertices, retains the order.
 */
void remove_duplicates(std::vector<Vertex> & vertices);

/*
 *  Describes the current status of a vertex during the branching.
 */
enum Status
{
    undecided,    //not yet decided
    discarded,    //not in VC
    taken,        //in VC
    folded        //removed by degree-2 reduction
};

typedef std::vector<Status> StatusVector;

/*
 *  The possible states when determining bipartiteness and a modulator.
 */
enum BipState
{
    non_existent,     //already taken
    unvisited,        //not yet explored
    left_side,        //left side of bipartition
    right_side,       //right side of bipartition
    destructive       //modulator to bipartite set (greedy)
};

/*
 *  Describes how a vertex has changed.
 */
enum ChangeType
{
    state_change,
    new_vertex_due_to_fold
};

/*
 *  Describes how and which vertex has changed.
 *  Used for dynamically changing the graph when moving between different branches,
 *  to avoid creating unneccessary copies. I.e., we are essentially keeping a stack
 *  of changes that are reverted after returning from some branch.
 */
struct Change
{
    ChangeType _type;
    Vertex _vertex;

    Change(Vertex v);
    Change(ChangeType type, Vertex v);
};

typedef std::vector<Change> ChangeVector;

/*
 *  Bookkeeping structure for folds, i.e. degree-2 reductions.
 *  Stores the ID of the newly created vertex and which vertices to
 *  choose for the vertex cover based on the state of the newly created vertex.
 */
struct FoldData
{
    //the vertex created by this fold
    Vertex _new_vertex;

    //the original vertex that must be chosen when we discard new_vertex
    Vertex _when_discarded;

    //the original vertices that must be chosen when we take new_vertex
    std::vector<Vertex> _when_taken;

    //create FoldData with the corresponding entries
    FoldData(Vertex new_vertex, Vertex when_discarded, std::vector<Vertex> & when_taken);
};

/*
 *  The graph on which the solver works. At every step we consider the graph induced by the undecided vertices.
 *  The graph is updated dynamically by the algorithm, only when we split into connected components will we
 *  create copies of of parts of the graph.
 */
struct Graph
{
    //Number of vertices, including the vertices that were created by degree 2 reduction
    Number _num_vertices;
    Number _num_edges;

    //The adjacency lists of all vertices
    std::vector<AdjacencyList> _adj;

    //the states of all vertices
    StatusVector _state;

    //best vertex cover seen so far
    VertexCover _best_vc;

    //size of best vertex cover seen so far
    Number _best_size;

    //whether we skip the next degree reduction, used after graph splits into several connected components
    //to avoid redundant computations
    bool _skip_degree_reduction_once;

    //Number of vertices that were not created by degree 2 reduction.
    Number _num_orig_vertices;

    //handling folds
    std::vector<FoldData> _folds;

    std::vector<Degree> _degrees;

    /*
     *  Adds a new vertex to the graph with ID _num_vertices, increases counters
     *  and vectors appropriately. The new vertex will be connected to the vertices
     *  in neighbors and no others. An appropriate change will be pushed onto changes.
     */
    void add_vertex(std::vector<Vertex> & neighbors, ChangeVector & changes);

    /*
     *  Removes duplicates in vertices, order possibly changes, and updates _degrees[v] for all v in vertices.
     *  Vertices with degree at most 2 are pushed onto the process stack.
     */
    void remove_duplicates_with_degree_update(std::vector<Vertex> & vertices, std::vector<Vertex>  & process_stack);

    /*
     *  Deletes the last vertex in the graph, used for undoing an add vertex operation.
     */
    void delete_last_vertex();

    /*
     *  Change the state of vertex v to new_state and keep track of it in the ChangeVector changes.
     */
    void change_state(Vertex v, Status new_state, ChangeVector & changes);

    /*
     *  Test whether the edge between v and w exists in the graph.
     */
    bool edge_exists(Vertex v, Vertex w);

    /*
     *  Push the neighbors of v that are undecided onto the vector vertices.
     */
    void add_undecided_neighbors_to_vector(Vertex v, std::vector<Vertex> & vertices);

    /*
    *	@return vector of all the connected components
    */
    std::vector<ConnectedComponent> compute_connected_components_BFS(std::vector<BipState> & bipartite_states,
                                                                     std::vector<bool> & is_bipartite);

    /*
     *  Perform a breadth-first-search on the graph to determine whether it's bipartite (if not compute a modulator)
     *  and determines the connected components
     */
    bool BFS_and_bipartite_test(Vertex v,
                                std::vector<BipState> & bipartite_states,
                                std::vector<ConnectedComponent> & cc,
                                Number i);

   /*
    *	Computes the degree of the vertex v, taking into account the current state.
    */
    Degree degree(Vertex v);

    /*
     *  Recompute all entries of _degrees.
     */
    void recompute_degrees();

    /*
     *  Computes a "weighted" degree of vertex v, where all neighbors with a degree smaller than small_degree
     *  count as 1.25 vertices.
     */
    double small_degree_neighbors(Vertex v, Degree small_degree);

    /*
    *	Branch on vertex with highest degree.
    *	Returns size of best solution found in this subtree.
    */
    Number branch(Number lb, std::vector<Vertex> branch_set);

    /*
     *  Branch on vertex v.
     *  Returns size of best solution found in this subtree.
     */
    Number branch_on_vertex(Vertex const v, ChangeVector & changes, Number lb, std::vector<Vertex> branch_set);

    /*
     *  Specialized branching for vertices with degree 3. When taking v, we ensure that at least 2 of its neighbors
     *  are discarded, the remaining neighbors will be taken.
     */
    Number branch_on_degree_three(Vertex v, Number lb);

    /*
     *  Similar, but subtly different, to branch_on_degree_three for degrees larger than 3. When taking v, we discard 2
     *  of its neighbors and the remaining neighbors stay undecided.
     */
    Number branch_on_small_neighborhood(Vertex v, Number lb);

    /*
    *	Computes the size of the (partial) vertex cover
    *	as currently indicated by the vertex states.
    */
    Number const vertex_cover_size();

    /*
    *	Stores the current vertex cover as specified by the state vector
    *	into _best_vc. For fold handling some changes are performed on the
     *	state vector, hence we must pass changes as argument.
    */
    void update_best_vertex_cover(StatusVector & tmp_state);

    /*
    *	Validates whether the vertex cover stored in _best_vc
    *	actually covers all edges.
    *	Does NOT check optimality.
    */
    bool const validate_vertex_cover();

   
    /*
    *	Applies reduction rules on the current graph.
    *	The graph modifications inflicted by the reduction rules
    *	are appended to changes.
    */
    void apply_reduction_rules_and_lower_bounds(ChangeVector & changes, Number & lb);

    /*
    * 	Applies reduction rule for degree 0 on the vertex v
    *
    */
    void apply_red_degree_zero(ChangeVector & changes, Vertex v);

    /*
    * 	Applies reduction rule for degree 1 on vertex v
    *   If adjacent vertices of v are of degree < 3 after the reduction, those vertices are pushed to the process_stack
    */
    void apply_red_degree_one(ChangeVector & changes, Vertex v, std::vector<Vertex>  & process_stack);

    /*
    * 	Applies reduction rule for degree 2 vertex v 
    *   If adjacent vertices of v are of degree < 3 after the reduction, those vertices are pushed to the process_stack
    */
    void apply_red_degree_two(ChangeVector & changes, Vertex v, std::vector<Vertex> & process_stack, Number & lb);

    /*
    *   Applies the dominance reduction rule on the graph
    *   Return true if the reduction rule could be used
    */
    bool apply_dominance_rule(ChangeVector & changes, Number & lb);

   
    /*
    *	Reverts the graph modifications stored in changes.
    *	changes will be cleared afterwards.
    */
    void undo_changes(ChangeVector & changes);

    /*
    *	To be used upon discarding the vertex v.
    *	All undecided neighbors of v are put into
    *	the current vertex cover.
    */
    void take_neighbors(Vertex v, ChangeVector & changes, Number & lb);

 
    /*
    * creates a map from the vertices of the connected component (with arbitrary numbers attachted to) to the numbers
    * 0, 1, 2, ..., |cc|-1
    */
    std::map<Vertex, Vertex> compute_reassignments_from_graph_to_cc(ConnectedComponent const & cc);

/*
    * creates a map from the vertices of the connected component (with numbers 0 to |cc|-1 ) to the old, original numbers
    * of the graph
    */
    std::map<Vertex, Vertex> compute_reassignments_from_cc_to_graph(ConnectedComponent const & cc);

    /*
    * Initialize cc_graph, which is a induced subgraph of graph, induced by the vertices in cc
    * In reassignment (size=_num_vertices) one stores the position of the vertex in its connected component
    */
    void copy_connected_component(Graph & cc_graph, ConnectedComponent const & cc, std::map<Vertex, Vertex> & reassignments);

    /*
     *  Compute optimal Vertex Covers on the Connected Components in all_cc and combine these with the
     *  current partial Vertex Cover to obtain a feasible solution. Updates best solution if we beat
     *  the previous one. Size of computed solution is returned. Changes must be passed to this function
     *  as we need to undo the changes from reduction rules before returning from this branch.
     */
    Number solve_on_connected_components(std::vector<ConnectedComponent> & all_cc, ChangeVector & changes);

    /*
     *  Solve vertex cover on bipartite graphs in polynomial time by reduction to maximum matching.
     *  The bipartition is given by bipartite_states. To solve maximum matching we use a slightly modified
     *  Hopcroft-Karp algorithm.
     */
    Number solve_bipartite_vertex_cover(std::vector<BipState> const & bipartite_states, ChangeVector & changes);

    /*
     *  First part of Hopcroft-Karp: Perform an alternating BFS starting at the exposed vertices on the left side
     *  of the bipartition. The current matching is given by matching_partner, where the endpoints of a matching
     *  edge point to each other. All exposed vertices have matching_partner = -1. The distances computed by the
     *  alternating BFS are stored in distance. The reachable exposed vertices on the right side of the bipartition
     *  are stored in dfs_startpoints.
     */
    void alternating_distances(std::vector<Vertex> & dfs_startpoints,
                               std::vector<Vertex> const & left_side,
                               std::vector<MaybeVertex> const & matching_partner,
                               std::vector<Number> & distance);

    /*
     *  Second part of Hopcroft-Karp: Perform a backwards DFS starting at dfs_startpoints. Every edge we take
     *  has to decrease the distance by exactly 1. Upon reaching an exposed vertex on the left side of the
     *  bipartition, we have found an augmenting path, which we use to enlargen the matching.
     */
    void backwards_dfs(std::vector<Vertex> const & dfs_startpoints,
                       std::vector<MaybeVertex> & matching_partner,
                       std::vector<Number> const & distance);

   /*
    *  Augmentation in Hopcroft-Karp algorithm.
    *  Taking the predecessors, found by backwards DFS, starting at v, we flip the matching along the resulting path.
    */
    void augment(Vertex v, std::vector<Vertex> & predecessor, std::vector<MaybeVertex> & matching_partner);

   /*
    *  Via Königs Theorem we transform the maximum matching to a minimum vertex cover. We dont need the explicit
    *  matching, the distances computed by the final alternating BFS suffice to get the vertex cover.
    */
    void matching_to_vertex_cover(std::vector<Vertex> & vc,
                                  std::vector<BipState> const & bipartite_states,
                                  std::vector<Number> const & distance);

   /*
    *  Greedily compute a maximal matching of the graph, which we use as a lower bound for the vertex cover number.
    */
    Number greedy_maximal_matching();

    /*
     *  Update the _degrees of the neighbors of v, when v was taken or discarded.
     */
    void update_neighbor_degrees(Vertex v);

    /*
     *  Greedily compute an odd cycle transversal, i.e. modulator to bipartite, and if it's small enough we prioritise
     *  branching on this modulator.
     */
    void compute_modulator_to_bipartite(std::vector<Vertex> & branch_set);

};

#endif
#pragma once

#include <stdint.h>

#include <vector>

/*
 * This is a module for manipulating general undirected graphs (vertices and
 * edges). Note: In this representation, a node is connected by other node with
 * exactly one edge.
 *
 */

namespace gmx {

// Data structure for graphs
class Graph {
public:
    /*
     * Here m is the number of vertices, and xadj[i] is the index in adj of the
     * first element in the ith adjency list. The adjacency lists are stored
     * contiguously in the array adj. This is a very compact way of storing
     * graphs.
     *
     * Example:
     *
     * 0 - 1
     * |
     * 2
     *
     * m = 3
     * xadj = [0, 3, 5, 7]
     * adj = [0, 1, 2, 0, 1, 0, 2]
     *
     */
    int nnodes;
    std::vector<int> xadj;
    std::vector<int> adj;

    /**
     * Constructs an empty graph.
     *
     */
    Graph();

    /**
     * Constructs a graph and reserves space for the nodes and edges.
     *
     * @param nodes The number of nodes of the graph.
     * @param edges The number of edges of the graph.
     */
    Graph(int nodes, int edges);

    /**
     * Get the number of nodes in the graph.
     *
     * @return int The number of nodes in the graph.
     */
    int num_nodes() const;

    /**
     * Get the number of edges in the graph.
     *
     * @return The number of edges in the graph.
     */
    int num_edges() const;

    /**
     * This function assumes that the graph is disjoint and has at least k
     * disjoint subgraphs. It assigns 1 or more subgraphs to each partition
     * while attempting to balance the number of nodes in each partition. The
     * process starts with vertex 0 and continues in order. This function
     * returns an array of num_nodes() vertices in which each element is an
     * integer between 0 and k-1, indicating the partition to which the node
     * belongs.
     *
     * @param k The desired number of partitions.
     * @return An array of length at least the number of nodes in the graph.
     * Each element of the array is an integer between 0 and k-1, indicating
     * the partition that the node belongs to.
     */
    std::vector<int> kway_partition_disjoint(int k) const;

    /**
     * Partitions the graph into k partitions heuristically minimizing the
     * number of edge cuts while attempting to balance the number of nodes in
     * each partition. The quality of the partitions produced by the greedy
     * algorithm used of this function heavily depends on the order of the
     * vertices. Having connected vertices with distant indices will likely
     * result in very unbalanced partitions with many edge cuts. This function
     * returns an array of num_nodes() vertices in which each element is an
     * integer between 0 and k-1, indicating the partition to which the node
     * belongs.
     * @param k The desired number of partitions.
     * @return An array of length at least the number of nodes in the graph.
     * Each element of the array is an integer between 0 and k-1, indicating
     * the partition that the node belongs to.
     */
    std::vector<int> kway_partition(int k) const;

    /**
     * Renumbers the vertices of a graph using a given permutation of length
     * the number of nodes in the graph.
     * The permutation is given as in MATLAB. Example:
     *  p = [2, 1, 0] Means that
     *  Old position 2 is now position 0
     *  Old position 1 is now position 1
     *  Old position 0 is now position 2
     *
     * @param perm The permutation array.
     * @param iperm The inverse permutation array. This is an optional
     * parameter, is an empty vector is provided, it will be computed in this
     * function.
     */
    void renumber_vertices(const std::vector<int> &perm, std::vector<int> &iperm);
};
};   // namespace gmx

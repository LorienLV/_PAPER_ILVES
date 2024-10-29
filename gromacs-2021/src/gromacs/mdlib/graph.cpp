#include "graph.h"

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#if 0
    #include <metis.h>
    #include <mtmetis.h>
    #include <omp.h>
#endif
namespace gmx {

Graph::Graph() : nnodes(0) {}

Graph::Graph(const int nodes, const int edges)
    : nnodes(nodes), xadj(nodes + 1), adj(edges) {
    xadj[0] = 0;
}

int Graph::num_nodes() const { return nnodes; }

int Graph::num_edges() const { return adj.size(); }

std::vector<int> Graph::kway_partition_disjoint(const int k) const {
    std::vector<int> ids(nnodes, -1);

    // Assign atoms to each partition.
    for (int p = 0, node = 0; p < k; ++p) {
        // Heuristic to assign similar number of nodes to each partition.
        const int remaining_nodes = nnodes - node;
        const int remaining_partitions = k - p;
        const int target_size_p = remaining_nodes / remaining_partitions;

        int size = 0;
        while (node < nnodes && size <= target_size_p) {
            if (ids[node] != -1) {
                ++node;
                continue;   // Already assigned.
            }

            ids[node] = p;
            ++size;

            std::function<void(const int)> add_neighbors =
                [&](const int other_node) -> void {
                for (int j = xadj[other_node]; j < xadj[other_node + 1]; ++j) {
                    const int neigh = adj[j];

                    if (ids[neigh] != -1) {
                        continue;   // Already assigned.
                    }

                    ids[neigh] = p;
                    ++size;

                    add_neighbors(neigh);
                }
            };
            add_neighbors(node);   // Keep all the nodes of a block in the same
                                   // partition.

            ++node;
        }
    }

    return ids;
}

std::vector<int> Graph::kway_partition(const int k) const {
    // Try to minimize the cuts, i.e., the number of nodes not assigned to
    // each partition p but connected to it (boundary nodes), while keeping
    // the size of the partition close to the target size.

    // Heuristic parameters.
    // The size of each partition p can be at most 1 + max_above_target times
    // the target size.
    constexpr double max_above_target = 0.05;
    // A cut penalizes as if we were assigning cut_penalty nodes above or below
    // the target size to each partition.
    constexpr int cut_penalty = 10;

    // Auxiliary vector and also the result.
    std::vector<int> ids(nnodes, -1);

    int node = 0;
    for (int p = 0; p < k; ++p) {
        // Heuristic to assign similar number of nodes to each partition.
        const int remaining_nodes = nnodes - node;
        const int remaining_partitions = k - p;
        const int target_size_p = remaining_nodes / remaining_partitions;

        const int max_target_size_p = target_size_p * (1.0 + max_above_target);

        // The best last node assigned to this partition and its score.
        std::pair<int, int> best_cut = {-1, std::numeric_limits<int>::max()};

        int size = 0;
        int ncuts = 0;
        // Note that we do not need to reset ids after each execution of this
        // loop since we are using p to mark the boundary nodes of the p
        // partition.
        while (node < nnodes && size <= max_target_size_p) {
            if (ids[node] == p) {
                // This node was previously a boundary node, but now it is part
                // of p.
                --ncuts;
            }
            else {
                ids[node] = p;
            }
            ++size;

            // Iterate over the neighbors of the node with idx > node.
            // The adjacency list is sorted, so we can use a binary search to
            // find the position of node in the list.
            auto binarySearch =
                [](const int *const vec, const int size, const int target) -> int {
                int left = 0;
                int right = size;

                while (left <= right) {
                    const int mid = left + (right - left) / 2;

                    if (vec[mid] == target) {
                        return mid;
                    }
                    else if (vec[mid] < target) {
                        left = mid + 1;
                    }
                    else {
                        right = mid - 1;
                    }
                }

                return -1;
            };

            // Idx of node in its neighbor list.
            const int node_adj_idx = binarySearch(adj.data() + xadj[node],
                                                  xadj[node + 1] - xadj[node],
                                                  node) +
                                     xadj[node];

            // Find the neighbors of the node with greater idx. These are
            // vertices that are not in the p partition but are connected to it.
            for (int k = node_adj_idx; k < xadj[node + 1]; ++k) {
                const int neigh = adj[k];

                // Check if we have already counted neigh as a boundary node.
                if (ids[neigh] < p) {
                    ids[neigh] = p;
                    ++ncuts;
                }
            }

            const int score = std::abs(size - target_size_p) + ncuts * cut_penalty;

            if (score < best_cut.second) {
                best_cut = {node, score};
            }

            ++node;
        }

        // Reset node the the correct idx based on the best cut.
        node = best_cut.first + 1;
    }
    // Assign the remaining nodes to the last partition.
    // This loop will not be executed in most cases.
    for (; node < nnodes; ++node) {
        ids[node] = k - 1;
    }

    return ids;
}

#if 0
/**
 * This is a wrapper for the k-way partitioning function from METIS.
 * After calling this function, each element of ids will be an integer
 * between 0 and k-1, indicating the partition that the node belongs to.
 *
 * @param k The desired number of partitions.
 * @return An array of length at least the number of nodes in the graph.
 * Each element of the array is an integer between 0 and k-1, indicating
 * the partition that the node belongs to.
 */
std::vector<int> Graph::kway_partition_mtmetis(const int k) const {
    std::cerr << "CALLING METIS\n";
    static_assert(std::is_same_v<unsigned int, mtmetis_vtx_type>,
                  "MTMETIS: mtmetis_vtx_type != int");

    std::vector<int> ids(nnodes);

    if (k == 1) {
        std::fill(ids.begin(), ids.end(), 0);
        return ids;
    }

    // Allocate space for METIS options
    double *options = mtmetis_init_options();

    // Set default options for METIS
    // METIS_SetDefaultOptions(options.data());
    // Try to enforce contiguous partitions
    // It can fail when constraining h-bonds. The graph may be
    // disconnected. options[METIS_OPTION_CONTIG] = 1;
    options[MTMETIS_OPTION_NTHREADS] = 1;

    mtmetis_vtx_type m = nnodes;
    mtmetis_vtx_type idx_one = 1;
    mtmetis_vtx_type k_idx = k;

    mtmetis_wgt_type objval;

    // Call the k-way partitioning function from METIS
    MTMETIS_PartGraphKway(&m,
                          &idx_one,
                          reinterpret_cast<const mtmetis_vtx_type *>(xadj.data()),
                          reinterpret_cast<const mtmetis_vtx_type *>(adj.data()),
                          NULL,
                          NULL,
                          NULL,
                          &k_idx,
                          NULL,
                          NULL,
                          options,
                          &objval,
                          reinterpret_cast<mtmetis_vtx_type *>(ids.data()));

    free(options);

    return ids;
}
#endif

#if 0
std::vector<int32_t> Graph::kway_partition_metis(const int k) {
    // Allocate space to store the identification number for each array.
    static_assert(std::is_same_v<int32_t, idx_t>, "METIS: idx_t != int32_t");
    std::vector<int32_t> id(nnodes);

    // METIS *will* crash if we try to do a 1-way partition.
    // A logical check is made to prevent this from happening
    if (k < 2) {
        std::fill(id.begin(), id.end(), 0);
        return id;
    }

    idx_t objval;

    // Allocate space for METIS options
    std::vector<idx_t> options(METIS_NOPTIONS);

    // Set default options for METIS
    METIS_SetDefaultOptions(options.data());
    // Try to enforce contiguous partitions
    // It can fail when constraining h-bonds. The graph may be
    // disconnected. options[METIS_OPTION_CONTIG] = 1;

    idx_t m = nnodes;
    idx_t idx_one = 1;
    idx_t k_idx = k;

    // Call the k-way partitioning function from METIS
    METIS_PartGraphKway(&m,
                        &idx_one,
                        xadj.data(),
                        adj.data(),
                        NULL,
                        NULL,
                        NULL,
                        &k_idx,
                        NULL,
                        NULL,
                        options.data(),
                        &objval,
                        id.data());

    return id;
}
#endif

void Graph::renumber_vertices(const std::vector<int> &perm,
                              std::vector<int> &iperm) {
    // Compute the inverse permutation if it is not provided.
    if (iperm.empty()) {
        iperm.resize(nnodes);
        for (int node = 0; node < nnodes; ++node) {
            iperm[perm[node]] = node;
        }
    }

    Graph perm_graph(nnodes, adj.size());

    for (int i = 0, edge = 0; i < nnodes; ++i) {
        const int node = perm[i];   // The ith node is the old perm[i] node.

        const int old_edge = edge;
        for (int k = xadj[node]; k < xadj[node + 1]; ++k, ++edge) {
            // The ith node is connected to iperm[adj[k]].
            perm_graph.adj[edge] = iperm[adj[k]];
        }
        // The adjacency list of node may not be in order.
        std::sort(perm_graph.adj.begin() + old_edge, perm_graph.adj.begin() + edge);

        perm_graph.xadj[i + 1] = edge;
    }

    *this = std::move(perm_graph);
}

};   // namespace gmx

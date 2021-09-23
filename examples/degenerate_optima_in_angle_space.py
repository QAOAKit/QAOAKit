import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from QAOAKit import get_3_reg_dataset_table, get_3_reg_dataset_table_row

"""
Compute typical degenerate optima for 3 regular QAOA on the ensemble of all 3 regular graphs
 of leq 16 vertices, using K-means clustering.
"""


"""
First, import all degenerate p=1 and 2 angles from the dataset
"""
# The pandas database of all 3 regular graphs
tab = get_3_reg_dataset_table().reset_index()

# Database is ordered by p=1...10 for each graph, so requesting every 10th yields p=1 angles
p1_degenerate_gamma = np.stack(
    tab[tab["p_max"] == 1]["all gamma (degenerate optima)"].to_numpy()
)
p1_degenerate_beta = (
    np.stack(tab[tab["p_max"] == 1]["all beta (degenerate optima)"].to_numpy()) % 0.5
)

p2_degenerate_gamma = np.stack(
    tab[tab["p_max"] == 2]["all gamma (degenerate optima)"].to_numpy()
)
p2_degenerate_beta = (
    np.stack(tab[tab["p_max"] == 2]["all beta (degenerate optima)"].to_numpy()) % 0.5
)

# Sanity check
assert p1_degenerate_gamma.shape == (4681, 4, 1)
assert p1_degenerate_beta.shape == (4681, 4, 1)
assert p2_degenerate_gamma.shape == (4681, 8, 2)
assert p2_degenerate_beta.shape == (4681, 8, 2)


"""
Flatten the data to a single coordinate
"""
p1_thetas = np.concatenate((p1_degenerate_beta, p1_degenerate_gamma), axis=2)
p2_thetas = np.concatenate((p2_degenerate_beta, p2_degenerate_gamma), axis=2)

p1_thetas = p1_thetas.reshape(4681 * 4, 2)
p2_thetas = p2_thetas.reshape(4681 * 8, 4)


"""
Poor-mans K-means: initialize cluster centers by the fixed angles,
 which are optimal for the smallest graph with no cycles leq 5. Conveneniently, this is a Cage graphs,
 specifically the heawood graph, which has a networkx generator =).
Then, assign each angle to a particular cluster by euclidian distance
Finally, optimize to minimize sum of squares variance ("inertia")
"""


# Look up the data for the Heawood graph
G = nx.heawood_graph()
row1 = get_3_reg_dataset_table_row(G, 1)
row2 = get_3_reg_dataset_table_row(G, 2)

# Format angles appropriately
p1_fixed_gamma = np.array(row1["all gamma (degenerate optima)"])
p1_fixed_beta = np.array(row1["all beta (degenerate optima)"]) % 0.5
p2_fixed_gamma = np.array(row2["all gamma (degenerate optima)"])
p2_fixed_beta = np.array(row2["all beta (degenerate optima)"]) % 0.5
p1_fixed_thetas = np.concatenate((p1_fixed_beta, p1_fixed_gamma), axis=1)
p2_fixed_thetas = np.concatenate((p2_fixed_beta, p2_fixed_gamma), axis=1)


# K means clustering: cluster by distance to naive center
# Establish distance between naive centers and each point
dists1 = (
    (p1_thetas[:, 0].reshape(4681 * 4, 1) - p1_fixed_thetas[:, 0].reshape(1, 4))
) ** 2
dists1 += (
    (p1_thetas[:, 1].reshape(4681 * 4, 1) - p1_fixed_thetas[:, 1].reshape(1, 4))
) ** 2
dists1 = np.sqrt(dists1)

dists2 = (
    (p2_thetas[:, 0].reshape(4681 * 8, 1) - p2_fixed_thetas[:, 0].reshape(1, 8))
) ** 2
dists2 += (
    (p2_thetas[:, 1].reshape(4681 * 8, 1) - p2_fixed_thetas[:, 1].reshape(1, 8))
) ** 2
dists2 += (
    (p2_thetas[:, 2].reshape(4681 * 8, 1) - p2_fixed_thetas[:, 2].reshape(1, 8))
) ** 2
dists2 += (
    (p2_thetas[:, 3].reshape(4681 * 8, 1) - p2_fixed_thetas[:, 3].reshape(1, 8))
) ** 2
dists2 = np.sqrt(dists2)

# Assign cluster labels by minimum distance
cluster_labels_1 = np.argmin(dists1, axis=1)
cluster_labels_2 = np.argmin(dists2, axis=1)


# Minimization of SoS error is equivalently setting the center point equal to the average
average1 = []
for ind in range(4):
    average1.append(np.average(p1_thetas[cluster_labels_1 == ind, :], 0))
average1 = np.array(average1)
average2 = []
for ind in range(8):
    average2.append(np.average(p2_thetas[cluster_labels_2 == ind, :], 0))
average2 = np.array(average2)


"""
Visualize clustering
"""

# p=1 only has two axes
plt.rcParams.update({"font.size": 16})
plt.figure(figsize=(8, 6))
plt.subplots_adjust(left=0.15, right=0.85, bottom=0.15)
ax = plt.subplot(1, 1, 1)

plt.scatter(
    p1_thetas[:, 0], p1_thetas[:, 1], c=cluster_labels_1, s=10, edgecolors="None"
)
plt.scatter(average1[:, 0], average1[:, 1], c="b", s=500, marker="+")

plt.xlabel("$\\beta/\\pi$")
plt.ylabel("$\\gamma/\\pi$")
plt.savefig("p=1_optima.pdf")


# p=2 has 4 axes, so we project along 6 orthogonal planes
labels = ["$\\beta_1/\\pi$", "$\\beta_2/\\pi$", "$\\gamma_1/\\pi$", "$\\gamma_2/\\pi$"]
for i0 in range(4):
    for i1 in range(i0 + 1, 4):
        plt.figure(figsize=(8, 6))
        plt.subplots_adjust(left=0.15, right=0.85, bottom=0.15)
        ax = plt.subplot(1, 1, 1)
        plt.scatter(
            p2_thetas[:, i0],
            p2_thetas[:, i1],
            c=cluster_labels_2,
            s=10,
            edgecolors="None",
        )
        plt.scatter(average2[:, i0], average2[:, i1], c="b", s=500, marker="+")

        plt.xlabel(labels[i0])
        plt.ylabel(labels[i1])

        plt.savefig("p=2_optima_{:0.0f}_{:0.0f}.pdf".format(i0, i1))


"""
Print centroid angles
"""
print("p=1 degenerate angles:")
for ind in range(4):
    print("Beta:", average1[ind, 0], "Gamma:", average1[ind, 1])

print("")
print("p=2 degenerate angles:")
for ind in range(8):
    print("Beta:", average2[ind, 0:2], "Gamma:", average2[ind, 2::])

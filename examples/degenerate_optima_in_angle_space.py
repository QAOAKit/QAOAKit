import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from matplotlib import rc

rc("font", **{"family": "serif", "serif": ["Times"]})
rc("text", usetex=True)

import seaborn as sns
import sklearn.cluster
from QAOAKit import get_3_reg_dataset_table


"""
Compute typical degenerate optima for 3 regular QAOA on the ensemble of all 3 regular graphs
 of leq 16 vertices, using K-means clustering.
"""


# The pandas database of all 3 regular graphs
tab = get_3_reg_dataset_table().reset_index()

# Get the degenerate angles
p1_degenerate_gamma = np.stack(
    tab[tab["p_max"] == 1]["all gamma (degenerate optima)"].to_numpy()
)
p1_degenerate_beta = np.stack(
    tab[tab["p_max"] == 1]["all beta (degenerate optima)"].to_numpy()
)

p2_degenerate_gamma = np.stack(
    tab[tab["p_max"] == 2]["all gamma (degenerate optima)"].to_numpy()
)
p2_degenerate_beta = np.stack(
    tab[tab["p_max"] == 2]["all beta (degenerate optima)"].to_numpy()
)

#
p1_degenerate_gamma = (p1_degenerate_gamma + 1) % 2 - 1
p2_degenerate_gamma = (p2_degenerate_gamma + 1) % 2 - 1
p1_degenerate_beta = (p1_degenerate_beta + 0.25) % 0.5 - 0.25
p2_degenerate_beta = (p2_degenerate_beta + 0.25) % 0.5 - 0.25

# Sanity check
assert p1_degenerate_gamma.shape == (4681, 4, 1)
assert p1_degenerate_beta.shape == (4681, 4, 1)
assert p2_degenerate_gamma.shape == (4681, 8, 2)
assert p2_degenerate_beta.shape == (4681, 8, 2)


# Flatten the data to a single coordinate
p1_thetas = np.concatenate((p1_degenerate_beta, p1_degenerate_gamma), axis=2)
p2_thetas = np.concatenate((p2_degenerate_beta, p2_degenerate_gamma), axis=2)
p1_thetas = p1_thetas.reshape(4681 * 4, 2)
p2_thetas = p2_thetas.reshape(4681 * 8, 4)


# Execute K means clustering
skl1 = sklearn.cluster.KMeans(n_clusters=4)
skl2 = sklearn.cluster.KMeans(n_clusters=8)
skl1.fit(p1_thetas)
skl2.fit(p2_thetas)


"""
Visualize clustering
"""
# p=1 only has two axes
# plt.rcParams.update({"font.size": 16})
plt.figure(figsize=(1.5 * 6.92654 / 2, 6.92654 / 2), dpi=500)
# plt.subplots_adjust(left=0.15, right=0.85, bottom=0.15)
ax = plt.subplot(1, 1, 1)

plt.scatter(
    p1_thetas[:, 0],
    p1_thetas[:, 1],
    c=skl1.labels_,
    s=10,
    edgecolors="None",
    rasterized=True,
)
plt.scatter(
    skl1.cluster_centers_[:, 0], skl1.cluster_centers_[:, 1], c="r", s=500, marker="+"
)

plt.xlabel("$\\beta/\\pi$")
plt.ylabel("$\\gamma/\\pi$")
plt.axis([-0.25, 0.25, -1, 1])
plt.tight_layout()
plt.savefig("p=1_optima.pdf")
plt.savefig("p=1_optima.png")


# p=2 has 4 axes, so we project along 6 orthogonal planes
labels = ["$\\beta_1/\\pi$", "$\\beta_2/\\pi$", "$\\gamma_1/\\pi$", "$\\gamma_2/\\pi$"]
axees = [[-0.25, 0.25], [-0.25, 0.25], [-1, 1], [-1, 1]]
for i0 in range(4):
    for i1 in range(i0 + 1, 4):
        plt.figure(figsize=(1.5 * 6.92654 / 2, 6.92654 / 2), dpi=500)
        plt.subplots_adjust(left=0.15, right=0.85, bottom=0.15)
        ax = plt.subplot(1, 1, 1)
        plt.scatter(
            p2_thetas[:, i0],
            p2_thetas[:, i1],
            c=skl2.labels_,
            s=10,
            edgecolors="None",
            rasterized=True,
        )
        plt.scatter(
            skl2.cluster_centers_[:, i0],
            skl2.cluster_centers_[:, i1],
            c="r",
            s=500,
            marker="+",
        )

        plt.xlabel(labels[i0])
        plt.ylabel(labels[i1])
        plt.axis(axees[i0] + axees[i1])
        plt.tight_layout()
        plt.savefig("p=2_optima_{:0.0f}_{:0.0f}.pdf".format(i0, i1))
        plt.savefig("p=2_optima_{:0.0f}_{:0.0f}.png".format(i0, i1))


"""
Print centroid angles
"""
print("p=1 degenerate angles:")
for ind in range(4):
    print(
        "Beta:", skl1.cluster_centers_[ind, 0], "Gamma:", skl1.cluster_centers_[ind, 1]
    )
print("")
print("p=1 RMS distance:", np.sqrt(skl1.inertia_ / (4681 * 4)))
print("")
print("p=2 degenerate angles:")
for ind in range(8):
    print(
        "Beta:",
        skl2.cluster_centers_[ind, 0:2],
        "Gamma:",
        skl2.cluster_centers_[ind, 2::],
    )
print("")
print("p=2 RMS distance:", np.sqrt(skl2.inertia_ / (4681 * 8)))

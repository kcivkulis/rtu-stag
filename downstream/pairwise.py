import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt


def bray_curtis(u, v):
    assert(len(u) == len(v))
    c = sum(min(x, y) for x, y in zip(u, v))
    return 1 - 2 * c / (sum(u) + sum(v))

genus_to_id = {}
genus_count = 0


def read_abundance_vector(filename):
    v = [0] * genus_count
    with open(filename) as f:
        first = True
        for line in f:
            if not first:
                data = line.rstrip().split('\t')
                genus = data[0]
                cnt = float(data[5])
                v[genus_to_id[genus]] = cnt
            first = False
    return v


for filename in os.listdir("./bracken_1"):
    if not filename.endswith(".bracken"):
        continue
    with open(os.path.join("bracken_1", filename)) as f:
        first = True
        for line in f:
            if not first:
                genus = line.rstrip().split('\t')[0]
                if not genus in genus_to_id:
                    genus_to_id[genus] = genus_count
                    genus_count += 1
            first = False
    

abundances_before = []
abundances_after = []

with open("patient_richness.csv") as f:
    first = True
    for line in f:
        if not first:
            data = line.rstrip().split(',')
            sample_name_before = data[0] + "_" + data[1] + "_" + "T1"
            sample_name_after = data[0] + "_" + data[1] + "_" + "T2"

            before = read_abundance_vector(os.path.join("bracken_1", sample_name_before + ".bracken"))
            after = read_abundance_vector(os.path.join("bracken_1", sample_name_after + ".bracken"))

            abundances_before.append([data[0], before])
            abundances_after.append([data[0], after])
        first = False

n = len(abundances_before)

d_before, d_after = np.zeros(shape=(n, n)), np.zeros(shape=(n, n))


for i in range(n):
    for j in range(n):
        d_before[i][j] = bray_curtis(abundances_before[i][1], abundances_before[j][1])
        d_after[i][j]  = bray_curtis(abundances_after[i][1], abundances_after[j][1])


fig, ax = plt.subplots(2, 1, figsize=(5, 10))

g1 = sns.heatmap(d_before, linewidth = 0.5, ax=ax[0], square=True)
g1.set(xlabel = None, ylabel = None)
g1.set(xticklabels = [], yticklabels = [])
g1.tick_params(left=False, bottom=False)
g1.set(title="Before")

g2 = sns.heatmap(d_after, linewidth = 0.5, ax=ax[1], square=True)
g2.set(xlabel = None, ylabel = None)
g2.set(xticklabels = [], yticklabels = [])
g2.tick_params(left=False, bottom=False)
g2.set(title="After")

plt.tight_layout()
plt.savefig('heatmap.png')

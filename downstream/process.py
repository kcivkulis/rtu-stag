import json

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from skbio.diversity import alpha


# Utils

def strip_otu_levels(list):
    return tuple(x[1] for x in list)


def map_sorted_samples(f, samples):
    res = {}
    for treatment in ["Control", "STD2", "STD3"]:
        res[treatment] = {}
        for time in ["before", "after"]:
            res[treatment][time] = {}
            for sample in samples[treatment][time]:
                res[treatment][time][sample.id] = f(sample)
    return res


def vectorize(sample, keys, attr):
    m = getattr(sample, attr)
    return [m[x] if x in m else 0 for x in keys]


def get_abundances_including_other(sample, category_names, attr):
    m = getattr(sample, attr)
    v = vectorize(sample, category_names, attr)
    v.append(sum(m.values()) - sum(v))
    return np.array(v)


def get_all_sample_keys(samples, attr):
    all_keys = set()
    for sample in samples:
        if getattr(sample, attr) is not None:
            for k in getattr(sample, attr):
                all_keys.add(k)
    return sorted(all_keys)

# End of utils


def parse_lists(lists, list_parser):
    result = {}
    for list in lists:
        x = list_parser(list[:-1])
        y = int(list[-1])
        result[x] = y
    return result


class Sample:
    def __init__(self, id, patient):
        self.id = id
        self.patient = patient
        with open("outputs/samples/" + id + ".json") as f:
            data = json.load(f)
            self.resistome = parse_lists(data["resistome"],
                                         lambda list: tuple(list))
            self.taxonomy = parse_lists(data["taxonomy"],
                                        lambda list: tuple(tuple(p) for p in list))
            self.read_count = int(data["read_count"])


def get_sample(sample_id, self):
    if type(sample_id) == str:
        return Sample(sample_id, self)
    else:
        return None


class Patient:
    def __init__(self, id, sample_id_before, sample_id_after, treatment):
        self.id = id
        self.sample_before = get_sample(sample_id_before, self)
        self.sample_after = get_sample(sample_id_after, self)
        self.treatment = treatment


def read_patient_data(filename):
    table = pd.read_csv(filename, index_col="patient_id")
    patients = []
    for index, row in table.iterrows():
        patients.append(Patient(index,
                                row["sample_id_before"],
                                row["sample_id_after"],
                                row["treatment"]))
    return patients


def get_all_samples(patients):
    samples = []
    for p in patients:
        if p.sample_before is not None:
            samples.append(p.sample_before)
        if p.sample_after is not None:
            samples.append(p.sample_after)
    return samples


def sort_samples(samples):
    sorted_samples = {}
    sorted_samples["Control"] = {}
    sorted_samples["STD2"] = {}
    sorted_samples["STD3"] = {}
    for treatment in sorted_samples:
        sorted_samples[treatment]["before"] = []
        sorted_samples[treatment]["after"] = []

    for s in samples:
        if s.patient.sample_before == s:
            time = "before"
        else:
            time = "after"
        sorted_samples[s.patient.treatment][time].append(s)
    return sorted_samples


# Drawings

def draw_resistome_amount(samples, filename):
    def f(samples):
        return list(filter(lambda x: x < 0.08,
                    map(lambda s: sum(s.resistome.values()) / s.read_count,
                        samples)))

    fig, ax = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(4, 6))

    for i, treatment in enumerate(["Control", "STD2", "STD3"]):
        for j, time in enumerate(["before", "after"]):
            ax[i, j].violinplot(f(samples[treatment][time]),
                                showmeans=True,
                                showextrema=False)

    fig.suptitle("Percentage of reads mapped to resistome")

    ax[1, 0].set_xticks([])

    ax[2, 0].set_xlabel("Before")
    ax[2, 1].set_xlabel("After")

    ax[0, 0].set_ylabel("Control")
    ax[1, 0].set_ylabel("STD2")
    ax[2, 0].set_ylabel("STD3")

    plt.savefig(filename)


def draw_beta_diversity(samples, filename):
    all_otus = get_all_sample_keys(samples, "taxonomy")
    m = np.array([vectorize(s, all_otus, "taxonomy") for s in samples])
    d = pdist(m, "jensenshannon")
    fig, ax = plt.subplots(figsize=(50, 10))
    dendrogram(linkage(d, method="average"),
               labels=[s.id for s in samples],
               ax=ax,
               orientation="top",
               color_threshold=0)
    ax.set_yticks([])
    plt.savefig(filename, bbox_inches="tight")


def draw_horizontal_bars(bar_widths, category_names, suptitle, filename):
    fig, axs = plt.subplots(2, 3, figsize=(60, 30))
    fig.suptitle(suptitle, fontsize=40)
    axs = axs.transpose()

    category_colors = plt.get_cmap('hsv')(np.linspace(0.05,
                                                      0.95,
                                                      len(category_names)))

    def draw_chart(d, labels, ax):
        d_cum = d.cumsum(axis=1)
        for i, (name, color) in enumerate(zip(category_names, category_colors)):
            widths = d[:, i]
            starts = d_cum[:, i] - widths
            ax.barh(labels, widths, left=starts, height=0.5, label=name, color=color)

    axs[0][0].set_ylabel("Before", fontsize=40)
    axs[0][1].set_ylabel("After", fontsize=40)

    for treatment, ax_col in zip(["Control", "STD2", "STD3"], axs):
        ax_col[0].set_title(treatment, fontsize=40)

        for time, ax in zip(["before", "after"], ax_col):
            d, labels = [], []
            for id in bar_widths[treatment][time]:
                labels.append(id)
                s = sum(bar_widths[treatment][time][id])
                d.append([x / s for x in bar_widths[treatment][time][id]])
            draw_chart(np.array(d), labels, ax)

    axs[0][0].legend(ncol=len(category_names), bbox_to_anchor=(0, 1.2), loc='lower left', fontsize='small')

    plt.savefig(filename)
    plt.close(fig)


def draw_enterotypes(samples, significant_otus, filename):
    def sample_to_vector(sample):
        return get_abundances_including_other(sample, significant_otus, "taxonomy")

    draw_horizontal_bars(map_sorted_samples(sample_to_vector, samples),
                         [strip_otu_levels(x[-2:]) for x in significant_otus] + ["Other"],
                         "Taxonomy",
                         filename)

# For tables


def print_attr_table(samples, output_file, key_pretifier, attr):
    keys = get_all_sample_keys(samples, attr)
    with open(output_file, "w") as f:
        print("sample_id", *(key_pretifier(k) for k in keys), sep=",", file=f)
        for s in samples:
            if getattr(s, attr) is not None:
                print(s.id, *vectorize(s, keys, attr), sep=",", file=f)


def print_taxonomy_table(samples, output_file):
    print_attr_table(samples, output_file, lambda k: k[-1][1], "taxonomy")


def print_resistome_table(samples, output_file):
    print_attr_table(samples, output_file, lambda r: r[-1], "resistome")


def get_significant_keys(samples, attr, offset):
    significant = []

    for sample in samples:
        m = getattr(sample, attr)
        if m is None:
            continue
        total_reads = sum(m.values())

        while True:
            s = sum(vectorize(sample, significant, attr))
            if s > offset * total_reads:
                break
            i_take = None
            for i in m:
                if i in significant:
                    continue
                if i_take is None or m[i] > m[i_take]:
                    i_take = i
            significant.append(i_take)

    return significant


def get_significant_otus(samples):
    return get_significant_keys(samples, "taxonomy", 0.1)


def get_significant_resistance_mechanisms(samples):
    return get_significant_keys(samples, "resistome", 0.7)


def calculate_enterotype(sample, significant_otus):
    if sample.taxonomy is None:
        return None
    total = sum(sample.taxonomy.values())
    max_otu = None
    for otu in filter(lambda otu: otu in sample.taxonomy, significant_otus):
        if max_otu is None or sample.taxonomy[otu] > sample.taxonomy[max_otu]:
            max_otu = otu
    if sample.taxonomy[max_otu] > total / 5:
        max_otu = max_otu[-1][1]
        if max_otu == "Bacteroides":
            return "A"
        elif max_otu == "Prevotella":
            return "B"
        else:
            return "C"
    return "D"


def get_patient_table(patients):
    def skip_if_none(f, x):
        return f(x) if x is not None else None

    def get_metric_applier(metric):
        return lambda sample: metric(list(sample.taxonomy.values()))

    def append_both(name, func):
        columns.append((name + "_before",
                        lambda patient: skip_if_none(func,
                                                     patient.sample_before)))
        columns.append((name + "_after",
                        lambda patient: skip_if_none(func,
                                                     patient.sample_after)))

    metrics = [
        (alpha.chao1, "chao1"),
        (alpha.shannon, "shannon"),
        (alpha.observed_otus, "observed"),
        (alpha.ace, "ace"),
        (alpha.berger_parker_d, "berger_parker"),
    ]

    columns = []

    for metric, metric_name in metrics:
        append_both(metric_name, get_metric_applier(metric))
    append_both("enterotype", lambda sample: sample.enterotype)
    append_both("read_count", lambda sample: sample.read_count)
    columns.append(("treatment", lambda patient: patient.treatment))

    all_data = []

    for patient in patients:
        all_data.append([f(patient) for _, f in columns])

    df = pd.DataFrame(columns=[name for name, _ in columns],
                      index=[p.id for p in patients],
                      data=all_data)
    df.index.name = "patient_id"
    return df


if __name__ == "__main__":
    all_patients = read_patient_data("outputs/metadata.csv")
    samples = get_all_samples(all_patients)
    sorted_samples = sort_samples(samples)

    print_taxonomy_table(samples, "outputs/tables/taxonomy.csv")
    print_resistome_table(samples, "outputs/tables/resistome.csv")

    significant_otus = get_significant_otus(samples)

    for sample in samples:
        sample.enterotype = calculate_enterotype(sample, significant_otus)

    get_patient_table(all_patients).to_csv("outputs/tables/patients.csv")

    draw_enterotypes(sorted_samples, significant_otus, "outputs/pictures/enterotypes.svg")
    draw_resistome_amount(sorted_samples, "outputs/pictures/resistome_change.png")
    draw_beta_diversity(samples, "outputs/pictures/beta_diversity.png")

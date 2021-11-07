import json

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.stats import wilcoxon
from skbio.diversity import alpha
import seaborn as sns


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


# AMRplusplus

def parse_amrplusplus_name(s):
    name = s.split('|')
    if name[-1] == "RequiresSNPConfirmation":
        name.pop()
    name.append(name[0])
    name.pop(0)
    return tuple(name)


def get_amr_gene_lengths():
    amr_gene_lengths = {}
    filename = "megares_full_database_v2.00.fasta"
    with open(filename) as f:
        name = None
        length = 0
        for line in f:
            line = line.rstrip()
            if line[0] == '>':
                if name is not None:
                    amr_gene_lengths[name] = length
                length = 0
                name = parse_amrplusplus_name(line[1:])
            else:
                length += len(line)
        amr_gene_lengths[name] = length
    return amr_gene_lengths


# END

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

def draw_violin_plots(patients, f, skip_condition, text_if_skipped, suptitle, ylabel, filename):
    treatment, time, value = [], [], []
    for patient in patients:
        if patient.sample_before is None or patient.sample_after is None:
            print("Skipping", patient.id, ": missing samples")
            continue
        if skip_condition(patient.sample_before):
            print("Skipping", patient.id, " (before):", text_if_skipped(patient.sample_before))
            continue
        if skip_condition(patient.sample_after):
            print("Skipping", patient.id, " (after):", text_if_skipped(patient.sample_after))
        t = patient.treatment
        if t == "STD2":
            t = "Amoxicillinum & Clarithromycinum"
        elif t == "STD3":
            t = "Amoxicillinum"
        treatment.append(t)
        treatment.append(t)
        time.append("Before")
        time.append("After")
        value.append(f(patient.sample_before))
        value.append(f(patient.sample_after))

    df = pd.DataFrame({"treatment": treatment, "time": time, "value": value})

    fig, ax = plt.subplots(figsize=(10, 3))
    sns.violinplot(x="treatment", y="value", hue="time", data=df, palette="muted", split=True, linewidth=0.5, inner="stick", ax=ax)
    fig.suptitle(suptitle)
    ax.set_ylabel(ylabel)
    ax.set_xlabel("")
    ax.get_legend().set_title("")
    fig.savefig(filename, bbox_inches="tight")


def draw_changes_before_after(patients, f, suptitle, ylabel, filename):
    fig, axs = plt.subplots(1, 3, sharey=True, figsize=(8, 15))
    for ax, treatment in zip(axs, ["Control", "STD2", "STD3"]):
        for p in patients:
            if p.treatment != treatment:
                continue
            if p.sample_before is None or p.sample_after is None:
                print("Skipping", p.id, ": missing samples")
                continue
            ax.plot([0, 1], [f(p.sample_before), f(p.sample_after)], color="black", marker="o", linewidth=0.2, markersize=0.5)

    for ax, treatment in zip(axs, ["Control", "Amoxicillinum & Clarithromycinum", "Amoxicillinum"]):
        ax.set_xticklabels([])
        ax.set_xticks([])
        ax.set_xlabel(treatment)
    fig.suptitle(suptitle)
    axs[0].set_ylabel(ylabel)
    fig.savefig(filename, bbox_inches="tight")


def draw_and_print_resistome_amount(patients, filename, filename_changes, filename_percentages):
    def f(sample):
        return sum(sample.resistome.values()) / sample.read_count * 100

    def skip_condition(sample):
        return f(sample) > 8

    def text_if_skipped(sample):
        return "Resistome is suspiciously large ({:2.2f})".format(f(sample))

    suptitle = "Percentage of reads mapped to resistome"
    ylabel = "Mapped reads, %"

    draw_violin_plots(patients, f, skip_condition, text_if_skipped, suptitle, ylabel, filename)
    draw_changes_before_after(patients, f, suptitle, ylabel, filename_changes)

    before, after = {}, {}
    before["Control"] = []
    before["STD2"] = []
    before["STD3"] = []
    after["Control"] = []
    after["STD2"] = []
    after["STD3"] = []

    for p in patients:
        if p.sample_before:
            before[p.treatment].append(f(p.sample_before))
        if p.sample_after:
            after[p.treatment].append(f(p.sample_after))

    with open(filename_percentages, "w") as out:
        for t in ["Control", "STD2", "STD3"]:
            print(t + " before", file=out)
            print(*before[t], file=out)
            print("mean:", np.mean(before[t]), "stddev:", np.std(before[t]), file=out)
            print(t + " after", file=out)
            print(*after[t], file=out)
            print("mean:", np.mean(after[t]), "stddev:", np.std(after[t]), file=out)


def draw_observed_amr_amount_single(patients, prefix_length, level_name, filename_prefix):
    def f(sample):
        return len({x[:prefix_length] for x in sample.resistome})

    suptitle = "Number of different observed AMR " + level_name
    ylabel = "Count"

    draw_violin_plots(patients, f, lambda sample: False, lambda sample: "", suptitle, ylabel, filename_prefix + "_violin.svg")
    draw_changes_before_after(patients, f, suptitle, ylabel, filename_prefix + "_change.svg")


def draw_observed_amr_amount(patients, filename_prefix):
    for prefix_count, level_name in [(1, "types"), (2, "classes"), (3, "mechanisms"), (4, "genes")]:
        draw_observed_amr_amount_single(patients, prefix_count, level_name, filename_prefix + "_" + level_name)


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


def draw_resistome_profile(samples):
    amr_lengths = get_amr_gene_lengths()

    def read_count_to_relative_abundance(sample):
        resistome = sample.resistome
        total_reads = sum(resistome[amr] for amr in resistome)
        m = {}
        for amr in resistome:
            m[amr] = resistome[amr] / (amr_lengths[amr] / 1e3 * total_reads / 1e6)
        s = sum(m[amr] for amr in m)
        for amr in resistome:
            m[amr] /= s

        return m

    abundance = map_sorted_samples(read_count_to_relative_abundance, samples)

    avg_profile = {}

    for group in ["Control", "STD2", "STD3"]:
        avg_profile[group] = {}

        for time in ["before", "after"]:
            n = len(abundance[group][time])
            m = {}
            for sample in abundance[group][time]:
                for amr in abundance[group][time][sample]:
                    if amr not in m:
                        m[amr] = 0
                    m[amr] += abundance[group][time][sample][amr]
            for amr in m:
                m[amr] /= n

            avg_profile[group][time] = m

    fig, ax = plt.subplots(figsize=(5, 20))

    for group in ["Control", "STD2", "STD3"]:
        for time in ["before", "after"]:
            new_map = {}
            for amr in avg_profile[group][time]:
                if amr[:3] not in new_map:
                    new_map[amr[:3]] = 0
                new_map[amr[:3]] += avg_profile[group][time][amr]

            print(new_map)

            print(new_map.keys())

            ax.barh(list(map(lambda s: s[-1], new_map.keys())), list(new_map.values()))
            plt.show()


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


def statistics_for_groups(patients):
    def amount_of_penicillin_binding_protein(sample):
        total = 0
        for r in sample.resistome:
            if r[2] == "Penicillin_binding_protein":
                total += sample.resistome[r]
        return total / sample.read_count * 100

    def relative_amount_of_penicillin_binding_protein(sample):
        total = 0
        for r in sample.resistome:
            if r[2] == "Penicillin_binding_protein":
                total += sample.resistome[r]
        return total / sum(sample.resistome.values()) * 100

    def relative_amount_of_amr_reads(sample):
        return sum(sample.resistome.values()) / sample.read_count * 100

    functions = [(amount_of_penicillin_binding_protein, "Amount of penicillin binding protein relative to all reads, %"),
                 (relative_amount_of_penicillin_binding_protein, "Amount of penicillin binding protein relative all resistome reads, %"),
                 (relative_amount_of_amr_reads, "AMR count relative to all reads, %")]

    for func, func_name in functions:
        print("Calculating paired wilcoxon for:", func_name)

        for treatment in ["Control", "STD2", "STD3"]:
            print("Group:", treatment)
            x = []
            y = []
            for p in patients:
                if p.treatment != treatment:
                    continue
                if p.sample_before is None or p.sample_after is None:
                    continue
                if func == relative_amount_of_amr_reads:
                    if func(p.sample_before) > 8 or func(p.sample_after) > 8:
                        continue
                    if func(p.sample_after) > 1.3:
                        print(p.id)
                x.append(func(p.sample_before))
                y.append(func(p.sample_after))

            print("Before: {:7f}".format(np.mean(x)), "After: {:7f}".format(np.mean(y)), wilcoxon(x, y))
        print()


# Rarefaction

def rarefy(v, f):
    X = []
    Y = []
    step = 5
    sizes = [sum(v) * f // 100 for f in range(step, 101, step)]
    arr = []
    for i, x in enumerate(v):
        arr += [i] * x
    np.random.shuffle(arr)

    cur = [0] * len(v)

    j = 0
    for i, x in enumerate(arr):
        cur[x] += 1
        if i + 1 == sizes[j]:
            X.append(sum(cur))
            Y.append(f(cur))
            j += 1
    return X, Y


def draw_rarefaction(samples, filename, filename_text):
    fig, ax = plt.subplots(figsize=(15, 15))
    ax.grid()
    ax.ticklabel_format(axis='x', style='plain')
    ax.xaxis.set_tick_params(rotation=90)
    ax.set_xlabel("Read count")
    ax.set_ylabel("Observed species")

    with open(filename_text, "w") as out:
        for s in samples:
            x, y = rarefy(list(s.taxonomy.values()), alpha.observed_otus)
            ax.plot(x, y, linewidth=0.12)
            ax.text(x[-1], y[-1] - 2, s.id, fontsize=2)
            print(s.id, file=out)
            for a, b in zip(x, y):
                print(a, b, file=out)

    ax.set_xlim(0, None)
    ax.set_ylim(0, None)

    plt.savefig(filename, bbox_inches='tight')
    plt.close(fig)


if __name__ == "__main__":
    all_patients = read_patient_data("outputs/metadata.csv")
    samples = get_all_samples(all_patients)
    sorted_samples = sort_samples(samples)

    draw_resistome_profile(sorted_samples)

    statistics_for_groups(all_patients)

    print_taxonomy_table(samples, "outputs/tables/taxonomy.csv")
    print_resistome_table(samples, "outputs/tables/resistome.csv")

    significant_otus = get_significant_otus(samples)

    for sample in samples:
        sample.enterotype = calculate_enterotype(sample, significant_otus)

    get_patient_table(all_patients).to_csv("outputs/tables/patients.csv")

    draw_enterotypes(sorted_samples, significant_otus, "outputs/pictures/enterotypes.svg")
    draw_and_print_resistome_amount(all_patients, "outputs/pictures/resistome_change.svg", "outputs/pictures/resistome_change_2.svg", "outputs/tables/resistome_percents.txt")
    draw_observed_amr_amount(all_patients, "outputs/pictures/observed_amr")
    draw_beta_diversity(samples, "outputs/pictures/beta_diversity.png")

    draw_rarefaction(samples, "outputs/pictures/taxonomy_rarefaction.svg", "outputs/tables/taxonomy_rarefaction.txt") # This is slow

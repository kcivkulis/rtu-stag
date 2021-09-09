import os
import subprocess
import json


outputs_dir = os.path.join(os.getenv("GROUP"), "projects_real", "ERAF184")
bracken_db_path = os.path.join(os.getenv("GROUP"),
                               "databases",
                               "full_ref_bafp",
                               "database150mers.kmer_distrib")


def get_all_sample_ids():
    with open(os.path.join(outputs_dir, "samples.txt")) as f:
        return [s.rstrip() for s in f]


def print_metadata(all_sample_ids, output_filename):
    m = {}
    for sample_id in all_sample_ids:
        patient_id, treatment, time = sample_id.split('_')
        if patient_id not in m:
            m[patient_id] = [treatment, "", ""]
        if time == "T1":
            m[patient_id][1] = sample_id
        elif time == "T2":
            m[patient_id][2] = sample_id
        else:
            raise RuntimeError("Wrong time in sample name")
    with open(output_filename, "w") as f:
        print("patient_id",
              "treatment",
              "sample_id_before",
              "sample_id_after",
              sep=',',
              file=f)
        for patient_id in m:
            print(patient_id, *(m[patient_id]), sep=',', file=f)


def get_kraken_filename(sample_id):
    sample_subdir = os.path.join(outputs_dir, sample_id)
    kraken_dir = os.path.join(sample_subdir, "kraken2")
    if os.path.isdir(kraken_dir):
        for s in ["1", "2", sample_id]:
            kraken_filename = os.path.join(kraken_dir, s + ".kreport")
            if os.path.isfile(kraken_filename):
                return kraken_filename
        raise AssertionError("Can't find .kreport file in kraken2 subdirectory")
    return None


def run_bracken(kraken_filename, output_name, kreport_name, level):
    proc = subprocess.run(["../../Bracken/src/est_abundance.py",
                           "-i", kraken_filename,
                           "-k", bracken_db_path,
                           "-o", output_name,
                           "-t", "1",
                           "-l", level,
                           "--out-report", kreport_name])
    if proc.returncode != 0:
        raise RuntimeError("Bracken failed")


def create_bracken_report(sample_id, bracken_output, kreport_output):
    kraken_filename = get_kraken_filename(sample_id)
    if kraken_filename is None:
        print(sample_id, " missing kraken2, skipping sample")
        return

    output_name = os.path.join("outputs", "bracken_output", bracken_output)
    kreport_name = os.path.join("outputs", "bracken_output", kreport_output)
    run_bracken(kraken_filename, output_name, kreport_name, "G")


def parse_kreport(filename):
    level_map = {'R': 0, 'D': 1, 'K': 2, 'P': 3, 'C': 4, 'O': 5, 'F': 6, 'G': 7}

    def is_higher(s1, s2):
        if level_map[s1[0]] != level_map[s2[0]]:
            return level_map[s1[0]] < level_map[s2[0]]
        return s1[1:] < s2[1:]

    full_name = []
    result = []

    with open(filename) as f:
        for line in f:
            entries = line.rstrip().split('\t')
            level = entries[3].lstrip(). rstrip()

            while full_name and not is_higher(full_name[-1][0], level):
                full_name.pop()

            full_name.append((level, entries[5].lstrip().rstrip()))

            if level == "G":
                result.append(full_name + [entries[1]])
    return result


def get_amrplusplus_filename(sample_id):
    f = os.path.join(outputs_dir,
                     sample_id,
                     "amrplusplus",
                     "ResistomeResults",
                     "AMR_analytics_matrix.csv")
    assert(os.path.isfile(f))
    return f


def parse_amrplusplus(filename):
    result = []
    with open(filename) as f:
        first = True
        for line in f:
            if first:
                first = False
                continue
            entries = line.rstrip().split(',')
            count = entries[1].rstrip(".0")
            full_name = entries[0].split('|')
            if full_name[-1] == "RequiresSNPConfirmation":
                full_name.pop()
            full_name.append(full_name[0])
            full_name.pop(0)
            result.append(full_name + [count])
    return result


def get_read_count(sample_id):
    p = os.path.join(outputs_dir, sample_id, "preprocessing_read_counts.txt")
    with open(p) as f:
        first = True
        for line in f:
            if first:
                first = False
                continue
            return line.rstrip().split("\t")[3]


def create_json_output(sample_id, kreport_name, amrplusplus_name, output):
    kreport_result = parse_kreport(kreport_name)
    amrplusplus_result = parse_amrplusplus(amrplusplus_name)
    result = {"taxonomy": kreport_result,
              "resistome": amrplusplus_result,
              "read_count": get_read_count(sample_id)}
    with open(output, "w") as f:
        print(json.dumps(result, sort_keys=True, indent=4), file=f)


if __name__ == "__main__":
    all_sample_ids = get_all_sample_ids()

    os.mkdir("outputs/bracken_output")
    os.mkdir("outputs/samples")
    print_metadata(all_sample_ids, "outputs/metadata.csv")

    for id in all_sample_ids:
        create_bracken_report(id, id + ".bracken", id + ".kreport")
        create_json_output(id,
                           "outputs/bracken_output/" + id + ".kreport",
                           get_amrplusplus_filename(id),
                           "outputs/samples/" + id + ".json")

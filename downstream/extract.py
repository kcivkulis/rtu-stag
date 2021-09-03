import os
import shutil
import subprocess


outputs_dir = os.path.join(os.getenv("GROUP"), "projects_real", "ERAF184")
bracken_db_path = os.path.join(os.getenv("GROUP"), "databases", "full_ref_bafp", "database150mers.kmer_distrib")


def get_all_sample_ids():
    with open(os.path.join(outputs_dir, "samples.txt")) as f:
        return [s.rstrip() for s in f]


def print_metadata(output_filename):
    m = {}
    for sample_id in get_all_sample_ids():
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
        print("patient_id", "treatment", "sample_id_before", "sample_id_after", sep = ',', file = f)
        for patient_id in m:
            print(patient_id, *(m[patient_id]), sep = ',', file = f)


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


def get_amrplusplus_filename(sample_id):
    sample_subdir = os.path.join(outputs_dir, sample_id)
    amrplusplus_dir = os.path.join(sample_subdir, "amrplusplus", "RunResistome")
    if os.path.isdir(amrplusplus_dir):
        for f in os.listdir(amrplusplus_dir):
            if f.endswith(".mechanism.tsv"):
                return os.path.join(amrplusplus_dir, f)
        raise AssertionError("Can't find .mechanism.tsv file im amrplusplus subdirectory")
    return None


def create_amrplusplus_krona(input, output):
    with open(input) as i:
        with open(output, "w") as o:
            first = True
            for line in i:
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
                print(count + "\t" + "\t".join(full_name), file=o)


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


def create_bracken_reports():
    os.mkdir("outputs/bracken_output")
    for sample_id in get_all_sample_ids():
        kraken_filename = get_kraken_filename(sample_id)
        if kraken_filename == None:
            print(sample_id, " missing kraken2, skipping sample")
            continue

        output_name = os.path.join("outputs", "bracken_output", str(sample_id) + ".bracken")
        kreport_name = os.path.join("outputs", "bracken_output", str(sample_id) + ".kreport")
        run_bracken(kraken_filename, output_name, kreport_name, "G")


def extract_amrplusplus_reports():
    target_dir = os.path.join("outputs", "amrplusplus_report")
    os.mkdir(target_dir)
    for sample_id in get_all_sample_ids():
        amrplusplus_filename = get_amrplusplus_filename(sample_id)
        if amrplusplus_filename == None:
            print(sample_id, " missing amrplusplus, skipping sample")
            continue
        shutil.copyfile(amrplusplus_filename, os.path.join(target_dir, str(sample_id) + ".tsv"))
        amrplusplus_base = os.path.dirname(os.path.dirname(amrplusplus_filename))
        create_amrplusplus_krona(os.path.join(amrplusplus_base, "ResistomeResults", "AMR_analytics_matrix.csv"),
                                 os.path.join(target_dir, str(sample_id) + ".krona"))


def create_amrplusplus_krona_condensed():
    control_before, control_after, std2_before, std2_after, std3_before, std3_after = {}, {}, {}, {}, {}, {}
    c_control_before, c_control_after, c_std2_before, c_std2_after, c_std3_before, c_std3_after = 0, 0, 0, 0, 0, 0

    for filename in os.listdir("amrplusplus/"):
        input = "amrplusplus/" + filename
        sid = filename.rstrip(".csv")
        total = 0
        with open(input) as i:
            first = True
            for line in i:
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
                total += int(count)

        with open(input) as i:
            first = True
            for line in i:
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
                x = "\t".join(full_name)
                p = float(count) / float(total)

                def a(res, x, p):
                    if x in res:
                        res[x] += p
                    else:
                        res[x] = p

                if sid.endswith("Control_T1"):
                    a(control_before, x, p)
                elif sid.endswith("Control_T2"):
                    a(control_after, x, p)
                elif sid.endswith("STD2_T1"):
                    a(std2_before, x, p)
                elif sid.endswith("STD2_T2"):
                    a(std2_after, x, p)
                elif sid.endswith("STD3_T1"):
                    a(std3_before, x, p)
                elif sid.endswith("STD3_T2"):
                    a(std3_after, x, p)

        if sid.endswith("Control_T1"):
            c_control_before += 1
        elif sid.endswith("Control_T2"):
            c_control_after += 1
        elif sid.endswith("STD2_T1"):
            c_std2_before += 1
        elif sid.endswith("STD2_T2"):
            c_std2_after += 1
        elif sid.endswith("STD3_T1"):
            c_std3_before += 1
        elif sid.endswith("STD3_T2"):
            c_std3_after += 1

    for resistome, cnt, name in [(control_before, c_control_before, "control_before"),
                                 (control_after, c_control_after, "control_after"),
                                 (std2_before, c_std2_before, "std2_before"),
                                 (std2_after, c_std2_after, "std2_after"),
                                 (std3_before, c_std3_before, "std3_before"),
                                 (std3_after, c_std3_after, "std3_after")]:
        s = 0
        with open(name + ".krona", "w") as o:
            for v in resistome:
                print(str(resistome[v] / cnt) + "\t" + v, file=o)
                s += (resistome[v] / cnt)
        print(name, cnt, s)


if __name__ == "__main__":
    create_amrplusplus_krona_condensed()
    print_metadata("outputs/metadata.csv")
    create_bracken_reports()
    extract_amrplusplus_reports()

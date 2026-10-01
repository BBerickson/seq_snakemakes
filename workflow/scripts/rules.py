### common helper scripts ###

import os
import re

def _cutadapt_summary(input, output):
    import glob
    # Derive the directory from the first input file
    cutadapt_dir = os.path.dirname(input[0])
    # Glob all stats files in that directory
    all_files = sorted(glob.glob(os.path.join(cutadapt_dir, "*_cutadapt_stats.txt")))
    
    with open(output, "w") as out:
        for file in all_files:
            name = os.path.basename(file)
            name = re.sub("_cutadapt_stats.txt", "", name)

            for line in open(file, "r"):
                match = re.search("Pairs | pairs |with adapter", line)

                if match:
                  line = line.strip()
                  mett = re.split(r'\s{2,}', line)
                  num = mett[1].strip()
                  num = re.sub("\,", "", num)
                  num = re.sub(" ", "", num)
                  met = re.search("[\w\(\) ]+:", line).group(0)
                  met = re.sub(":", "", met)
                  met = re.sub(" ", "_", met)

                  out.write("%s\t%s\t%s\n" % (name, met, num))

def _clumpify_summary(input, output, pairedmap):
    if isinstance(pairedmap, str):
        import ast
        try:
            pairedmap = ast.literal_eval(pairedmap)
        except:
            print("Could not convert pairedmap from string to dict")
            return
    import glob
    # Derive the directory from the first input file
    clumpify_dir = os.path.dirname(input[0])
    # Glob all stats files in that directory
    all_files = sorted(glob.glob(os.path.join(clumpify_dir, "*_clumpify.log")))
    
    with open(output, "w") as out:
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_clumpify.log", "", name)
        
            for line in open(file, "r"):
                match = re.search("Reads In: |Duplicates Found:", line)
                if match:
                    num = re.search("[0-9,]{2,}", line)
                    num = int(re.sub(",", "", num.group(0)))
                    num_final = num // 2 if pairedmap[name] else num
                    met = re.search("[\w\(\) ]+:", line).group(0)
                    met = re.sub(":", "", met).strip().replace(" ", "_")
                    
                    out.write("%s\t%s\t%s\n" % (name, met, num_final))

def _bbduk_summary(input, output):
    import glob
    # Derive the directory from the first input file
    bbduk_dir = os.path.dirname(input[0])
    # Glob all stats files in that directory
    all_files = sorted(glob.glob(os.path.join(bbduk_dir, "*_bbduk.log")))
    
    with open(output, "w") as out:
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_bbduk.log", "", name)
            for line in open(file, "r"):
                if "Total Removed:" in line:
                    # Split by tab
                    mett = line.strip().split("\t")
                    
                    # Handle case where there might not be tabs (malformed log)
                    if len(mett) < 2:
                        # Try splitting by multiple spaces as fallback
                        mett = re.split(r'\s{2,}', line.strip())
                    
                    if len(mett) >= 2:
                        num = mett[1].strip()
                        num = re.sub(" reads ", "", num)
                        met = re.search(r"[\w\(\) ]+:", line).group(0)
                        met = re.sub(":", "", met)
                        met = met.strip()
                        met = re.sub(" ", "_", met)
                        
                        out.write("%s\t%s\t%s\n" % (name, met, num))
                    break  # Only process first matching line per file

def _bowtie_summary(input, output, index):
    with open(output, "w") as out:
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_" + index + "_bowtie_stats.txt", "", name)

            for line in open(file, "r"):
                match = re.search("overall alignment rate", line)

                if match:
                    num = re.search("[0-9.%]+", line)
                    num = num.group(0).strip()
                    met = re.search("[a-z\s]+", line).group(0)
                    met = met.strip()
                    met = re.sub(" ", "_", met)

                    out.write("%s\t%s\t%s\t%s\n" % (index, name, met, num))
                else:
                    line  = re.sub("; of these:", "", line.strip())
                    line  = re.sub(" \([0-9\.%]+\)", "", line)
                    words = line.split(" ")
                    num   = words[0]
                    met   = words[1:]
                    met   = " ".join(met)

                    out.write("%s\t%s\t%s\t%s\n" % (index, name, met, num))

def _hisatPE_summary(input, output, index):
    with open(output, "w") as out:
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_" + index + "_hisat_stats.txt", "", name)

            for line in open(file, "r"):
                match = re.search("overall alignment rate", line)

                if match:
                    num = re.search("[0-9.%]+", line)
                    num = num.group(0).strip()
                    met = re.search("[a-z\s]+", line).group(0)
                    met = met.strip()
                    met = re.sub(" ", "_", met)

                    out.write("%s\t%s\t%s\t%s\n" % (index, name, met, num))
                else:
                    line  = re.sub("; of these:", "", line.strip())
                    line  = re.sub(" \([0-9\.%]+\)", "", line)
                    words = line.split(" ")
                    num   = words[0]
                    met   = words[1:]
                    met   = " ".join(met)

                    out.write("%s\t%s\t%s\t%s\n" % (index, name, met, num))

def _dedup_summary(input, output):
    with open(output, "w") as out:
          metrics = [
            "Input Reads: [0-9]+",
            "Number of reads out: [0-9]+",
            "Total number of positions deduplicated: [0-9]+",
            "Mean number of unique UMIs per position: [0-9\.]+",
            "Max. number of unique UMIs per position: [0-9]+"
          ]

          for file in input:
              name  = os.path.basename(file)
              name  = re.sub("_dedup_stats.txt", "", name)

              for line in open(file, "r"):
                  for metric in metrics:
                      met = re.search(metric, line)

                      if met:
                          met = met.group(0)
                          num = re.search("[0-9\.]+$", met).group(0)
                          met = re.sub(": [0-9\.]+$", "", met)

                          out.write("%s\t%s\t%s\n" % (name, met, num))

                      
def _trim_galore_summary(input, output):
    import glob
    # Derive the directory from the first input file
    trimgalore_dir = os.path.dirname(input[0])
    # Glob all stats files in that directory
    all_files = sorted(glob.glob(os.path.join(trimgalore_dir, "*_trimming_report.txt")))

    with open(output, "w") as out:
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_trimming_report.txt", "", name)

            for line in open(file, "r"):
                # Trim Galore's report embeds a cutadapt-style summary block --
                # same wording cutadapt itself uses, just with Ribo-seq's
                # single-end phrasing ("Total reads processed", not "Total
                # read pairs processed").
                match = re.search(
                    r"Total reads processed|Reads with adapters|Reads written \(passing filters\)",
                    line,
                )
                if match:
                    line = line.strip()
                    mett = re.split(r'\s{2,}', line)
                    if len(mett) < 2:
                        continue
                    num = mett[1].strip()
                    num = re.sub(",", "", num)
                    num = re.sub(" ", "", num)
                    met = re.search(r"[\w\(\) ]+:", line).group(0)
                    met = re.sub(":", "", met)
                    met = re.sub(" ", "_", met.strip())

                    out.write("%s\t%s\t%s\n" % (name, met, num))


def _rseqc_summary(input, output):
    # read_distribution.py writes a fixed-format block of
    # "<group>  <Total_bases>  <Tag_count>  <Tags/Kb>" rows (CDS_Exons,
    # 5'UTR_Exons, 3'UTR_Exons, Introns, TSS_up_1kb, TES_down_1kb, etc.) --
    # pull out Tags/Kb per group, one row per sample.
    with open(output, "w") as out:
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_read_distribution.txt", "", name)

            in_table = False
            for line in open(file, "r"):
                if line.startswith("Group"):
                    in_table = True
                    continue
                if line.startswith("="):
                    if in_table:
                        break
                    continue
                if in_table:
                    fields = line.split()
                    if len(fields) == 4:
                        group, total_bases, tag_count, tags_per_kb = fields
                        out.write("%s\t%s\t%s\n" % (name, group, tags_per_kb))


def _sortmerna_summary(input, output):
    with open(output, "w") as out:
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_sortmerna.log", "", name)
            for line in open(file, "r"):
                m = re.search(r"Total reads passing E-value threshold[^:]*:\s*([0-9,]+)", line)
                if m:
                    out.write("%s\trRNA_reads\t%s\n" % (name, m.group(1).replace(",", "")))


def _strandedness_summary(input, output):
    import json
    with open(output, "w") as out:
        out.write("sample\texpected_format\tcompatible_fraction\n")
        for file in input:
            name = os.path.basename(file)
            name = re.sub("_lib_format_counts.json", "", name)
            if os.path.getsize(file) == 0:
                out.write("%s\tskipped (no SALMON_INDEX)\tNA\n" % name)
                continue
            with open(file) as f:
                counts = json.load(f)
            expected = counts.get("expected_format", "NA")
            frac = counts.get("num_compatible_fragments", 0)
            total = counts.get("num_assigned_fragments", 1) or 1
            out.write("%s\t%s\t%.3f\n" % (name, expected, frac / total))


def _extract_star_log_info(input_files, output_file, index):
    patterns = {
        "Number of input reads": ("reads", re.compile(r"Number of input reads \|\s+(.*)")),
        "Uniquely mapped reads number": ("aligned exactly 1 time", re.compile(r"Uniquely mapped reads number \|\s+(.*)")),
        "Number of reads mapped to multiple loci": ("aligned >1 times", re.compile(r"Number of reads mapped to multiple loci \|\s+(.*)")),
        "Uniquely mapped reads %": ("overall_alignment_rate", re.compile(r"Uniquely mapped reads % \|\s+(.*)%"))
    }

    with open(output_file, "w") as out:
        for file in input_files:
            name = os.path.basename(file).replace("_Log.final.out", "")
            with open(file, "r") as f:
                for line in f:
                    for key, (metric, pattern) in patterns.items():
                        match = pattern.search(line)
                        if match:
                            num = match.group(1)
                            if key == "Uniquely mapped reads %":
                                num += "%"
                            out.write("%s\t%s\t%s\t%s\n" % (index, name, metric, num))
                            break


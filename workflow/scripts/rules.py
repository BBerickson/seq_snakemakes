### common helper scripts ###

import os
import re

def _cutadapt_summary(input, output):
    with open(output, "w") as out:
        for file in input:
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


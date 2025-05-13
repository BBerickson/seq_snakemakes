### common helper scripts ###

import os
import re

def _extract_star_log_info(input_files, output_file):
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
                            out.write("{}\t{}\t{}\n".format(name, metric, num))
                            break


### common helper scripts ###

import os
import re
import gzip
import glob
import pandas as pd
import itertools
import subprocess

# estimates the amount of memory based on file size
def get_file_size_gb(filepath):
    """Returns the file size in gigabytes (GB)."""
    return os.path.getsize(filepath) / (1024 ** 3)

def memory_estimator(input_files, multiplier, minsize=1, maxsize=250, unit='GB'):
    """
    Estimates memory usage based on total file size.

    Parameters:
        input_files (list): List of file paths.
        multiplier (float): Scaling factor for estimated memory.
        minsize (int): Minimum memory in GB.
        maxsize (int): Maximum memory in GB.
        unit (str): 'GB' (default) or 'MB' for output unit.

    Returns:
        int: Estimated memory usage in the specified unit.
    """
    total_size_gb = sum(get_file_size_gb(f) for f in input_files)
    estimated_gb = total_size_gb * multiplier
    estimated_gb = max(estimated_gb, minsize)
    estimated_gb = min(estimated_gb, maxsize)

    if unit.upper() == 'MB':
        return round(estimated_gb * 1024)
    return round(estimated_gb)


# returns literal or ref for bbduk
def _bbduk_adapter_param(barcodes):
    if os.path.isfile(barcodes):
        return "ref={}".format(barcodes)
    else:
        if "," in barcodes:
            adapters = [a.strip() for a in barcodes.split(",") if a.strip()]
        else:
            adapters = [a.strip() for a in barcodes.split() if a.strip()]
        return "literal={}".format(','.join(adapters))

    
# gets number for norm fraction sample
def _get_norm_fraction(wildcards, index_type, filename):
    sample = wildcards.sample + "_" + index_type
    num = 0
    with open(filename, "r") as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) == 3 and parts[0] == sample:
                num = parts[2]
                break
    return float(num)
        
# Find all fastqs matching sample name in provided directories
def _find_fqs(sample, dirs):
    fq_pat   = ".*" + sample + r".+\.(fastq|fq)\.gz$"
    fq_paths = []
    
    for dir in dirs:
        all_files = glob.glob(os.path.abspath(os.path.join(dir, "*.gz")))
        paths     = [f for f in all_files if re.match(fq_pat, f)]
    
        for path in paths:
            fq_paths.append(path)

    if not fq_paths:
        sys.exit("ERROR: no fastqs found for " + fq_pat + ".")

    return fq_paths

# Determine the suffix (e.g. fastq.gz) for a list of fastqs matching a single
# sample name the expectation is that all fastqs in the list will have the same suffix
def _get_fq_sfx(fqs):
    sfx = []
    
    for fq in fqs:
        if re.search(r"\.fastq\.gz$", fq):
            sf = ".fastq.gz"
    
        if re.search(r"\.fq\.gz$", fq):
            sf = ".fq.gz"
    
        sfx.append(sf)
    
    sfx = set(sfx)
    
    if len(sfx) > 1:
        sys.exit("ERROR: Multiple fastqs found for a sample.")
    
    sfx = list(sfx)[0]

    return sfx

# For the fastq suffix (e.g. fastq.gz) get the complete suffix for both reads
# (e.g. _R1_001.fastq.gz)
def _get_full_fq_sfxs(sfx, paired=True):
    if sfx == ".fastq.gz":
        if paired:
            fq_sfx = ["_" + x + "_001" + sfx for x in ["R1", "R2"]]
        else:
            fq_sfx = ["_R1_001" + sfx]

    elif sfx == ".fq.gz":
        if paired:
            fq_sfx = ["_" + x + sfx for x in ["1", "2"]]
        else:
            fq_sfx = ["_1" + sfx]

    else:
        fq_sfx = [sfx]

    return fq_sfx

# Get the fastq file for both reads for the provided sample name
def _get_fqs(sample, dirs, link_dir, full_name = False, paired=True):

    fq_paths = _find_fqs(sample, dirs)
    
    sfx = _get_fq_sfx(fq_paths)

    sfxs = _get_full_fq_sfxs(sfx,paired)

    # Find matching fastqs and create symlinks 
    fastqs = []
    
    for full_sfx in sfxs:
        fq_pat = ".*/" + sample + ".*" + full_sfx
    
        fq = [f for f in fq_paths if re.match(fq_pat, f)]
    
        # Check for duplicate paths
        if not fq:
            sys.exit("ERROR: no fastqs found for " + fq_pat + ".")
    
        if len(fq) > 1:
            sys.exit("ERROR: Multiple fastqs found for " + fq_pat + ".")
    
        fq = fq[0]
    
        fastq = os.path.basename(fq)

        # Create symlinks
        # Using subprocess since os.symlink requires a target file instead of
        # target directory
        fq_lnk = link_dir + "/" + fastq
    
        if not os.path.exists(fq_lnk):
            cmd = "ln -s " + fq + " " + link_dir
    
            if cmd != "":
                subprocess.run(cmd, shell = True)
    
        # Return fastq path or just the fastq name
        if full_name:
            fastqs.append(fq_lnk)
    
        else:
            fastqs.append(fastq)
    
    return fastqs

# Find all bigwigs matching sample name in provided directories
def _find_bws(sample, dirs):
    bw_pat   = ".*" + sample + r".+\.bw$"
    bw_paths = []
    
    for dir in dirs:
        all_files = glob.glob(os.path.abspath(os.path.join(dir, "*.bw")))
        paths     = [f for f in all_files if re.match(bw_pat, f)]
        bw_paths.extend(paths)

    if not bw_paths:
        sys.exit("ERROR: no bigwigs found for " + sample + ".")

    return bw_paths

# Determine the suffix (e.g. fw/rev pos/neg) for a list of bigwigs matching a single
# sample name the expectation is that all bigwigs in the list will have the same suffix
def _get_bw_sfx(bws):
    sfx = set()

    for bw in bws:
        if re.search(r"_fw\.bw$", bw):
            sfx.add("_fw.bw")
        elif re.search(r"_pos\.bw$", bw):
            sfx.add("_pos.bw")
        else:
            sfx.add(".bw")
    
    return sfx.pop()

# For the bigwig suffix (e.g. .bw, fw/rev, pos/neg) get the complete suffix for both reads
def _get_full_bw_sfxs(sfx, paired=True):
    if sfx in ["_fw.bw", "_pos.bw"] and paired:
        return [f"_{x}.bw" for x in ["fw", "rev"]] if sfx == "_fw.bw" else [f"_{x}.bw" for x in ["pos", "neg"]]
    return [sfx]


# Get the fastq file for both reads for the provided sample name
def _get_bws(sample, dirs, link_dir, paired=True):
    bw_paths = _find_bws(sample, dirs)
    sfx = _get_bw_sfx(bw_paths)
    sfxs = _get_full_bw_sfxs(sfx,paired)

    # Find matching fastqs and create symlinks 
    bws = []
    for full_sfx in sfxs:
        bw_pat = ".*/" + sample + ".*" + full_sfx
        bw = [f for f in bw_paths if re.match(bw_pat, f)]
    
        # Check for duplicate paths
        if not bw:
            sys.exit("ERROR: no bigwigs found for " + sample + ".")
        if len(bw) > 1:
            sys.exit("ERROR: Multiple bigwigs found for " + sample + ": " + ", ".join(bw) + ".")
        bw = bw[0]
        bigwig = os.path.basename(bw)

        # Create symlinks
        # Using subprocess since os.symlink requires a target file instead of
        # target directory
        bw_lnk = link_dir + "/" + bigwig
    
        if not os.path.exists(bw_lnk):
            cmd = "ln -s " + bw + " " + link_dir
            if cmd != "":
                subprocess.run(cmd, shell = True)
        # Return bw path
        bws.append(bw_lnk)
    
    return bws
  
# build color sample name dictonary
def _get_colors(sample_key, color):
    if len(sample_key) >= len(color):
      color.extend(["0,0,0"]*len(set(sample_key)))
    res = {}
    for key in sample_key:
      for value in color:
        res[key] = value
        color.remove(value)
        break  
    return(res)

# grab color from dict
def _get_col(sample, cols_dict):
    return cols_dict.get(sample, "0,0,0")
        
# build bamCoverage scaleFactor sample dictonary
def _get_norm_scale(samples, sample_key, norm_type, index_sample):
    res = {}
    for key in sample_key:
        if norm_type.lower() in ["subsample", "none", "rpkm", "cpm", "bpm", "rpgc", "c"] or norm_type.isspace():
            res[key] = 1
        else:
            norm_type_with_index = norm_type + "_" + index_sample
            path = PROJ + "/counts/"
            sample = samples[key][0]
            filename = glob.glob(path + sample + "*_count.txt")
            if filename:
                filename = filename[0]
                with open(filename, "r") as f:
                    for line in f:
                        if norm_type_with_index in line:
                            num = line.strip().split()
                            if float(num[1]) != 0:
                                res[key] = 1000000 / float(num[1])
                            else:
                                res[key] = 1
                            break
                    else:
                        res[key] = "NA"
            else:
                res[key] = "NA"
    return res

# grab bowtie options for multimap cleaning 
def _extract_k_option(bowtie2, orientation): 
    # Extract the -k value using regex
    match = re.search(r'-k\s*(\d+)', bowtie2)
    if not match:
        return ""  # Return empty string if -k is not found

    # Normalize the output to "-k <value>"
    k_value = match.group(1)
    result = f"-k {k_value}"

    # Add --paired-end if orientation doesn't match R1 or R2
    if orientation not in {"R1", "R2"}:
        result += " --paired-end"
    
    return result



# grab normalization options for bamCoverage for each sample 
def _get_norm(newnam, samples, sample_key, norm_type, index_sample): 
    NORMS_DICT = _get_norm_scale(samples, sample_key, norm_type, index_sample)
    if norm_type in ["RPKM","CPM","BPM","RPGC"]:
      results = "--normalizeUsing " + norm_type
    else:
      results = "--normalizeUsing None"
    if newnam in NORMS_DICT:
      value = (NORMS_DICT[newnam])
      results = results + " --scaleFactor " + str(value)
    return results
    
# file nameing based on normalization and filter options    
def _get_normtype(normUsing, norm_type, blacklist, orientation):
    word = "_norm_" + norm_type
    
    match = re.search(r"--Offset\s+(-?\d+)", normUsing)
    if match:
      if orientation == "R2R1" or orientation == "R2":
        orent = -1
      else:
        orent = 1
      num = int(match.group(1)) * orent
      if num == -1:
        offset = "_3end"
      elif num == 1:
        offset = "_5end"
      else:
        offset = "_offset_" + str(num)
      word = word + offset
    if re.search(r"\S", blacklist):
      word = word + "_BL"
    return " ".join(word.split())
  
# set orientation for deeptools bamCovrage stranded data
def _get_bamCov_strand(mytype, orientation):
    if orientation == "R2R1" or orientation == "R2":
      return mytype
    elif mytype == "forward":
      return "reverse"
    else:
      return "forward"
    
# set featureCounts options based on orientation and .gtf/.saf ref file
def _get_featCout(orientation):
  # set options for read orientation
    if orientation == "R2R1" or orientation == "R2":
      results = "-s 2 "
    elif orientation == "R1R2" or orientation == "R1":
      results = "-s 1 "
    else:
      results = "-s 0 "
  # set paired end options
    if orientation != "R2" and orientation != "R1":
      results = results + "-p -C --countReadPairs "
  # Options for GTF file
    results = results + "-F GTF --extraAttributes 'gene_name,gene_biotype' -t gene -O "
    return results

# file nameing based on normalization, Matrix options, and genelist   
def _get_matrixtype(normUsing, computeMatrix,genelist):
    if genelist != "":
      genelist = "_" + genelist
    matchu = re.search(r"--upstream (\w+)", computeMatrix)
    if matchu:
      value = int(matchu.group(1))
      result = str(value / 1000) + "k_"
    else:
      result = "0k_"
    matchu5 = re.search(r"--unscaled5prime (\w+)", computeMatrix)
    if matchu5:
      value = int(matchu5.group(1))
      result = result + str(value / 1000) + "k_"
    matchb = re.search(r"--regionBodyLength (\w+)", computeMatrix)
    if matchb:
      value = int(matchb.group(1))
      result = result + str(value / 1000) + "k_"
    matchu3 = re.search(r"--unscaled3prime (\w+)", computeMatrix)
    if matchu3:
      value = int(matchu3.group(1))
      result = result + str(value / 1000) + "k_"
    matchd = re.search(r"--downstream (\w+)", computeMatrix)
    if matchd:
      value = int(matchd.group(1))
      result = result + str(value / 1000) + "k_"
    matchbin = re.search(r"--binSize (\w+)", computeMatrix)
    if matchbin:
      value = matchbin.group(1)
      result = result + value + "bin"
    else:
      result = result + "0bin"
    message = result + normUsing + genelist
    return message


# builds list of matrixtypes
def _get_all_matrixtypes(regions, normUsing, matrix_args, genelist):
  results = []
  for region in regions:
        if region == "PI" or region == "EI":
          computeMatrix = ""
        else:
          computeMatrix = matrix_args[f"region{region}"]  # Retrieve the specific CMD_PARAM
          
        result = _get_matrixtype(normUsing, computeMatrix, genelist)
        results.append(result)
    
  return results  # Return all results

# controls and sets output if subsample normalzation is set
def _get_bampath(bampath):
    if bampath == "subsample":
        word = "bams_sub"
    else:
        word = "bams"
    return word

## RGB conversion helpers
#
def _rgb2hex(samples,group,cols_dict):
    hex_colors = []
    for group in samples:
      if group in cols_dict:
        results = cols_dict[group]
      else:
        results = "0,0,0"
      myrgb = tuple(map(int, results.split(",")))
      myhex = "#{:02x}{:02x}{:02x}".format(*myrgb)
      hex_colors.append(myhex)
    # Join hex colors with space separator
    hex_colors = hex_colors + hex_colors
    return " ".join(hex_colors)
    
#    
def _rgb2hexplus(samples,group,cols_dict):
    hex_colors = []
    hex_colors2 = []
    for group in samples:
      if group in cols_dict:
        results = cols_dict[group]
      else:
        results = "0,0,0"
      myrgb = tuple(map(int, results.split(",")))
      myhex = "white,#{:02x}{:02x}{:02x}".format(*myrgb)
      myhex2 = "#{:02x}{:02x}{:02x},white".format(*myrgb)
      hex_colors.append(myhex)
      hex_colors2.append(myhex2)
    # Join hex colors with space separator
    hex_colors = hex_colors + hex_colors2
    return " ".join(hex_colors)

#
def _rgb2hexplus2(samples,group,cols_dict):
    hex_colors = []
    hex_colors2 = []
    for group in samples:
      if group in cols_dict:
        results = cols_dict[group]
      else:
        results = "0,0,0"
      myrgb = tuple(map(int, results.split(",")))
      myhex = "white,#{:02x}{:02x}{:02x}".format(*myrgb)
      myhex2 = "#{:02x}{:02x}{:02x},white".format(*myrgb)
      hex_colors.append(myhex)
      hex_colors2.append(myhex2)
    # Join hex colors with space separator
    hex_colors = hex_colors2 + hex_colors
    return " ".join(hex_colors)

#

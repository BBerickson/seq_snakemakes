# ===== Functions for pipeline =================================================

import os
import re


# Remove duplicated items from list while maintaining order of list
def _get_uniq_list(list_in):
    seen = set()
    out = []
    
    for item in list_in:
        if item not in seen:
            seen.add(item)
            out.append(item)

    return(out)


# Define persistent dict for subsampling
# load pytools within the function since the container will have this
# module but not the user
def _get_dict(dict_name, dict_dir):
    from pytools.persistent_dict import PersistentDict 
    d = PersistentDict(dict_name, container_dir = dict_dir)

    return d


# Clear persistent dict
def _clear_dict(dict_name, dict_dir):
    d = _get_dict(dict_name, dict_dir)

    d.clear()

    return d


# Define persistent dict for per gene subsampling
def _get_gene_sub_dict(group, region, dict_dir):
    dict_name = "GENE_SUB_DICT_" + group + "_" + region
    dict_dir  = dict_dir + "/" + dict_name

    os.makedirs(dict_dir, exist_ok = True)

    gene_sub_dict = _get_dict(dict_name, dict_dir)

    return gene_sub_dict


# Clear persistent dict for per gene subsampling
def _clear_gene_sub_dict(group, region, dict_dir):
    gene_sub_dict = _get_gene_sub_dict(group, region, dict_dir)

    gene_sub_dict.clear()

    return gene_sub_dict


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
# sample name the expectation is that all fastqs in the list will have the
# same suffix
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
        sys.exit("ERROR: Multiple fastqs found for " + sample + ".")
    
    sfx = list(sfx)[0]

    return sfx


# For the fastq suffix (e.g. fastq.gz) get the complete suffix for both reads
# (e.g. _R1_001.fastq.gz)
def _get_full_fq_sfxs(sfx):
    if sfx == ".fastq.gz":
        fq_sfx = ["_" + x + "_001" + sfx for x in ["R1", "R2"]]

    elif sfx == ".fq.gz":
        fq_sfx = ["_" + x + sfx for x in ["1", "2"]]

    else:
        fq_sfx = sfx

    return fq_sfx


# Get the fastq file for both reads for the provided sample name
def _get_fqs(sample, dirs, link_dir, full_name = False):

    fq_paths = _find_fqs(sample, dirs)
    
    sfx = _get_fq_sfx(fq_paths)

    sfxs = _get_full_fq_sfxs(sfx)

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


# Function to retrieve fastq paths
def _get_fq_paths(wildcards):
    fqs = _get_fqs(wildcards.sample, RAW_DATA, FASTQ_DIR, full_name = True)

    return fqs


# set options for cutadapt 
def _get_ca_options(orientation, umi_pattern):
  # Check conditions to determine Cutadapt options
  if re.search(r'N{6,}', umi_pattern):
      if orientation == "R1":
          ca_opts = "--pair-filter=first --discard-untrimmed "
      else:
          ca_opts = "--pair-filter=any --discard-untrimmed "
  else:
      ca_opts = ""

  return ca_opts


#!/usr/bin/env python

# Import modules
import subprocess as sp
import sys, re, copy, os, codecs, gzip
from argparse import ArgumentParser, Action as ArgParseAction
from collections import OrderedDict, defaultdict
from statistics import median as stat_median, mean as stat_mean
from typing import Optional
import pysam

# updated ggplot size to linewidth BE 2023-01-10
# 2024-12-18 added --filter_duplicates
# 2025-xx-xx modernized for Python 3.10+ and R 4+

__version__ = "1.2.0"


def get_version() -> str:
    """Return version information."""
    return f"ggsashimi v{__version__}"


def define_options() -> ArgumentParser:

    class DebugInfoAction(ArgParseAction):

        def __init__(self, option_strings, dest, **kwargs):
            super().__init__(option_strings, dest, nargs=0, **kwargs)

        def __call__(self, parser, namespace, values, option_string=None):
            returncode = 0
            try:
                get_debug_info()
            except sp.CalledProcessError as CPE:
                print(f"ERROR: {CPE.output.strip().decode('utf-8')}")
                returncode = CPE.returncode
            except Exception as e:
                print(f"ERROR: {e}")
                returncode = 1
            finally:
                parser.exit(returncode)

    parser = ArgumentParser(description="Create sashimi plot for a given genomic region")
    parser.add_argument("-b", "--bam", type=str, required=True,
        help="""
        Individual bam file or file with a list of bam files.
        In the case of a list of files the format is tsv:
        1col: id for bam file,
        2col: path of bam file,
        3+col: additional columns
        """)
    parser.add_argument("-c", "--coordinates", type=str, required=True,
        help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
    parser.add_argument("-o", "--out-prefix", type=str, dest="out_prefix", default="sashimi",
        help="Prefix for plot file name [default=%(default)s]")
    parser.add_argument("-S", "--out-strand", type=str, dest="out_strand", default="both",
        help="Only for --strand other than 'NONE'. Choose which signal strand to plot: <both> <plus> <minus> [default=%(default)s]")
    parser.add_argument("-M", "--min-coverage", type=int, default=1, dest="min_coverage",
        help="Minimum number of reads supporting a junction to be drawn [default=1]")
    parser.add_argument("-j", "--junctions-bed", type=str, dest="junctions_bed", default="",
        help="Junction BED file name [default=no junction file]")
    parser.add_argument("-g", "--gtf",
        help="Gtf file with annotation (only exons is enough)")
    parser.add_argument("--filter_duplicates", action="store_true",
        help="Filter out transcripts with identical coordinates")
    parser.add_argument("-s", "--strand", default="NONE", type=str,
        help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
    parser.add_argument("--shrink", action="store_true",
        help="Shrink the junctions by a factor for nicer display [default=%(default)s]")
    parser.add_argument("-O", "--overlay", type=int,
        help="Index of column with overlay levels (1-based)")
    parser.add_argument("-A", "--aggr", type=str, default="",
        help="""Aggregate function for overlay: <mean> <median> <mean_j> <median_j>.
                Use mean_j | median_j to keep density overlay but aggregate junction counts [default=no aggregation]""")
    parser.add_argument("-C", "--color-factor", type=int, dest="color_factor",
        help="Index of column with color levels (1-based)")
    parser.add_argument("--alpha", type=float, default=0.5,
        help="Transparency level for density histogram [default=%(default)s]")
    parser.add_argument("-P", "--palette", type=str,
        help="Color palette file. tsv file with >=1 columns, where the color is the first column. Both R color names and hexadecimal values are valid")
    parser.add_argument("-L", "--labels", type=int, dest="labels", default=1,
        help="Index of column with labels (1-based) [default=%(default)s]")
    parser.add_argument("--fix-y-scale", default=False, action="store_true", dest="fix_y_scale",
        help="Fix y-scale across individual signal plots [default=%(default)s]")
    parser.add_argument("--height", type=float, default=2,
        help="Height of the individual signal plot in inches [default=%(default)s]")
    parser.add_argument("--ann-height", type=float, default=1.5, dest="ann_height",
        help="Height of annotation plot in inches [default=%(default)s]")
    parser.add_argument("--width", type=float, default=10,
        help="Width of the plot in inches [default=%(default)s]")
    parser.add_argument("--base-size", type=float, default=14, dest="base_size",
        help="Base font size of the plot in pch [default=%(default)s]")
    parser.add_argument("-F", "--out-format", type=str, default="pdf", dest="out_format",
        help="Output file format: <pdf> <svg> <png> <jpeg> <tiff> [default=%(default)s]")
    parser.add_argument("-R", "--out-resolution", type=int, default=300, dest="out_resolution",
        help="Output file resolution in PPI (pixels per inch). Applies only to raster output formats [default=%(default)s]")
    parser.add_argument("--debug-info", action=DebugInfoAction,
        help="Show several system information useful for debugging purposes [default=%(default)s]")
    parser.add_argument("--version", action="version", version=get_version())
    return parser


def parse_coordinates(c: str) -> tuple[str, int, int]:
    c = c.replace(",", "")
    chr_ = c.split(":")[0]
    start, end = c.split(":")[1].split("-")
    # Convert to 0-based
    start, end = int(start) - 1, int(end)
    return chr_, start, end


def count_operator(
    CIGAR_op: str,
    CIGAR_len: int,
    pos: int,
    start: int,
    end: int,
    a: list,
    junctions: dict,
) -> int:

    # Match
    if CIGAR_op == "M":
        for i in range(pos, pos + CIGAR_len):
            if i < start or i >= end:
                continue
            ind = i - start
            a[ind] += 1

    # Insertion or Soft-clip
    if CIGAR_op in ("I", "S"):
        return pos

    # Deletion — advance position below
    # Junction
    if CIGAR_op == "N":
        don = pos
        acc = pos + CIGAR_len
        if don > start and acc < end:
            junctions[(don, acc)] = junctions.get((don, acc), 0) + 1

    pos = pos + CIGAR_len
    return pos


def flip_read(s: str, samflag: int) -> int:
    if s in ("NONE", "SENSE"):
        return 0
    if s == "ANTISENSE":
        return 1
    if s == "MATE1_SENSE":
        return 0 if int(samflag) & 64 else 1
    if s == "MATE2_SENSE":
        return 1 if int(samflag) & 64 else 0
    return 0


def read_bam(f: str, c: str, s: str) -> tuple[dict, dict]:
    chr_, start, end = parse_coordinates(c)

    a: dict[str, list] = {"+": [0] * (end - start)}
    junctions: dict[str, OrderedDict] = {"+": OrderedDict()}
    if s != "NONE":
        a["-"] = [0] * (end - start)
        junctions["-"] = OrderedDict()

    samfile = pysam.AlignmentFile(f)

    for read in samfile.fetch(chr_, start, end):
        if read.is_unmapped:
            continue

        samflag = read.flag
        read_start = read.reference_start + 1
        CIGAR = read.cigarstring

        if CIGAR is None:
            continue

        # Ignore reads with exotic CIGAR operators
        if any(op in CIGAR for op in ("H", "P", "X", "=")):
            continue

        read_strand = ["+", "-"][flip_read(s, samflag) ^ bool(int(samflag) & 16)]
        if s == "NONE":
            read_strand = "+"

        CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
        CIGAR_ops = re.split(r"\d+", CIGAR)[1:]

        pos = read_start
        for n, CIGAR_op in enumerate(CIGAR_ops):
            CIGAR_len = int(CIGAR_lens[n])
            pos = count_operator(
                CIGAR_op, CIGAR_len, pos, start, end,
                a[read_strand], junctions[read_strand],
            )

    samfile.close()
    return a, junctions


def get_bam_path(index: str, path: str) -> str:
    if os.path.isabs(path):
        return path
    base_dir = os.path.dirname(index)
    return os.path.join(base_dir, path)


def read_bam_input(f: str, overlay, color, label):
    if f.endswith(".bam"):
        bn = f.strip().split("/")[-1].removesuffix(".bam")
        yield bn, f, None, None, bn
        return
    with codecs.open(f, encoding="utf-8") as openf:
        for line in openf:
            line_sp = line.strip().split("\t")
            bam = get_bam_path(f, line_sp[1])
            overlay_level = line_sp[overlay - 1] if overlay else None
            color_level = line_sp[color - 1] if color else None
            label_text = line_sp[label - 1] if label else None
            yield line_sp[0], bam, overlay_level, color_level, label_text


def prepare_for_R(a: list, junctions: dict, c: str, m: int) -> tuple:
    _, start, _ = parse_coordinates(args.coordinates)

    x = list(i + start for i in range(len(a)))
    y = a

    dons, accs, yd, ya, counts = [], [], [], [], []

    for (don, acc), n in junctions.items():
        if n < m:
            continue
        dons.append(don)
        accs.append(acc)
        counts.append(n)
        yd.append(a[don - start - 1])
        ya.append(a[acc - start + 1])

    return x, y, dons, accs, yd, ya, counts


def intersect_introns(data):
    data = sorted(data)
    if not data:
        return
    it = iter(data)
    a, b = next(it)
    for c, d in it:
        if b > c:
            b = min(b, d)
            a = max(a, c)
        else:
            yield a, b
            a, b = c, d
    yield a, b


def shrink_density(x: list, y: list, introns) -> tuple[list, list]:
    new_x, new_y = [], []
    shift = 0
    start = 0
    for a, b in introns:
        end = x.index(a) + 1
        new_x += [int(i - shift) for i in x[start:end]]
        new_y += y[start:end]
        start = x.index(b)
        l = b - a
        shift += l - l ** 0.7
    new_x += [int(i - shift) for i in x[start:]]
    new_y += y[start:]
    return new_x, new_y


def shrink_junctions(dons: list, accs: list, introns) -> tuple[dict, list, list]:
    new_dons, new_accs = [0] * len(dons), [0] * len(accs)
    real_introns: dict = {}
    shift_acc = 0
    shift_don = 0
    s: set = set()
    junctions = list(zip(dons, accs))
    for a, b in introns:
        l = b - a
        shift_acc += l - int(l ** 0.7)
        real_introns[a - shift_don] = a
        real_introns[b - shift_acc] = b
        for i, (don, acc) in enumerate(junctions):
            if a >= don and b <= acc:
                if (don, acc) not in s:
                    new_dons[i] = don - shift_don
                    new_accs[i] = acc - shift_acc
                else:
                    new_accs[i] = acc - shift_acc
                s.add((don, acc))
        shift_don = shift_acc
    return real_introns, new_dons, new_accs


def read_palette(f: Optional[str]) -> tuple:
    palette = ("#ff0000", "#00ff00", "#0000ff", "#000000")
    if f:
        with open(f) as openf:
            palette = tuple(line.split("\t")[0].strip() for line in openf)
    return palette


def read_gtf(f: str, c: str, filter_duplicates: bool = False):
    exons: OrderedDict = OrderedDict()
    transcripts: OrderedDict = OrderedDict()
    genes: OrderedDict = OrderedDict()
    chr_, start, end = parse_coordinates(c)
    end = end - 1

    opener = gzip.open(f, "rt") if f.endswith(".gz") else open(f)
    with opener as openf:
        for line in openf:
            if line.startswith("#"):
                continue
            el_chr, _, el, el_start, el_end, _, strand, _, tags = line.strip().split("\t")
            if el_chr != chr_:
                continue
            if el not in ("transcript", "exon"):
                continue
            transcript_id_match = re.findall(r'transcript_id\s*("[^"]+")', tags)
            gene_name_match = re.findall(r'gene_name\s*("[^"]+")', tags)
            if not transcript_id_match or not gene_name_match:
                continue
            transcript_id = transcript_id_match[0]
            gene_name = gene_name_match[0]
            el_start, el_end = int(el_start) - 1, int(el_end)
            strand = '"' + strand + '"'

            if el == "transcript":
                if el_end > start and el_start < end:
                    transcripts[transcript_id] = max(start, el_start), min(end, el_end), strand
                    genes[transcript_id] = gene_name
                continue

            if el == "exon":
                if start < el_start < end or start < el_end < end:
                    exons.setdefault(transcript_id, []).append(
                        (max(el_start, start), min(end, el_end), strand)
                    )

    if filter_duplicates:
        filtered_transcripts: OrderedDict = OrderedDict()
        filtered_exons: defaultdict = defaultdict(list)
        filtered_genes: OrderedDict = OrderedDict()
        seen_coords: set = set()

        for transcript_id, coords in exons.items():
            coords_tuple = tuple(sorted(coords))
            if coords_tuple not in seen_coords:
                seen_coords.add(coords_tuple)
                filtered_transcripts[transcript_id] = transcripts[transcript_id]
                filtered_exons[transcript_id] = coords
                filtered_genes[transcript_id] = genes[transcript_id]

        return filtered_transcripts, filtered_exons, filtered_genes
    else:
        return transcripts, exons, genes


def make_introns(transcripts, exons, intersected_introns=None) -> dict:
    new_transcripts = copy.deepcopy(transcripts)
    new_exons = copy.deepcopy(exons)
    introns: OrderedDict = OrderedDict()

    if intersected_introns:
        for tx, (tx_start, tx_end, strand) in new_transcripts.items():
            total_shift = 0
            for a, b in intersected_introns:
                l = b - a
                shift = l - int(l ** 0.7)
                total_shift += shift
                for i, (exon_start, exon_end, _strand) in enumerate(exons.get(tx, [])):
                    new_exon_start, new_exon_end = new_exons[tx][i][:2]
                    if a < exon_start:
                        if b > exon_end:
                            if i == len(exons[tx]) - 1:
                                total_shift = total_shift - shift + (exon_start - a) * (1 - int(l ** -0.3))
                            shift = (exon_start - a) * (1 - int(l ** -0.3))
                            new_exon_end = new_exons[tx][i][1] - shift
                        new_exon_start = new_exons[tx][i][0] - shift
                    if b <= exon_end:
                        new_exon_end = new_exons[tx][i][1] - shift
                    new_exons[tx][i] = (new_exon_start, new_exon_end, _strand)
            tx_start = min(tx_start, sorted(new_exons.get(tx, [[sys.maxsize]]))[0][0])
            new_transcripts[tx] = (tx_start, tx_end - total_shift, strand)

    for tx, (tx_start, tx_end, strand) in new_transcripts.items():
        intron_start = tx_start
        ex_end = 0
        for ex_start, ex_end, strand in sorted(new_exons.get(tx, [])):
            intron_end = ex_start
            if tx_start < ex_start:
                introns.setdefault(tx, []).append((intron_start, intron_end, strand))
            intron_start = ex_end
        if tx_end > ex_end:
            introns.setdefault(tx, []).append((intron_start, tx_end, strand))

    return {"transcripts": new_transcripts, "exons": new_exons, "introns": introns}


def gtf_for_ggplot(annotation: dict, start: int, end: int, arrow_bins: int, genes: dict) -> str:
    arrow_space = int((end - start) / arrow_bins)
    s = """
    # data table with exons
    ann_list = list(
        "exons" = data.table(),
        "introns" = data.table()
    )
    gene_list = list(
        "exons" = data.table(),
        "introns" = data.table()
    )
    """

    if annotation["exons"]:
        s += """
        ann_list[['exons']] = data.table(
            tx = rep(c(%(tx_exons)s), c(%(n_exons)s)),
            start = c(%(exon_start)s),
            end = c(%(exon_end)s),
            strand = c(%(strand)s)
        )
        gene_list[['exons']] = data.table(
            tx = c(%(tx_gene2)s),
            gene = c(%(tx_gene)s)
        )
        ann_list[['exons']] <- left_join(ann_list[['exons']], gene_list[['exons']], by = "tx")
        ann_list[['exons']] <- dplyr::mutate(ann_list[['exons']], tx=paste(tx, gene, sep=":"))
        """ % {
            "tx_exons": ",".join(annotation["exons"].keys()),
            "tx_gene": ",".join(genes.values()),
            "tx_gene2": ",".join(genes.keys()),
            "n_exons": ",".join(map(str, map(len, annotation["exons"].values()))),
            "exon_start": ",".join(map(str, (v[0] for vs in annotation["exons"].values() for v in vs))),
            "exon_end": ",".join(map(str, (v[1] for vs in annotation["exons"].values() for v in vs))),
            "strand": ",".join(map(str, (v[2] for vs in annotation["exons"].values() for v in vs))),
        }

    if annotation["introns"]:
        s += """
        ann_list[['introns']] = data.table(
            tx = rep(c(%(tx_introns)s), c(%(n_introns)s)),
            start = c(%(intron_start)s),
            end = c(%(intron_end)s),
            strand = c(%(strand)s)
        )
        gene_list[['introns']] <- data.table(
            tx = c(%(tx_gene2)s),
            gene = c(%(tx_gene)s)
        )
        ann_list[['introns']] <- left_join(ann_list[['introns']], gene_list[['introns']], by = "tx")
        ann_list[['introns']] <- dplyr::mutate(ann_list[['introns']], tx=paste(tx, gene, sep=":"))
        # Create data table for strand arrows
        txarrows = data.table()
        introns = ann_list[['introns']]
        # Add right-pointing arrows for plus strand
        if ("+" %%in%% introns$strand && nrow(introns[strand=="+" & end-start>5, ]) > 0) {
            txarrows = rbind(
                txarrows,
                introns[strand=="+" & end-start>5, list(
                    seq(start+4, end, by=%(arrow_space)s)-1,
                    seq(start+4, end, by=%(arrow_space)s)
                    ), by=.(tx, start, end)
                ]
            )
        }
        # Add left-pointing arrows for minus strand
        if ("-" %%in%% introns$strand && nrow(introns[strand=="-" & end-start>5, ]) > 0) {
            txarrows = rbind(
                txarrows,
                introns[strand=="-" & end-start>5, list(
                    seq(start, max(start+1, end-4), by=%(arrow_space)s),
                    seq(start, max(start+1, end-4), by=%(arrow_space)s)-1
                    ), by=.(tx, start, end)
                ]
            )
        }
        """ % {
            "tx_introns": ",".join(annotation["introns"].keys()),
            "tx_gene": ",".join(genes.values()),
            "tx_gene2": ",".join(genes.keys()),
            "n_introns": ",".join(map(str, map(len, annotation["introns"].values()))),
            "intron_start": ",".join(map(str, (v[0] for vs in annotation["introns"].values() for v in vs))),
            "intron_end": ",".join(map(str, (v[1] for vs in annotation["introns"].values() for v in vs))),
            "strand": ",".join(map(str, (v[2] for vs in annotation["introns"].values() for v in vs))),
            "arrow_space": arrow_space,
        }

    s += """
    gtfp = ggplot()
    if (length(ann_list[['introns']]) > 0) {
        gtfp = gtfp + geom_segment(data=ann_list[['introns']], aes(x=start, xend=end, y=tx, yend=tx), linewidth=0.3)
        gtfp = gtfp + geom_segment(data=txarrows, aes(x=V1, xend=V2, y=tx, yend=tx), arrow=arrow(length=unit(0.02,"npc")))
    }
    if (length(ann_list[['exons']]) > 0) {
        gtfp = gtfp + geom_segment(data=ann_list[['exons']], aes(x=start, xend=end, y=tx, yend=tx), linewidth=5, alpha=1)
    }
    gtfp = gtfp + scale_y_discrete(expand=c(0, 0.5))
    gtfp = gtfp + scale_x_continuous(expand=c(0, 0.25))
    gtfp = gtfp + coord_cartesian(xlim = c(%s, %s))
    gtfp = gtfp + labs(y=NULL)
    gtfp = gtfp + theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank())
    """ % (start, end)

    return s


def setup_R_script(h: float, w: float, b: float, label_dict: dict) -> str:
    s = """
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(dplyr))
    suppressPackageStartupMessages(library(grid))
    suppressPackageStartupMessages(library(gridExtra))
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(gtable))
    suppressPackageStartupMessages(library(labeling))

    scale_lwd = function(r) {
        lmin = 0.1
        lmax = 4
        return( r*(lmax-lmin)+lmin )
    }

    base_size = %(b)s
    height = ( %(h)s + base_size*0.352777778/67 ) * 1.02
    width = %(w)s
    theme_set(theme_bw(base_size=base_size))
    theme_update(
        plot.margin = unit(c(15,15,15,15), "pt"),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(linewidth=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(angle=0, vjust=0.5)
    )

    labels = list(%(labels)s)

    density_list = list()
    junction_list = list()

    """ % {
        "h": h,
        "w": w,
        "b": b,
        "labels": ",".join(f'"{id}"="{lab}"' for id, lab in label_dict.items()),
    }
    return s


def make_R_lists(
    id_list: list,
    d: dict,
    overlay_dict: dict,
    aggr: str,
    intersected_introns,
) -> str:
    s = ""
    aggr_f = {
        "mean": stat_mean,
        "median": stat_median,
    }
    id_list = id_list if not overlay_dict else list(overlay_dict.keys())
    shrinked_introns: dict = {}

    for k in id_list:
        shrinked_introns_k: dict = {}
        shrinked_intronsid: dict = {}
        x, y, dons, accs, yd, ya, counts = [], [], [], [], [], [], []

        if not overlay_dict:
            x, y, dons, accs, yd, ya, counts = d[k]
            if intersected_introns:
                x, y = shrink_density(x, y, intersected_introns)
                shrinked_introns_k, dons, accs = shrink_junctions(dons, accs, intersected_introns)
                shrinked_introns.update(shrinked_introns_k)
        else:
            for id_ in overlay_dict[k]:
                xid, yid, donsid, accsid, ydid, yaid, countsid = d[id_]
                if intersected_introns:
                    xid, yid = shrink_density(xid, yid, intersected_introns)
                    shrinked_intronsid, donsid, accsid = shrink_junctions(donsid, accsid, intersected_introns)
                    shrinked_introns.update(shrinked_intronsid)
                x += xid
                y += yid
                dons += donsid
                accs += accsid
                yd += ydid
                ya += yaid
                counts += countsid
            if aggr and "_j" not in aggr:
                x = d[overlay_dict[k][0]][0]
                y = list(map(aggr_f[aggr], zip(*(d[id_][1] for id_ in overlay_dict[k]))))
                if intersected_introns:
                    x, y = shrink_density(x, y, intersected_introns)

        s += """
        density_list[["%(id)s"]] = data.frame(x=c(%(x)s), y=c(%(y)s))
        junction_list[["%(id)s"]] = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))
        """ % {
            "id": k,
            "x": ",".join(map(str, x)),
            "y": ",".join(map(str, y)),
            "dons": ",".join(map(str, dons)),
            "accs": ",".join(map(str, accs)),
            "yd": ",".join(map(str, yd)),
            "ya": ",".join(map(str, ya)),
            "counts": ",".join(map(str, counts)),
        }

    if intersected_introns:
        s += """
        coord_dict = data.frame(shrinked=c(%(shrinked_introns_keys)s), real=c(%(shrinked_introns_values)s))
        intersected_introns = data.frame(real_x=c(%(intersected_introns_x)s), real_xend=c(%(intersected_introns_xend)s))
        """ % {
            "shrinked_introns_keys": ",".join(map(str, shrinked_introns.keys())),
            "shrinked_introns_values": ",".join(map(str, shrinked_introns.values())),
            "intersected_introns_x": ",".join(str(coord[0]) for coord in intersected_introns),
            "intersected_introns_xend": ",".join(str(coord[1]) for coord in intersected_introns),
        }
    return s


def plot(R_script: str) -> None:
    p = sp.Popen("R --vanilla --slave", shell=True, stdin=sp.PIPE)
    p.communicate(input=R_script.encode("utf-8"))


def colorize(d: dict, p: tuple, color_factor) -> str:
    levels = list(OrderedDict.fromkeys(d.values()).keys())
    n = len(levels)
    if n > len(p):
        p = (p * n)[:n]
    if color_factor:
        s = "color_list = list(%s)\n" % (
            ",".join(f'"{k}"="{p[levels.index(v)]}"' for k, v in d.items())
        )
    else:
        s = "color_list = list(%s)\n" % (
            ",".join(f'"{k}"="grey"' for k in d)
        )
    return s


def get_debug_info() -> None:
    import platform
    system = platform.system()
    info: OrderedDict = OrderedDict()
    info["OS"] = f"{system}-{platform.machine()}"
    if system == "Linux":
        release = sp.check_output(["lsb_release", "-ds"])
        info["Distro"] = release.strip().decode("utf-8")
    info["Python"] = platform.python_version()
    print(get_version())
    print("")
    maxlen = max(map(len, info.keys()))
    for k, v in info.items():
        print(f"{k:{maxlen}}: {v:>}")
    print("")

    r_exec = ";".join([
        "library(ggplot2)",
        "library(dplyr)",
        "library(grid)",
        "library(gridExtra)",
        "library(data.table)",
        "library(gtable)",
        "sessionInfo()",
    ])
    r_command = f"R --vanilla --slave -e '{r_exec}'"
    r_info = sp.check_output(r_command, shell=True, stderr=sp.STDOUT)
    print(r_info.strip().decode("utf-8"))


def r_template_fill(template: str, mapping: dict) -> str:
    """Fill an R script template using <<key>> placeholders.
    This avoids any interaction with Python's %-formatting, so R operators
    like %%, %in%, and %% 2 are passed through untouched.
    """
    result = template
    for key, val in mapping.items():
        result = result.replace(f"<<{key}>>", str(val))
    return result


# ---------------------------------------------------------------------------
# R script body (R 4+ compatible)
# Uses <<key>> placeholders instead of %(key)s to avoid conflicts with
# R's % operators (%%in%%, %% 2, etc.)
# ---------------------------------------------------------------------------

FIX_Y_SCALE_BLOCK = r"""
maxheight   = max(unlist(lapply(density_list, function(df){ max(df$y) })))
breaks_y    = labeling::extended(0, maxheight, m=4)
maxheight_j = 0
for (bam_index in 1:length(density_list)) {
    id        = names(density_list)[bam_index]
    d         = data.table(density_list[[id]])
    junctions = data.table(junction_list[[id]])
    row_i = c()
    if (nrow(junctions) > 0) {
        junctions$jlabel = as.character(junctions$count)
        junctions = setNames(
            junctions[, .(max(y), max(yend),
                          round(<<aggr>>(count)),
                          paste(jlabel, collapse=",")),
                      keyby=.(x, xend)],
            names(junctions)
        )
        if ("<<aggr>>" != "") {
            junctions = setNames(
                junctions[, .(max(y), max(yend),
                              round(<<aggr>>(count)),
                              round(<<aggr>>(count))),
                          keyby=.(x, xend)],
                names(junctions)
            )
        }
        junctions = distinct(junctions, jlabel, .keep_all=TRUE)
        junctions = arrange(junctions, y, yend)
        row_i = 1:nrow(junctions)
    }
    for (i in row_i) {
        j = as.numeric(junctions[i, 1:5])
        if ("<<aggr>>" != "") {
            j[3] = ifelse(length(d[x==j[1]-1, y])==0, 0,
                          max(as.numeric(d[x==j[1]-1, y])))
            j[4] = ifelse(length(d[x==j[2]+1, y])==0, 0,
                          max(as.numeric(d[x==j[2]+1, y])))
        }
        if (i %% 2 != 0) {
            set.seed(mean(j[3:4]))
            maxheight_j = max(maxheight_j, max(j[3:4]) * runif(1, 1.5, 2.8))
        }
    }
}
"""

COORD_DICT_BLOCK = r"""
if (exists('coord_dict')) {
    all_pos_shrinked  = do.call(rbind, density_list)$x
    s2r = merge(intersected_introns, coord_dict, by.x='real_xend', by.y='real')
    s2r = merge(s2r, coord_dict, by.x='real_x', by.y='real', suffixes=c('_xend','_x'))
    breaks_x_shrinked = labeling::extended(min(all_pos_shrinked), max(all_pos_shrinked), m=5)
    breaks_x  = c()
    out_range = c()
    for (b in breaks_x_shrinked) {
        iintron = FALSE
        for (j in 1:nrow(s2r)) {
            l = s2r[j, ]
            if (b >= l$shrinked_x && b <= l$shrinked_xend) {
                p    = (b - l$shrinked_x) / (l$shrinked_xend - l$shrinked_x)
                realb = round(l$real_x + p * (l$real_xend - l$real_x))
                breaks_x = c(breaks_x, realb)
                iintron  = TRUE
                break
            }
        }
        if (!iintron) {
            if (b <= min(s2r$shrinked_x)) {
                l = s2r[which.min(s2r$shrinked_x), ]
                if (any(b == all_pos_shrinked)) {
                    s     = l$shrinked_x - b
                    realb = l$real_x - s
                    breaks_x = c(breaks_x, realb)
                } else {
                    out_range <- c(out_range, which(breaks_x_shrinked == b))
                }
            } else if (b >= max(s2r$shrinked_xend)) {
                l = s2r[which.max(s2r$shrinked_xend), ]
                if (any(b == all_pos_shrinked)) {
                    s     = b - l$shrinked_xend
                    realb = l$real_xend + s
                    breaks_x = c(breaks_x, realb)
                } else {
                    out_range <- c(out_range, which(breaks_x_shrinked == b))
                }
            } else {
                delta       = b - s2r$shrinked_xend
                delta[delta < 0] = Inf
                l     = s2r[which.min(delta), ]
                s     = b - l$shrinked_xend
                realb = l$real_xend + s
                breaks_x = c(breaks_x, realb)
            }
        }
    }
    if (length(out_range)) {
        breaks_x_shrinked = breaks_x_shrinked[-out_range]
    }
}
"""

R_BODY_TEMPLATE = r"""
pdf(NULL)   # suppress the blank PDF that ggplotGrob can produce

<<fix_y_scale_block>>

<<coord_dict_block>>

density_grobs = list()

for (bam_index in 1:length(density_list)) {

    id   = names(density_list)[bam_index]
    d    = data.table(density_list[[id]])
    junctions = data.table(junction_list[[id]])

    # ---- Density bar plot ----
    gp = ggplot(d) +
        geom_bar(aes(x, y), width=1, position='identity', stat='identity',
                 fill=color_list[[id]], alpha=<<alpha>>)

    if (bam_index == 1) {
        gp = gp + ggtitle("<<coordinates>>")
    }
    gp = gp + labs(y=labels[[id]])

    if (exists('coord_dict')) {
        gp = gp + scale_x_continuous(expand=c(0, 0.25),
                                     breaks=breaks_x_shrinked, labels=breaks_x)
    } else {
        gp = gp + scale_x_continuous(expand=c(0, 0.25))
    }

    if (!<<fix_y_scale>>) {
        maxheight  = max(d[['y']])
        breaks_y   = labeling::extended(0, maxheight, m=4)
        min_ymid_bottom = 0
        if (nrow(junctions) > 0) {
            min_ymid_bottom = -1.2 * maxheight
        }
        gp = gp + scale_y_continuous(breaks=breaks_y,
                                     limits=c(min_ymid_bottom, NA))
    } else {
        min_ymid_bottom = 0
        if (nrow(junctions) > 0) {
            min_ymid_bottom = -1.2 * maxheight
        }
        gp = gp + scale_y_continuous(breaks=breaks_y,
                                     limits=c(min_ymid_bottom, max(maxheight, maxheight_j)))
    }

    # ---- Junction arcs ----
    row_i = c()
    if (nrow(junctions) > 0) {
        junctions$jlabel = as.character(junctions$count)
        junctions = setNames(
            junctions[, .(max(y), max(yend),
                          round(<<aggr>>(count)),
                          paste(jlabel, collapse=",")),
                      keyby=.(x, xend)],
            names(junctions)
        )
        if ("<<aggr>>" != "") {
            junctions = setNames(
                junctions[, .(max(y), max(yend),
                              round(<<aggr>>(count)),
                              round(<<aggr>>(count))),
                          keyby=.(x, xend)],
                names(junctions)
            )
        }
        junctions = distinct(junctions, x, xend, y, yend, jlabel, .keep_all=TRUE)
        junctions = arrange(junctions, x, xend)
        row_i = 1:nrow(junctions)
    }

    for (i in row_i) {
        j_tot_counts = sum(junctions[['count']])
        j = as.numeric(junctions[i, 1:5])

        if ("<<aggr>>" != "") {
            j[3] = ifelse(length(d[x==j[1]-1, y])==0, 0, max(as.numeric(d[x==j[1]-1, y])))
            j[4] = ifelse(length(d[x==j[2]+1, y])==0, 0, max(as.numeric(d[x==j[2]+1, y])))
        }

        xmid = round(mean(j[1:2]), 1)
        set.seed(mean(j[3:4]))
        ymid = max(maxheight - max(j[3:4]), max(j[3:4])) * runif(1, 1.3, 1.8)

        lwd       = scale_lwd(j[5] / j_tot_counts) * 2.5
        curve_par = gpar(lwd=lwd, col=color_list[[id]])

        nss = i
        if (nss %% 2 == 0) {   # bottom arc
            ymid = -runif(1, 0.5, 1.2) * maxheight
            curve = xsplineGrob(x=c(0,0,1,1), y=c(1,0,0,0), shape=1, gp=curve_par)
            gp = gp + annotation_custom(curve, j[1], xmid, 0, ymid)
            curve = xsplineGrob(x=c(1,1,0,0), y=c(1,0,0,0), shape=1, gp=curve_par)
            gp = gp + annotation_custom(curve, xmid, j[2], 0, ymid)
        } else {                # top arc
            curve = xsplineGrob(x=c(0,0,1,1), y=c(0,1,1,1), shape=1, gp=curve_par)
            gp = gp + annotation_custom(curve, j[1], xmid, j[3], ymid)
            curve = xsplineGrob(x=c(1,1,0,0), y=c(0,1,1,1), shape=1, gp=curve_par)
            gp = gp + annotation_custom(curve, xmid, j[2], j[4], ymid)
        }

        gp = gp + annotate("label", x=xmid, y=ymid,
                           label=as.character(junctions[i, 6]),
                           vjust=0.5, hjust=0.5,
                           label.padding=unit(0.01, "lines"),
                           label.size=NA, size=(base_size * 0.35))
    }

    # ---- Grob assembly ----
    gpGrob = ggplotGrob(gp)
    gpGrob$layout$clip[gpGrob$layout$name == "panel"] <- "off"

    if (bam_index == 1) {
        # ggplot2 >= 3.4 stable grob positions (no vs offset needed)
        maxWidth     = gpGrob$widths[2] + gpGrob$widths[3]
        maxYtextWidth = gpGrob$widths[3]
        xaxisGrob    = gtable_filter(gpGrob, "axis-b", trim=FALSE)
        xaxisGrob$heights[8] = gpGrob$heights[1]
        x.axis.height = gpGrob$heights[7] + gpGrob$heights[1]
    }

    # Remove x axis from density plots (added once after the loop)
    kept_names = gpGrob$layout$name[gpGrob$layout$name != "axis-b"]
    gpGrob = gtable_filter(gpGrob, paste(kept_names, sep="", collapse="|"), trim=FALSE)

    maxWidth      = grid::unit.pmax(maxWidth,      gpGrob$widths[2] + gpGrob$widths[3])
    maxYtextWidth = grid::unit.pmax(maxYtextWidth, gpGrob$widths[3])
    density_grobs[[id]] = gpGrob
}

# Append x-axis grob, then optional annotation grob
density_grobs[["xaxis"]] = xaxisGrob

if (<<has_gtf>> == 1) {
    gtfGrob = ggplotGrob(gtfp)
    maxWidth = grid::unit.pmax(maxWidth, gtfGrob$widths[2] + gtfGrob$widths[3])
    density_grobs[['gtf']] = gtfGrob
}

# Align all grob widths
for (id in names(density_grobs)) {
    density_grobs[[id]]$widths[1] <- (
        density_grobs[[id]]$widths[1] +
        maxWidth -
        (density_grobs[[id]]$widths[2] + maxYtextWidth)
    )
    density_grobs[[id]]$widths[3] <- maxYtextWidth
}

heights = unit.c(
    unit(rep(<<signal_height>>, length(density_list)), "in"),
    x.axis.height,
    unit(<<ann_height>> * <<has_gtf>>, "in")
)

argrobs = arrangeGrob(grobs=density_grobs, ncol=1, heights=heights)

if ("<<out_format>>" == "tiff") {
    ggsave("<<out>>", plot=argrobs, device="tiff",
           width=width, height=height, units="in",
           dpi=<<out_resolution>>, compression="lzw", limitsize=FALSE)
} else {
    ggsave("<<out>>", plot=argrobs, device="<<out_format>>",
           width=width, height=height, units="in",
           dpi=<<out_resolution>>, limitsize=FALSE)
}

dev.log = dev.off()
"""


if __name__ == "__main__":

    strand_dict  = {"plus": "+", "minus": "-"}
    strand_dict2 = {"+": "plus", "-": "minus"}

    parser = define_options()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    if args.aggr and not args.overlay:
        print("ERROR: Cannot apply aggregate function if overlay is not selected.")
        sys.exit(1)

    palette = read_palette(args.palette)

    bam_dict: dict    = {"+": OrderedDict()}
    overlay_dict: dict = OrderedDict()
    color_dict: dict  = OrderedDict()
    id_list: list     = []
    label_dict: dict  = OrderedDict()

    if args.strand != "NONE":
        bam_dict["-"] = OrderedDict()

    junctions_list: list = []
    if args.junctions_bed:
        junctions_list = []

    for id_, bam, overlay_level, color_level, label_text in read_bam_input(
        args.bam, args.overlay, args.color_factor, args.labels
    ):
        if not os.path.isfile(bam):
            continue
        a, junctions = read_bam(bam, args.coordinates, args.strand)

        # Python 3-safe zero-coverage check
        if list(a.keys()) == ["+"] and all(v == 0 for v in a["+"]):
            print(f"WARN: Sample {id_} has no reads in the specified area.")
            continue

        id_list.append(id_)
        label_dict[id_] = label_text

        for strand in a:
            if (
                args.strand == "NONE"
                or args.out_strand == "both"
                or strand == strand_dict[args.out_strand]
            ):
                if args.junctions_bed:
                    for (k0, k1), v in junctions[strand].items():
                        if v >= args.min_coverage:
                            junctions_list.append(
                                "\t".join([
                                    args.coordinates.split(":")[0],
                                    str(k0), str(k1), id_, str(v), strand,
                                ])
                            )
            bam_dict[strand][id_] = prepare_for_R(
                a[strand], junctions[strand], args.coordinates, args.min_coverage
            )

        if color_level is None:
            color_dict.setdefault(id_, id_)
        if overlay_level is not None:
            overlay_dict.setdefault(overlay_level, []).append(id_)
            label_dict[overlay_level] = overlay_level
            color_dict.setdefault(overlay_level, overlay_level)
        if overlay_level is None:
            color_dict.setdefault(id_, color_level)

    if not bam_dict["+"]:
        print("ERROR: No available bam files.")
        sys.exit(1)

    if args.junctions_bed:
        jbed_path = args.junctions_bed if args.junctions_bed.endswith(".bed") else args.junctions_bed + ".bed"
        with open(jbed_path, "w") as jbed:
            jbed.write("\n".join(sorted(junctions_list)))

    if args.gtf:
        transcripts, exons, genes = read_gtf(
            args.gtf, args.coordinates, filter_duplicates=args.filter_duplicates
        )

    if args.out_format not in ("pdf", "png", "svg", "tiff", "jpeg"):
        print(
            f"ERROR: Provided output format '{args.out_format}' is not available. "
            "Please select among 'pdf', 'png', 'svg', 'tiff' or 'jpeg'"
        )
        sys.exit(1)

    for strand in bam_dict:

        # Output file name
        if args.out_prefix.endswith((".pdf", ".png", ".svg", ".tiff", ".tif", ".jpeg", ".jpg")):
            out_split = os.path.splitext(args.out_prefix)
            if (
                args.out_format == out_split[1][1:]
                or args.out_format == "tiff" and out_split[1] in (".tiff", ".tif")
                or args.out_format == "jpeg" and out_split[1] in (".jpeg", ".jpg")
            ):
                args.out_prefix = out_split[0]
                out_suffix = out_split[1][1:]
            else:
                out_suffix = args.out_format
        else:
            out_suffix = args.out_format

        out_prefix = args.out_prefix
        if args.strand != "NONE":
            if args.out_strand != "both" and strand != strand_dict[args.out_strand]:
                continue
            out_prefix = f"{args.out_prefix}_{strand_dict2[strand]}"

        # Shrink introns if requested
        intersected_introns = None
        if args.shrink:
            introns = (v for vs in bam_dict[strand].values() for v in zip(vs[2], vs[3]))
            intersected_introns = list(intersect_introns(introns))

        # Plot height
        bam_height = args.height * len(id_list)
        if args.overlay:
            bam_height = args.height * len(overlay_dict)
        if args.gtf:
            bam_height += args.ann_height

        # Build R script
        R_script = setup_R_script(bam_height, args.width, args.base_size, label_dict)
        R_script += colorize(color_dict, palette, args.color_factor)

        arrow_bins = 50
        if args.gtf:
            annotation = make_introns(transcripts, exons, intersected_introns)
            x = list(bam_dict[strand].values())[0][0]
            if args.shrink:
                x, _ = shrink_density(x, x, intersected_introns)
            R_script += gtf_for_ggplot(annotation, x[0], x[-1], arrow_bins, genes)

        R_script += make_R_lists(id_list, bam_dict[strand], overlay_dict, args.aggr, intersected_introns)

        aggr_r = args.aggr.rstrip("_j") or "mean"  # fallback so R block is valid syntax

        # Fill FIX_Y_SCALE_BLOCK with its own <<aggr>> placeholder before embedding
        fix_y_block = r_template_fill(FIX_Y_SCALE_BLOCK, {"aggr": aggr_r}) if args.fix_y_scale else ""

        R_script += r_template_fill(R_BODY_TEMPLATE, {
            "fix_y_scale_block": fix_y_block,
            "coord_dict_block":  COORD_DICT_BLOCK,
            "alpha":             args.alpha,
            "coordinates":       args.coordinates,
            "aggr":              aggr_r,
            "fix_y_scale":       "TRUE" if args.fix_y_scale else "FALSE",
            "has_gtf":           float(bool(args.gtf)),
            "signal_height":     args.height,
            "ann_height":        args.ann_height,
            "out":               f"{out_prefix}.{out_suffix}",
            "out_format":        args.out_format,
            "out_resolution":    args.out_resolution,
        })

        if os.getenv("GGSASHIMI_DEBUG") is not None:
            with open("R_script", "w") as r:
                r.write(R_script)
        else:
            plot(R_script)

    sys.exit(0)

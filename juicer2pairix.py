#!/usr/bin/env python
import argparse
import pysam
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description='convert juicer format to pairix (unsorted or sorted).')
    parser.add_argument('--juicer_file', '-j', type=Path, required=True,
                        help='Path to the juicer file to parse.')
    parser.add_argument('--out_file', '-o', required=True, type=str,
                        help='Path of file to write, TSV of "pairs" format')
    parser.add_argument('--sort', '-s', required=False, action="store_true", default=False,
                        help='Whether to sort the file as expected by pairix (takes longer)')
    parser.add_argument('--bam_file', "-b", required=True, type=Path,
                        help="starting bam file, for header info.")
    args = parser.parse_args()
    return args


def strand_convert(strand: str):
    return "+" if strand is "0" else "-"


def get_header_from_bam(bamfile: Path):
    with pysam.AlignmentFile(bamfile, "rb") as bamhandle:
        ref_names = bamhandle.references
        ref_lens = bamhandle.lengths
        out_header = "## pairs format v1.0\n#sorted: chr1-chr2-pos1-pos2\n#shape: upper triangle\n#genome_assembly: hg38"
        for name, length in zip(ref_names, ref_lens):
            out_header += f"\n#chromsize:\t{name}\t{length}\n"
        out_header += "#columns:\treadID\tchr1\tpos1\tchr2\tpos2\tstrand1\tstrand2\n    "

    return out_header


def process_line(line: str):
    """handle juicer long format, convert to Pairix expected
    Juicer long:
    # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2> <mapq1> <cigar1> <sequence1> <mapq2> <cigar2> <sequence2> <readname1> <readname2>
    Pairix:
    ## pairs format v1.0
    #columns: readID chr1 pos1 chr2 pos2 strand1 strand2 <column_name> <column_name>
    readID, chr1, pos1, chr2, pos2, strand1, strand2, [optional:] frag1, frag2
    """
    fields = line.split("\t")
    str1 = strand_convert(fields[0])
    str2 = strand_convert(fields[4])
    chr1 = fields[1]
    chr2 = fields[5]
    pos1 = fields[2]
    pos2 = fields[6]
    read_id = fields[14]
    pairix_line = f"{read_id}\t{chr1}\t{pos1}\t{chr2}\t{pos2}\t{str1}\t{str2}"
    return pairix_line


def main():
    args = parse_args()
    with open(args.juicer_file) as juicer_file:
        with open(args.out_file, "w") as out_file:
            header = get_header_from_bam(args.bam_file)
            out_file.write(header)
            for line in juicer_file:
                print(process_line(line), file=out_file)

    if args.sort:
        print(
            "to sort the output (NYI in this script), I suggest the following commands (linux):"
            "# get header"
            "grep '#' output_file.pairs > output_file.pairs.sorted"
            "# sort output"
            "LC_ALL=C sort -k2,2 -k3,3n -k4,4 -k5,5n --parallel=12 output_file.pairs \""
            " >> output_file.pairs.sorted"
        )


if __name__ == "__main__":
    main()


#!/usr/bin/env python
from __future__ import print_function
import pysam
import argparse
from natsort import humansorted
from collections import defaultdict
import numpy as np

def parse_args():
    parser = argparse.ArgumentParser(description='Parse BAM component filter arguments.')
    parser.add_argument('--bam_file', '-b', type=str, required=True,
                        help='Path to the BAM file(s) to parse. (comma-separated if multiple)')
    parser.add_argument('--out_file', '-o', required=True, type=str, default='counts.tsv',
                        help='Path of file to write, a TSV of format <contig> <contig> <link_count>. Default: counts.tsv')
    parser.add_argument('--out_fmt', '-f', required=False, type=str, default="juicer",
                        help='Format of file to write. Default: juicer (short format). '
                             'others: counts (contigname1 contigname2 counts), '
                             'index_counts (index1 index2 counts) npy (NYI), matrix (NYI).')
    parser.add_argument('--mapq_filter', '-q', required=False, type=int, default=0)
    args = parser.parse_args()
    return vars(args)


# parse bam file to links
def parse_bam_to_link_counts(bamfiles, out_fmt="counts", mapq_filter=0):
    '''Parse the bam file to make it into a Hi-C graph among contig nodes with read edges. (Edges are unweighted).
    Also pull length info for the contigs to filter on later.

    Args:
        bamfiles ([str]): List of bamfile paths to parse (can be length 1).
        out_fmt (str): what the expected output format is.

    Returns:
        {str: {str: bool}}: Representation of the Hi-C graph.
        {str: int}: Mapping of contig names to their sequence length.

    '''
    is_juicer = False
    if out_fmt == "juicer":
        is_juicer = True

    net = defaultdict(lambda: defaultdict(int))
    if len(bamfiles) == 0:
        RuntimeError("No BAM files provided!")
    if not is_juicer:
        print("parsing bam(s)")
    for bamfile in bamfiles:
        bam = pysam.AlignmentFile(bamfile, 'rb')
        num = 0
        num_filtered = 0
        is_first_read = True
        refs = bam.references
        lengths = bam.lengths
        ref_lens = {}
        ref1 = None
        ref2 = None
        mate = ''
        juicer_line = "dummy_entry"
        for read in bam:
            num += 1
            if is_first_read:
                is_first_read = False
                read1_id = read.query_name
                # note that this does not account for mapq of mate read, which appears to be
                # non-trivial to access in pysam (?). So adding iff mate also passes (below).
                mate = [read.next_reference_start, read.next_reference_id]
                ref1 = refs[read.reference_id]
                ref_lens[ref1] = lengths[read.reference_id]
                ref2 = refs[read.next_reference_id]
                ref_lens[ref2] = lengths[read.next_reference_id]

                if read.is_duplicate or read.mapping_quality < mapq_filter:
                    num_filtered += 1
                    continue
                if is_juicer:
                    # juicer short format is as follows:
                    # <str1> <chr1> <pos1> <frag1> <str2> <chr2> <pos2> <frag2>
                    juicer_line = "\t".join(
                                     ["0" if not read.is_reverse else "1",  # strand
                                     ref1,
                                     str(read.reference_start),
                                     "0",  # frag1
                                     "0" if not read.mate_is_reverse else "1",  # strand
                                     ref2,
                                     str(read.next_reference_start),
                                     "1"]  # frag2
                                     )

                if not is_juicer and num % 10000000 == 0:
                    print("parsed {0} reads from file {1}, {2} filtered out".format(
                          num, bamfile, num_dupes)
                          )

            else:
                is_first_read = True
                are_paired = read.query_name == read1_id
                if not are_paired:
                    raise ValueError("BAM file is not sorted by read name!! (samtools sort -n)")
                if mate != [read.reference_start, read.reference_id]:
                    raise RuntimeError("Mates don't match!! (samtools sort -n). Alternately, filtering may have removed reads rather than read pairs?\n"
                                       "offending reads:", read1_id, read.query_name, mate, read.reference_start, read.reference_id)

                # duplicating filter for forward read
                if read.is_duplicate or read.mapping_quality < mapq_filter:
                    num_filtered += 1
                    continue

                # if BOTH reads pass the filter, then increment/print appropriately
                net[ref1][ref2] += 1
                net[ref2][ref1] += 1
                if is_juicer:
                    print(juicer_line)

    if not is_juicer:
        print("parsed {0} reads from file {1}, filtered out {2} on dupes and mapq".format(
              num, bamfile, num_filtered))
        if num_filtered == 0:
            print("number of filtered reads is zero- this may be because the duplicates flag was not"
                  " set in BAM (e.g. by samblaster)")
    return net, ref_lens

def write_net_as_counts(net, outfile):
    '''Write out all those counts from the net dict of defaultdicts to a file.'''
    print("writing output to {0}".format(outfile))
    with open(outfile, 'w') as out:
        for contig1 in humansorted(net.keys()):
            this_dict = net[contig1]
            for contig2 in humansorted(this_dict.keys()):
                outstr = "\t".join([contig1, contig2, str(this_dict[contig2])])
                out.write(outstr + "\n")

def write_net_as_counts_index(net, outfile):
    '''Write out all those counts from the net dict of defaultdicts to a file.
    Using index names to accommodate tools like EVR that expect such.'''
    print("writing counts data to {0}".format(outfile))
    keys = humansorted(net.keys())
    with open(outfile, 'w') as out:
        for contig1 in keys:
            idx_print1 = str(keys.index(contig1))
            this_dict = net[contig1]
            for contig2 in keys:
                idx_print2 = str(keys.index(contig2))
                outstr = "\t".join([idx_print1, idx_print2, str(this_dict[contig2])])
                out.write(outstr + "\n")
    # write out contig names separately for tools that require indices...
    with open(outfile+".names", "w") as outnames:
        outnames.write("\n".join(keys) + "\n")

def write_net_as_npy(net, outfile):
    '''Write out all those counts from the net dict as an npy matrix that can be used by e.g. pastis.
    Requires taking numpy dependency, may or mayn't be worth it.
    '''
    UserWarning("zeroing out diagonal!!!! pastis requires this")
    contigs = humansorted(net.keys())
    counts = np.zeros(shape=(len(contigs), len(contigs)), dtype=int, order="F")
    keys = humansorted(net.keys())
    for contig1 in keys:
        idx1 = keys.index(contig1)
        this_dict = net[contig1]
        for contig2 in this_dict.keys():
            if contig1 == contig2:
                continue
            idx2 = keys.index(contig2)
            counts[idx1, idx2] = this_dict[contig2]
    counts.dump(outfile)
    np.save(outfile, counts)

def write_net_as_matrix(net, outfile):
    '''Write out all those counts from the net dict as a flat matrix file'''
    pass

def read_net_from_counts(counts_file):
    '''Read a previously parsed network from a counts file. Assumes a 3-column file. Not used when run as a script.'''
    print("reading counts into a net from the file {0}.".format(counts_file))
    with open(counts_file) as file:
        net = defaultdict(lambda: defaultdict(int))
        for line in file:
            fields = line.strip().split()
            ref1 = fields[0]
            ref2 = fields[1]
            count = fields[2]
            net[ref1][ref2] = count
            net[ref2][ref1] = count
    return net

def main():
    c_args = parse_args()
    bams = c_args["bam_file"].split(",")
    net, ref_lens = parse_bam_to_link_counts(bams, c_args["out_fmt"], c_args["mapq_filter"])
    #print("writing out data")
    if c_args["out_fmt"] == "counts":
        write_net_as_counts(net=net, outfile=c_args["out_file"])
    elif c_args["out_fmt"] == "index_counts":
        write_net_as_counts_index(net=net, outfile=c_args["out_file"])
    elif c_args["out_fmt"] == "npy":
        write_net_as_npy(net=net, outfile=c_args["out_file"])
    if c_args["out_fmt"] == "juicer":
        pass  # don't need to write this, going to STDOUT

if __name__ == "__main__":
    main()

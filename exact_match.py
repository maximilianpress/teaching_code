#!/usr/bin/env python
"""
Maximilian Press
copyright (C) 2020

This is a very simple mapper that takes:

1) multifasta of query sequences
2) multifasta of reference sequences

and attempts to find exact matches (allowing Ns) of the queries in the references.

Probably very slow for many queries/large references!!!

------

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
from __future__ import print_function
import sys
import re
import gzip
from Bio import SeqIO
from contextlib import contextmanager
from collections import defaultdict

USAGE = "\nexact_match.py <REFERENCE_FASTA> <QUERY_FASTA>\n"
TAB = str.maketrans("ACGTN", "TGCAN")

def rev_comp(seq):
    '''Reverse complement a sequence

    :param seq: str- all ACGTN
    :return: str- the reverse complement
    '''
    return seq.translate(TAB)[::-1]

def seq_to_pat(seq):
    """Make a seq consisting of ACGTN into a regex pattern."""
    seq = seq.upper()
    ambig = 0
    seq_pat = ''
    for nt in seq:
        if nt == "N":
            ambig += 1
            continue
        elif nt != "N" and ambig > 0:
            seq_pat += "[ACGTN]{" + str(ambig) + '}'
            ambig = 0
        seq_pat += nt
    return re.compile(seq_pat)

def format_line(match_dict_entry):
    """Given a match dict entry, make it into a tab-delim line for writing to file."""
    return "{query}\t{ref}\t{start}\t{end}\t{orientation}\t{match_seq}".format(**match_dict_entry)

def print_results(match_dict):
    """Print any results in a match dict"""
    print("query\tref\tmatch_start\tmatch_end\torientation\tmatch_seq")
    for query in match_dict:
        # query is a list of dicts now
        for match in match_dict[query]:
            print(format_line(match))

def find_pats(ref_str, query_str, query_id, seq_id, match_dict, orientation="+"):
    """Mutates match_dict with any matches of the query

    Args:
        ref_str (str): string of a reference sequence
        query_str (str): string of a query sequence
        query_id (str): id of the query seq
        seq_id (str): id of the reference seq
        match_dict (defaultdict(str: [dict])): a dict of matches, to be MUTATED
        orientation (str): "+" or "-", namely whether the query is reverse complemented or not.

    """
    query_pat = seq_to_pat(query_str)
    query_iter = query_pat.finditer(ref_str)
    for match in query_iter:
        coords = match.span()
        match_dict[query_id].append(
            {
             "query" : query_id,
             "ref" : seq_id,
             "start" : coords[0],
             "end" : coords[1],
             "orientation" : orientation,
             "match_seq": ref_str[int(coords[0]):int(coords[1])]
             }
        )

@contextmanager
def open_seq_file(query_seq_file):
    """Figure out GZ status of file, return a Biopython SeqIO.parse iterator over records.
    Separate function to allow context management of opening.

    Args:
        query_seq_file (str): handle of the sequence file.

    Returns:
        filehandle: a successfully opened file.
    """
    if query_seq_file.endswith("gz"):
        handle = gzip.open(query_seq_file, "rt")
    else:
        handle = open(query_seq_file, "r")
    try:
        yield handle
    finally:
        handle.close()

def parse_seq_file(query_seq_file, query_handle):
    """figure out FASTA/Q format and use appropriate SeqIO parser.

    Args:
        query_seq_file (str): handle of the sequence file (for determining filename ending).
        query_handle (str): the actual file handle to read

    Returns:
        iterator over SeqRecords, as from SeqIO.parse
    """
    non_gz = query_seq_file.replace(".gz", "")
    if non_gz.endswith("fa") or non_gz.endswith("fasta"):
        seq_iter = SeqIO.parse(query_handle, "fasta")
    elif non_gz.endswith("fq") or non_gz.endswith("fastq"):
        seq_iter = SeqIO.parse(query_handle, "fastq")
    else:
        e = ValueError("{} not a recognized input file format!! File ending must be one of:\n"
                       ".fa, .fasta, .fq, .fastq (.gz ok)")
        raise e
    return seq_iter

def find_exact_matches(query_seq_file, ref_seq_file, should_rc=True):
    """Iterate over a query seq file, finding exact matches in the ref seq file for each query.

    Arguments are the file paths, and a bool for whether reverse complement should be considered.
    """
    match_dict = defaultdict(lambda: list())
    ref = SeqIO.parse(ref_seq_file, "fasta")
    for seq in ref:
        ref_str = str(seq.seq).upper()
        with open_seq_file(query_seq_file) as query_handle:
            queries = parse_seq_file(query_seq_file=query_seq_file, query_handle=query_handle)
            print("searching for matches in", seq.id)
            for query in queries:
                query_str = str(query.seq).upper()
                # note that find_pats() mutates match_dict
                find_pats(ref_str=ref_str, query_str=query_str, query_id=query.id, seq_id=seq.id,
                          match_dict=match_dict)
                if should_rc:
                    rc = rev_comp(query_str)
                    find_pats(ref_str=ref_str, query_str=rc, query_id=query.id, seq_id=seq.id,
                              match_dict=match_dict, orientation="-")

    return match_dict

def main():
    if len(sys.argv) != 3:
        print(USAGE)
        sys.exit()

    ref_seq_file = sys.argv[1]
    query_seq_file = sys.argv[2]
    match_dict = find_exact_matches(query_seq_file=query_seq_file,
                                    ref_seq_file=ref_seq_file,
                                    should_rc=True)
    print_results(match_dict)

if __name__ == "__main__":
    main()

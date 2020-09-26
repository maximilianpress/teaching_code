#!/usr/bin/env python
"""
Maximilian Press

# this script is 
# free software under the GNU General Public License (GPL). That 
# is, you can do whatever you want with it,
# so long as you don't sell or copyright it.  
# https://www.gnu.org/licenses/gpl-3.0.en.html

"""
import sys
import re
from Bio import SeqIO
from collections import defaultdict

TAB = str.maketrans("ACGTN", "TGCAN")

def rev_comp(seq):
    '''Reverse complement a sequence

    :param seq: str- all ACGTN
    :return: str- the reverse complement
    '''
    return seq.translate(TAB)[::-1]

def seq_to_pat(seq):
    seq = seq.upper()
    ambig = 0
    seq_pat = ''
    for nt in seq:
        if nt == "N":
            ambig += 1
            continue
        elif nt != "N" and ambig > 0:
            seq_pat += "[ACGT]{" + str(ambig) + '}'
            ambig = 0
        seq_pat += nt
    return re.compile(seq_pat)

def format_line(match_dict_entry):
    return "{query}\t{ref}\t{start}\t{end}\t{orientation}\t{match_seq}".format(**match_dict_entry)

def print_results(match_dict):
    print("query\tref\tmatch_start\tmatch_end\torientation\tmatch_seq")
    for query in match_dict:
        # query is a list of dicts now
        for match in match_dict[query]:
            print(format_line(match))

def find_exact_matches(query_seq_file, ref_seq_file):
    match_dict = defaultdict(lambda: list())
    queries = SeqIO.parse(open(query_seq_file), "fasta")
    for query in queries:
        query_str = str(query.seq)
        rc = rev_comp(query_str)
        query_pat = seq_to_pat(query_str)
        rc_pat = seq_to_pat(rc)
        ref = SeqIO.parse(ref_seq_file, "fasta")
        for seq in ref:
            ref_str = str(seq.seq)
            query_iter = query_pat.finditer(ref_str)
            rc_iter = rc_pat.finditer(ref_str)
            for match in query_iter:
                coords = match.span()
                match_dict[query.id].append(
                    {
                        "query" : query.id,
                        "ref" : seq.id,
                        "start" : coords[0],
                        "end" : coords[1],
                        "orientation" : "+",
                        "match_seq": ref_str[int(coords[0]):int(coords[1])]
                    }
                )

            for match in rc_iter:
                match_dict[query.id].append(
                    {
                        "query": query.id,
                        "ref": seq.id,
                        "start": match.span()[0],
                        "end": match.span()[1],
                        "orientation": "-",
                        "match_seq": rev_comp(ref_str[match.span()])
                    }
                )
        ref.close()
    return match_dict

def main():
    ref_seq_file = sys.argv[1]
    query_seq_file = sys.argv[2]
    match_dict = find_exact_matches(query_seq_file=query_seq_file, ref_seq_file=ref_seq_file)
    #print(match_dict)
    print_results(match_dict)

if __name__ == "__main__":
    main()

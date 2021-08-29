# builds heavily upon previous answer, e.g. find_base_in_alignment()!
import sys
import pysam
from typing import Optional
from collections import defaultdict

def find_base_in_alignment(alignment: pysam.AlignedSegment,
                           pos: int, 
                           bam_stores_revcomp: bool = False) -> Optional[str]:
    idx_q = 0
    idx_r = pos - alignment.reference_start
    if bam_stores_revcomp:
        seq = alignment.query_sequence
    else:
        seq = alignment.get_forward_sequence()
    
    if seq is None:
        return None
    
    for op, l in alignment.cigartuples:
        ref_consumed = op in {0, 2, 3, 7, 8}
        query_consumed = op in {0, 1, 4, 7, 8}
        
        if ref_consumed:
            idx_r -= l
        if query_consumed:
            idx_q += l
        
        if idx_r < 0:
            if query_consumed:
                # base is in query between idx_q-l , idx_q
                base = seq[idx_q + idx_r - 1]
                return base
            else:
                # position has been deleted
                return None

def two_variant_cooccurrence(bam: pysam.AlignmentFile, 
                               chrom: str, 
                               coord1 :int, 
                               coord2: int) \
                               -> Optional[defaultdict]:
    allele_cooccurrence = defaultdict(lambda: defaultdict(int))
    all_alns = pysam.AlignmentFile(bam)
    alns = all_alns.fetch(chrom, coord1, coord2)
    ref_seq = ""
    
    for read in alns:
      aln_pos = read.get_reference_positions()
      # not sure how this would happen but good to omit such just in case
      if not coord1 in aln_pos and coord2 in aln_pos:
        continue
      allele1 = find_base_in_alignment(read, coord1)
      allele2 = find_base_in_alignment(read, coord2)
      if allele1 is None or allele2 is None:
          continue
      allele_cooccurrence[allele1][allele2] += 1

    return allele_cooccurrence
    
def write_table(coocc):
    print("position1_allele", "position2_allele", "count")
    for i in coocc:
        for j in coocc[i]:
            print(i, j, coocc[i][j])
    
def main():
    bam_file = open(sys.argv[1])
    chrom = sys.argv[2]
    coord1 = int(sys.argv[3])
    coord2 = int(sys.argv[4])
    coocc = two_variant_cooccurrence(bam_file, chrom, coord1, coord2)
    write_table(coocc)
        
if __name__ == "__main__":
    main()
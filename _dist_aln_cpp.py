# __Auther__ == 'Haowen Xu'
# __Date__ == '02-06-2020'

import smith_waterman_dp
import numpy as np

def local_pairwise_align(seq1, seq2, gap_open_penalty,
                         gap_extend_penalty, substitution_matrix):
    kMaxLen = 1000;
    seq1 = np.array(seq1, dtype=np.int32)
    seq2 = np.array(seq2, dtype=np.int32)
    #mapping = {chr(i+64): i for i in range(1,27)}
    #reverse_mapping = {i: chr(i+64) for i in range(1,27)}
    #reverse_mapping[-1] = '-'
    #seq1 = np.array([mapping[seq1[i]] for i in range(len(seq1))], dtype=np.int32)
    #seq2 = np.array([mapping[seq2[i]] for i in range(len(seq2))], dtype=np.int32)
    aligned_seq1, aligned_seq2, score, seq1_start, seq1_end, seq2_start, seq2_end = \
        smith_waterman_dp.local_pairwise_align(
            seq1, seq2,
            gap_open_penalty,
            gap_extend_penalty,
            substitution_matrix,
            kMaxLen,
            kMaxLen)
    #aligned_seq1 = [reverse_mapping[aligned_seq1[i]] for i in range(len(aligned_seq1)) if aligned_seq1[i]!=0]
    #aligned_seq2 = [reverse_mapping[aligned_seq2[i]] for i in range(len(aligned_seq2)) if aligned_seq2[i]!=0]
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]
    start_end_positions = [(seq1_start, seq1_end),
                           (seq2_start, seq2_end)]
    return aligned_seq1, aligned_seq2, score, start_end_positions



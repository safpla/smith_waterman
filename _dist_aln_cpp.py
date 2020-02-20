# __Auther__ == 'Haowen Xu'
# __Date__ == '02-06-2020'

import smith_waterman_dp
import numpy as np

def local_pairwise_align(seq1, seq2, gap_open_penalty,
                         gap_extend_penalty, substitution_matrix):
    kMaxLen = 1000;
    seq1 = np.array(seq1, dtype=np.int32)
    seq2 = np.array(seq2, dtype=np.int32)
    aligned_seq1, aligned_seq2, score, seq1_start, seq1_end, seq2_start, seq2_end = \
        smith_waterman_dp.local_pairwise_align(
            seq1, seq2,
            gap_open_penalty,
            gap_extend_penalty,
            substitution_matrix,
            kMaxLen,
            kMaxLen)
    # reverse and removing padding
    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq1 = aligned_seq1[aligned_seq1 != -2]
    aligned_seq2 = aligned_seq2[::-1]
    aligned_seq2 = aligned_seq2[aligned_seq2 != -2]

    start_end_positions = [(seq1_start, seq1_end),
                           (seq2_start, seq2_end)]
    return aligned_seq1, aligned_seq2, score, start_end_positions



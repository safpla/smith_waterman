import os, sys
root_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(root_path)
import time
import random
import numpy as np
import _dist_aln_cpp
import _dist_aln

def generate_sequences(n_seqs, len_seq):
    seqs = []
    reverse_mapping = {i: chr(i+64) for i in range(1,27)}
    for _ in range(n_seqs):
        seq = []
        for _ in range(len_seq):
            seq.append(reverse_mapping[random.randint(1,5)])
        seqs.append(seq)
    return seqs


def test_correctness(seqs):
    gap_open_penalty = 3.
    gap_extend_penalty = 1.
    for i in range(len(seqs)):
        for j in range(i+1,len(seqs)):
            seq1 = seqs[i];
            seq2 = seqs[j];
            S = np.zeros([len(seq2), len(seq1)], dtype=np.float)
            for n in range(len(seq2)):
                for m in range(len(seq1)):
                    if seq2[n] == seq1[m]:
                        S[n,m] = 3.
                    else:
                        S[n,m] = -1.
            aligned1, aligned2, score, start_end_positions = \
                _dist_aln.local_pairwise_align(
                    seq1, seq2, gap_open_penalty,
                    gap_extend_penalty, S)
            py_score = score

            aligned1, aligned2, score, start_end_positions = \
                _dist_aln_cpp.local_pairwise_align(
                    seq1, seq2, gap_open_penalty,
                    gap_extend_penalty, S)
            cpp_score = score
            if py_score != cpp_score:
                print('mismatch')
                return
    print('correctness test passed')


def test_speedup(seqs):
    gap_open_penalty = 3.
    gap_extend_penalty = 1.
    start_time = time.time()
    for _ in range(10):
        for i in range(len(seqs)):
            for j in range(i+1,len(seqs)):
                seq1 = seqs[i];
                seq2 = seqs[j];
                S = np.zeros([len(seq2), len(seq1)], dtype=np.float)
                for n in range(len(seq2)):
                    for m in range(len(seq1)):
                        if seq2[n] == seq1[m]:
                            S[n,m] = 3.
                        else:
                            S[n,m] = -1.
                aligned1, aligned2, score, start_end_positions = \
                    _dist_aln.local_pairwise_align(
                        seq1, seq2, gap_open_penalty,
                        gap_extend_penalty, S)
    py_duration = time.time() - start_time
    print('python uses {}'.format(py_duration))
    start_time = time.time()
    for _ in range(10):
        for i in range(len(seqs)):
            for j in range(i+1,len(seqs)):
                seq1 = seqs[i];
                seq2 = seqs[j];
                S = np.zeros([len(seq2), len(seq1)], dtype=np.float)
                for n in range(len(seq2)):
                    for m in range(len(seq1)):
                        if seq2[n] == seq1[m]:
                            S[n,m] = 3.
                        else:
                            S[n,m] = -1.
                aligned1, aligned2, score, start_end_positions = \
                    _dist_aln_cpp.local_pairwise_align(
                        seq1, seq2, gap_open_penalty,
                        gap_extend_penalty, S)
    py_duration = time.time() - start_time
    print('c++ uses {}'.format(py_duration))


if __name__ == "__main__":
    seqs = generate_sequences(2, 500)
    test_correctness(seqs)
    #test_speedup(seqs)


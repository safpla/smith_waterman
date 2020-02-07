# C++ implementation of Smith-Waterman Algorithm callable from Python.
A simple implementation of the Smith-Waterman algorithm for local alignment of sequences written in C++ and wrapped as a python callable function.

## How to Use
* Use SWIG to generate python module wrapper.
```bash
swig -c++ -python smith_waterman_dp.i
```

* Compile
```bash
g++ -fPIC -shared smith_waterman_dp_wrap.cxx -o _smith_water_dp.so -I/home/safpla/.local/lib/python2.7/site-packages/numpy/core/include/  -I/usr/include/python2.7/ -lpython2.7
```
Replace those two path with yours.

* Call in python
```python
import smith_waterman_dp
aligned_seq1, aligned_seq2, score, seq1_start, seq1_end, seq2_start, seq2_end = \
        smith_waterman_dp.local_pairwise_align(
            seq1,
            seq2,
            gap_open_penalty,
            gap_extend_penalty,
            substitution_matrix,
            kMaxLen,
            kMaxLen)
```
To make the input and output format compatible with the original python implementation of Smith-Waterman given in _dist_aln.py, I wrote another wrapper in _dist_aln_cpp.py. You can use it in the following way:
```python
import _dist_aln_cpp

aligned1, aligned2, score, start_end_positions = \
    _dist_aln_cpp.local_pairwise_align(
        seq1,
        seq2,
        gap_open_penalty,
        gap_extend_penalty,
        substitution_matrix)

```



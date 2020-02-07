/* smith_waterman_dp.i */
%module smith_waterman_dp
%{
#define SWIG_FILE_WITH_INIT
#include "smith_waterman_dp.cpp"
%}
%include "numpy.i"
%init %{
import_array();
%}
%apply (int DIM1, int* IN_ARRAY1) {(int n, int *seq1)};
%apply (int DIM1, int* IN_ARRAY1) {(int m, int *seq2)};
%apply (int DIM1, int DIM2, double* IN_ARRAY2) {(int len_seq2, int len_seq1, double *substitution_matrix)};
%apply (int DIM1, int* ARGOUT_ARRAY1) {(int aligned_n, int *aligned_seq1)};
%apply (int DIM1, int* ARGOUT_ARRAY1) {(int aligned_m, int *aligned_seq2)};
%apply double *OUTPUT {double *score};
%apply int *OUTPUT {int *seq1_start};
%apply int *OUTPUT {int *seq1_end};
%apply int *OUTPUT {int *seq2_start};
%apply int *OUTPUT {int *seq2_end};

%include "smith_waterman_dp.cpp"

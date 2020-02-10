#include <iostream>
#include <string>
using namespace std;

#define kN 1000
#define kM 1000
#define kPenalizeTerminalGaps 0

void local_pairwise_align(int n, int *seq1, int m, int *seq2,
        double gap_open_penalty, double gap_extend_penalty,
        int len_seq2, int len_seq1, double *substitution_matrix,
        int aligned_n, int *aligned_seq1,
        int aligned_m, int *aligned_seq2,
        double *score,
        int *seq1_start, int *seq1_end, int *seq2_start, int *seq2_end){

    // Compute score and traceback matrices.
    const int AEND = 0;
    const int MATCH = 1;
    const int VGAP = 2;
    const int HGAP = 3;
    const int UNINIT = -1;
    const double kNewAlignmentScore = 0.0;
    const int kGapCharacter = -1;
    int i;
    int j;

    // score_matrix initialization.
    double score_matrix[kN][kM] = {0.0};
    // traceback_matrix initialization.
    int traceback_matrix[kN][kM];
    traceback_matrix[0][0] = AEND;
    for (i = 1; i < len_seq2 + 1; i++){
        traceback_matrix[i][0] = VGAP;
    }
    for (i = 1; i < len_seq1 + 1; i++){
        traceback_matrix[0][i] = HGAP;
    }
    for (i = 1; i < len_seq2 + 1; i++){
        for (j = 1; j < len_seq1 + 1; j++){
            traceback_matrix[i][j] = UNINIT;
        }
    }
    double substitution_score;
    double best_score;
    double diag_score;
    double up_score;
    double left_score;
    double best_score_in_history = 0.0; // largest in score_matrix
    *seq1_end = 0;
    *seq2_end = 0;
    int traceback;
    // Start dynamic programming
    for (i = 1; i < len_seq2 + 1; i++){
        for (j = 1; j < len_seq1 + 1; j++){
            // compute the score for a match/mismatch
            substitution_score = substitution_matrix[(i-1)*len_seq1 + j-1];
            diag_score = score_matrix[i-1][j-1] + substitution_score;

            // compute the score for adding a gap in aln2 (vertical)
            if (!kPenalizeTerminalGaps && (j == len_seq1))
                up_score = score_matrix[i-1][j];
            else if (traceback_matrix[i-1][j] == VGAP)
                up_score = score_matrix[i-1][j] - gap_extend_penalty;
            else
                up_score = score_matrix[i-1][j] - gap_open_penalty;

            // compute the score for adding a gap in aln1 (horizontal)
            if (!kPenalizeTerminalGaps && (i == len_seq2))
                left_score = score_matrix[i][j-1];
            else if (traceback_matrix[i][j-1] == HGAP)
                left_score = score_matrix[i][j-1] - gap_extend_penalty;
            else
                left_score = score_matrix[i][j-1] - gap_open_penalty;

            best_score = kNewAlignmentScore;
            traceback = AEND;
            
            if (left_score > best_score){
                best_score = left_score;
                traceback = HGAP;
            }
            if (diag_score > best_score){
                best_score = diag_score;
                traceback = MATCH;
            }
            if (up_score > best_score){
                best_score = up_score;
                traceback = VGAP;
            }
            score_matrix[i][j] = best_score;
            traceback_matrix[i][j] = traceback;
            if (best_score > best_score_in_history){
                best_score_in_history = best_score;
                *seq1_end = j;
                *seq2_end = i;
            }
        }
    }
    *score = best_score_in_history;
    //_print_matrix(score_matrix, len_seq2+1, len_seq1+1);
    // Trace back
    // Initialize aligned_seq1, aligned_seq2.
    for (i = 0; i < aligned_n; i++){
        aligned_seq1[i] = 0;
    }
    for (i = 0; i < aligned_m; i++){
        aligned_seq2[i] = 0;
    }
    int current_row = *seq2_end;
    int current_col = *seq1_end;
    int aligned_seq1_pointer = 0;
    int aligned_seq2_pointer = 0;
    int current_value;
    do {
        current_value = traceback_matrix[current_row][current_col];
        switch (current_value)
        {
            case MATCH:
                aligned_seq1[aligned_seq1_pointer] = seq1[current_col-1];
                aligned_seq2[aligned_seq2_pointer] = seq2[current_row-1];
                current_row -= 1;
                current_col -= 1;
                aligned_seq1_pointer += 1;
                aligned_seq2_pointer += 1;
                break;
            case VGAP:
                aligned_seq1[aligned_seq1_pointer] = kGapCharacter;
                aligned_seq2[aligned_seq2_pointer] = seq2[current_row-1];
                current_row -= 1;
                aligned_seq1_pointer += 1;
                aligned_seq2_pointer += 1;
                break;
            case HGAP:
                aligned_seq1[aligned_seq1_pointer] = seq1[current_col-1];
                aligned_seq2[aligned_seq2_pointer] = kGapCharacter;
                current_col -= 1;
                aligned_seq1_pointer += 1;
                aligned_seq2_pointer += 1;
                break;
            default:
                break;
        }
    } while (current_value != AEND);

    *seq1_start = current_col;
    *seq2_start = current_row;
    *seq1_end -= 1;
    *seq2_end -= 1;
}


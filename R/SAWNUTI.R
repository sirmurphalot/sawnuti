source("R/helper_functions/local.R")
source("R/helper_functions/global.R")

sawnuti <- function(.string1, .string2, .times1, .times2, .global = T,
                    .alpha, .match_function, .gap_penalty) {
  #' Compare Sequences with Non-Uniform Time Intervals
  #'
  #' @description The SAWNUTI algorithm performs sequence comparison for finite
  #' sequences of discrete events with non-uniform time intervals.  Further
  #' description of the algorithm can be found in the paper
  #' Comparing Finite Sequences of Discrete Events with Non-Uniform Time Intervals
  #' by Flynt, Murph, and King
  #'
  #'
  #' @param string1 character. Character of representations of each element in the first sequence
  #' seperated by a single white space.
  #' @param string2 character. Character of representations of each element in the first sequence
  #' seperated by a single white space.
  #' @param times1 character. Numerical time values between each value in the first sequence.
  #' Should be written as a single character string, where each value is seperated by a single white
  #' space.
  #' @param times2 character. Numerical time values between each value in the second sequence.
  #' Should be written as a single character string, where each value is seperated by a single white
  #' space.
  #' @param global logical. Logical representing whether a user wishes to perform a global or local
  #' alignment.  T = global, F = local.
  #' @param alpha double. Time interval penalty bias.  Weights the influence of time on the alignment
  #' calculation.
  #' @param match_function R function. Score given for alignment of particular values.  Must be able
  #' to take in two values from the universe of possible events, and return a numerical score
  #' for the alignment of those two elements.  Simpliest is a constant score returned, can also be
  #' implemented as a matrix look-up.
  #' @param gap_penalty numerical. This implementation only allows for constant gap penalties.
  #' @usage sawnuti(string1, string2, times1, times2, global = T, alpha,
  #' match_function, gap_penalty)
  #' @return Vector containing 'similarity score', 'alignment', and 'scoring matrix'.
  #' @details Zero-one scale of all time values in the entire dataset is assumed before calculation.
  #' See Examples for possible formatting of the match_function.
  #' @section PREPROCESSING ASSUMPTION:
  #' Zero-one scale of all time values in the entire dataset is assumed before calculation.  Since
  #' this preprocessing step must be done with ALL possible time intervals, it cannot be done in
  #' this alignment which only takes in two particular observations.
  #' @references Flynt A., King B. R., Murph A. C. (2018). Comparing Finite Sequences of
  #' Discrete Events with Non-Uniform Time Intervals.
  #' @export
  #' @examples
  #' > matchFunction = function(a,b){ifelse(a==b, 1, -1)}
  #'
  #' > sawnuti(string1 = "a b c", string2 = "d b c", "1 2 3",
  #' "3 2 1", alpha = 1, match_function = matchFunction, gap_penalty = 1)
  #'
  #' $ScoreingMatrix
  #'   [,1] [,2] [,3] [,4]
  #' [1,]    0   -3   -5   -6
  #' [2,]   -3   -1   -4   -5
  #' [3,]   -1   -1    0   -2
  #' [4,]   -4   -4   -3   -1
  #'
  #' $AlignmentScore
  #' [1] "-1"
  #'
  #' $Alignment
  #'   [,1] [,2] [,3]
  #' [1,] "d"  "b"  "c"
  #' [2,] "|"  "|"  "|"
  #' [3,] "a"  "b"  "c"
  #' @importFrom base paste

  if(global) {
    matricies = .propogate_matrix_global(string_1 = .string1,
                                        string_2 = .string1,
                                        times_1 = .times1,
                                        times_2 = .times2,
                                        gap_score = .gap_penalty, align_score = .match_function,
                                        time_proportion = .alpha)
    score = .derive_score_global(matricies)
    alignment = .print_comparison_global(matricies, .string1, .string2)
  } else {
    matricies = .propogate_matrix_local(string_1 = .string1,
                                       string_2 = .string1,
                                       times_1 = .times1,
                                       times_2 = .times2,
                                       gap_score = .gap_penalty, align_score = .match_function,
                                       time_proportion = .alpha)
    score = .derive_score_local(matricies)
    alignment = .print_comparison_local(matricies, .string1, .string2)
  }
  results = list(matricies[[1]], score, alignment)
  names(results) = c("ScoreingMatrix", "AlignmentScore", "Alignment")
  return(results)
}

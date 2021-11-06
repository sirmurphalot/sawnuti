################################################
# Global SAWNUTI Algorithm                     #
# Author: Alexander Murph, Bucknell University #
# Last updated: 6 / 11 / 2018                  #
################################################

.propogate_matrix_global = function(string_1, string_2, times_1, times_2,
                            gap_score, align_score,
                            time_proportion) {

  # Split the times and strings into vectors for easy access
  times_1 = as.numeric(strsplit(times_1, ' ',fixed = T)[[1]])
  times_2 = as.numeric(strsplit(times_2, ' ',fixed = T)[[1]])

  string1_list = c('NA', strsplit(string_1, ' ',fixed = T)[[1]])
  string2_list = c('NA', strsplit(string_2, ' ',fixed = T)[[1]])

  # Set up an empty matrix to put scores into, set up an empty dataframe
  # to hold directions of the propagation.
  value_matrix = matrix(0,
                        nrow = (length(string1_list)),
                        ncol = (length(string2_list)))

  source_row_matrix = as.data.frame(value_matrix)
  source_col_matrix = as.data.frame(value_matrix)
  source_gap_score = as.data.frame(value_matrix)
  need_constant_gap = as.data.frame(value_matrix)

  # Use binary to mark when a gap is and is not needed.  In the SAWNUTI implementation,
  # the gap value is only added in a single time for a run of gaps -- all other penaltys
  # come from the time interval penalty function.
  need_constant_gap[1,1] = 1

  # Iterate through the matrix, applying the SAWNUTI
  for(j in 0:(length(string2_list)-1)){

    for(i in 0:(length(string1_list)-1)) {

      score = align_score(string1_list[i+1], string2_list[j+1])

      # For the global approach, the first row and column involves the time penalties
      # of adding in a gap of a given length to the beginning of either sequence.
      if(j == 0 && i == 0){
        next
      } else if (j == 0) {
        above_time_score = .update_time_gap_global(i, j+1, value_matrix, source_row_matrix,
                                           source_col_matrix, times_1, times_2,
                                           time_proportion, gap_score,
                                           source_gap_score, "above")
        above_value = value_matrix[i, j+1] - gap_score * need_constant_gap[i, j+1] +
          source_gap_score[i, j+1] * time_proportion - above_time_score[[1]] * time_proportion
        max_score = above_value
        origin = 3

        # Record our findings in our matrix and data frame.
        value_matrix[i+1,j+1] = max_score
        source_row_matrix[i+1,j+1] = switch(origin, i, i+1, i)
        source_col_matrix[i+1,j+1] = switch(origin, j, j, j+1)
        source_gap_score[i+1, j+1] = above_time_score[[1]]
        need_constant_gap[i+1, j+1] = switch(origin, 1, 0, 0)
      } else if (i == 0) {
        left_time_score = .update_time_gap_global(i+1, j, value_matrix, source_row_matrix,
                                          source_col_matrix, times_1, times_2,
                                          time_proportion, gap_score,
                                          source_gap_score, "left")
        left_value = value_matrix[i+1, j] - gap_score* need_constant_gap[i+1, j] +
          source_gap_score[i+1, j] * time_proportion - left_time_score[[1]] * time_proportion
        max_score = left_value
        origin = 2

        # Record our findings in our matrix and data frame.
        value_matrix[i+1,j+1] = max_score
        source_row_matrix[i+1,j+1] = switch(origin, i, i+1, i)
        source_col_matrix[i+1,j+1] = switch(origin, j, j, j+1)
        source_gap_score[i+1, j+1] = left_time_score[[1]]
        need_constant_gap[i+1, j+1] = switch(origin, 1, 0, 0)
      } else {

        # Score all of the possibilities, and find the max.  Record also which
        # possibility gave you the maximum score.  This 'origin' index is used
        # for the eventual calculation of the Traceback Phase.
        diag_time_score = .update_time_gap_global(i, j, value_matrix, source_row_matrix,
                                          source_col_matrix, times_1, times_2,
                                          time_proportion, gap_score,
                                          source_gap_score, "diag")
        diag_value = value_matrix[i,j] + score + source_gap_score[i, j] * time_proportion -
          diag_time_score[[1]] * time_proportion

        left_time_score = .update_time_gap_global(i+1, j, value_matrix, source_row_matrix,
                                          source_col_matrix, times_1, times_2,
                                          time_proportion, gap_score,
                                          source_gap_score, "left")
        left_value = value_matrix[i+1, j] - gap_score* need_constant_gap[i+1, j] +
          source_gap_score[i+1, j] * time_proportion - left_time_score[[1]] * time_proportion

        above_time_score = .update_time_gap_global(i, j+1, value_matrix, source_row_matrix,
                                           source_col_matrix, times_1, times_2,
                                           time_proportion, gap_score,
                                           source_gap_score, "above")
        above_value = value_matrix[i, j+1] - gap_score * need_constant_gap[i, j+1] +
          source_gap_score[i, j+1] * time_proportion - above_time_score[[1]] * time_proportion

        max_score = max(diag_value, left_value, above_value)
        origin = which.max(c(diag_value, left_value, above_value))

        # Record our findings in our matrix and data frame.
        value_matrix[i+1,j+1] = max_score
        source_row_matrix[i+1,j+1] = switch(origin, i, i+1, i)
        source_col_matrix[i+1,j+1] = switch(origin, j, j, j+1)
        if(left_time_score[[2]] || above_time_score[[2]]) {
          origin = 1
        }
        source_gap_score[i+1, j+1] = switch(origin, 0,
                                            left_time_score[[1]], above_time_score[[1]])
        need_constant_gap[i+1, j+1] = switch(origin, 1, 0, 0)

      }
    }
  }

  return(list(value_matrix, source_row_matrix, source_col_matrix))
}

.update_time_gap_global = function(row_index, column_index, value_matrix, source_row_matrix,
                           source_col_matrix, times_1, times_2, time_proportion, gap_score,
                           source_gap_score, direction) {
  # Calculates the gap/time penalty for a given value in the scoring matrix.  Considers the
  # current longest gap length in the calculation.

  if(row_index == 1 || column_index == 1) {
    # This means we are dealing with a value alongside the edge.
    if(direction == "left") {
      row_index = row_index - 1
    } else if (direction == "above") {
      column_index = column_index - 1
    }
    return( list(abs(sum(times_1[1:row_index]) - sum(times_2[1:column_index])), 0 ))
  } else if ((source_row_matrix[row_index, column_index] < row_index) &&
             (source_col_matrix[row_index, column_index] < column_index)) {
    # This means this value came from our diagonal direction.
    if(direction == "left") {
      return(list(times_2[column_index], 0))
    } else if (direction == "above") {
      return(list(times_1[row_index], 0))
    } else {
      return( list(abs(times_1[row_index] - times_2[column_index]), 0 ))
    }
  } else if (source_row_matrix[row_index, column_index] < row_index) {
    # In this case, this value came from a downward movement, meaning an extended gap in the y-direction.
    if(direction == "left") {
      # This means that our best choice is a 'zigzag' movement.  So, we need to have the algorithm
      # reset the gap score, since we are now going to deal with a gap in the other sequence.
      return(list(abs(source_gap_score[row_index,column_index] - times_2[column_index]), 1) )
    } else if (direction == "above") {
      return(list(source_gap_score[row_index,column_index] + times_1[row_index], 0))
    } else {
      return( list(abs((source_gap_score[row_index,column_index] +
                          times_1[row_index]) - times_2[column_index]), 0) )
    }
  } else if (source_col_matrix[row_index, column_index] < column_index) {
    # In this case, this value came from a rightward movement, meaning an extended gap in the x-direction.
    if(direction == "left") {
      return( list(source_gap_score[row_index,column_index] + times_2[column_index], 0) )
    } else if (direction == "above") {
      # This means that our best choice is a 'zigzag' movement.  So, we need to have the algorithm
      # reset the gap score, since we are now going to deal with a gap in the other sequence.
      return( list(abs(source_gap_score[row_index,column_index] - times_1[row_index]), 1) )
    } else {

      return( list(abs((source_gap_score[row_index,column_index] +
                          times_2[column_index]) - times_1[row_index]), 0) )
    }
  }

}


.find_substring_1_global = function(row_index, col_index, matricies) {
  # Taking in the row and column of the maximum value, we can get
  # the index for each element of the substring.
  if( (row_index == 1) || (col_index == 1) ) {
    return( NA )
  } else {

    next_row = matricies[[2]][row_index, col_index]
    next_col = matricies[[3]][row_index, col_index]

    if ( (next_row < row_index) && (next_col < col_index) ) {
      value = col_index
    } else if (next_col < col_index) {
      value = col_index
    } else {
      value = NA
    }

    return( c(.find_substring_1_global(next_row, next_col, matricies), value) )

  }
}

.find_substring_2_global = function(row_index, col_index, matricies) {
  # Taking in the row and column of the maximum value, we can get
  # the index for each element of the substring.
  if( (row_index == 1) || (col_index == 1) ) {
    return( NA )
  } else {
    next_row = matricies[[2]][row_index, col_index]
    next_col = matricies[[3]][row_index, col_index]
    if ( (next_row < row_index) && (next_col < col_index) ) {
      value = row_index
    } else if (next_col < col_index) {
      value = NA
    } else {
      value = row_index
    }

    return( c(.find_substring_2_global(next_row, next_col, matricies), value) )

  }
}

.print_comparison_global = function(global_matricies, string_1, string_2) {
  # Taking in all of the "source indices" recording throughout the scoring phase,
  # uses the find_substring functions to print out the alignment used.

  matricies = global_matricies
  indicies = list(nrow(matricies[[1]]), ncol(matricies[[2]]))

  string1_list = strsplit(string_1, ' ',fixed = T)[[1]]
  string2_list = strsplit(string_2, ' ',fixed = T)[[1]]

  sub_string1 = string2_list[.find_substring_1_global(indicies[[1]],
                                              indicies[[2]], matricies) -1]
  sub_string2 = string1_list[.find_substring_2_global(indicies[[1]],
                                              indicies[[2]], matricies) -1]

  align_string = rep(" ", times = length(sub_string2))
  align_string = ifelse(is.na(sub_string1), " ", "|")
  align_string = ifelse((is.na(sub_string2) | align_string == " "), " ", "|")

  sub_string1 = ifelse(is.na(sub_string1), '-', sub_string1)
  sub_string2 = ifelse(is.na(sub_string2), '-', sub_string2)

  alignment_matrix = matrix(data = c(sub_string1[-1],
                                     align_string[-1],
                                     sub_string2[-1]),
                            nrow = 3, ncol = (length(sub_string1) - 1),
                            byrow = T)

  return(alignment_matrix)

}

.derive_score_global = function(global_matricies) {
  # Determines the similarity score from the scoring matrix.

  matricies = global_matricies
  nrow_mat = nrow(matricies[[1]])
  ncol_mat = ncol(matricies[[1]])

  return(as.character(matricies[[1]][nrow_mat, ncol_mat]))
}


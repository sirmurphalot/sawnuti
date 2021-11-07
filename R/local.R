################################################
# Local SAWNUTI Algorithm                      #
# Author: Alexander Murph, Bucknell University #
# Last updated: 6 / 11 / 2018                  #
################################################

propogate_matrix_local = function(string_1, string_2, times_1, times_2,
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

  # Iterate through the matrix, apply the Smith-Waterman Algorithm
  for(j in 1:(length(string2_list)-1)){

    for(i in 1:(length(string1_list)-1)) {

      score = align_score(string1_list[i+1], string2_list[j+1])

      # Score all of the possibilities, and find the max.  Record also which
      # possibility gave you the maximum score.  This 'origin' index is used
      # for the eventual calculation of the Traceback Phase.
      diag_time_score = update_time_gap_local(i, j, value_matrix, source_row_matrix,
                                        source_col_matrix, times_1, times_2,
                                        time_proportion, gap_score,
                                        source_gap_score, "diag")
      diag_value = value_matrix[i,j] + score + source_gap_score[i, j] * time_proportion -
        diag_time_score[[1]] * time_proportion

      left_time_score = update_time_gap_local(i+1, j, value_matrix, source_row_matrix,
                                        source_col_matrix, times_1, times_2,
                                        time_proportion, gap_score,
                                        source_gap_score, "left")
      left_value = value_matrix[i+1, j] - gap_score* need_constant_gap[i+1, j] +
        source_gap_score[i+1, j] * time_proportion - left_time_score[[1]] * time_proportion

      above_time_score = update_time_gap_local(i, j+1, value_matrix, source_row_matrix,
                                         source_col_matrix, times_1, times_2,
                                         time_proportion, gap_score,
                                         source_gap_score, "above")
      above_value = value_matrix[i, j+1] - gap_score * need_constant_gap[i, j+1] +
        source_gap_score[i, j+1] * time_proportion - above_time_score[[1]] * time_proportion

      max_score = max(diag_value, left_value, above_value, na.rm = TRUE)
      origin = which.max(c(diag_value, left_value, above_value, 0))

      if( origin == 4 ) {
        origin = 1
      }

      # Record our findings in our matrix and data frame.
      value_matrix[i+1,j+1] = ifelse(max_score < 0, 0, max_score)
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

  return(list(value_matrix, source_row_matrix, source_col_matrix))
}

update_time_gap_local = function(row_index, column_index, value_matrix, source_row_matrix,
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


find_substring_1_local = function(row_index, col_index, matricies) {
  # Taking in the row and column of the maximum value, we can get
  # the index for each element of the substring.
  matricies = matricies
  if( matricies[[1]][row_index,col_index] == 0 ) {
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

    return( c(find_substring_1_local(next_row, next_col, matricies), value) )

  }
}

find_substring_2_local = function(row_index, col_index, matricies) {
  # Taking in the row and column of the maximum value, we can get
  # the index for each element of the substring.
  matricies = matricies
  if( matricies[[1]][row_index,col_index] == 0 ) {
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

    return( c(find_substring_2_local(next_row, next_col, matricies), value) )

  }
}


print_comparison_local = function(local_matricies, string_1, string_2) {
  # Taking in all of the "source indices" recording throughout the scoring phase,
  # uses the find_substring functions to print out the alignment used.
  indicies = which(local_matricies[[1]] == max(local_matricies[[1]],
                                               na.rm = TRUE), arr.ind = TRUE)
  string1_list = strsplit(string_1, ' ',fixed = T)[[1]]
  string2_list = strsplit(string_2, ' ',fixed = T)[[1]]

  sub_string1 = string2_list[find_substring_1_local(indicies[[1]],
                                              indicies[[2]], local_matricies) -1]
  sub_string2 = string1_list[find_substring_2_local(indicies[[1]],
                                              indicies[[2]], local_matricies) -1]

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

derive_score_local = function(local_matricies) {
  # Determines the similarity score from the scoring matrix.

  nrow_mat = which(local_matricies[[1]] == max(local_matricies[[1]]), arr.ind = TRUE)[1]
  ncol_mat = which(local_matricies[[1]] == max(local_matricies[[1]]), arr.ind = TRUE)[2]

  return(as.character(local_matricies[[1]][nrow_mat, ncol_mat]))
}

the <- new.env(parent = emptyenv())
the$number_of_states <- 3
the$all_possible_states <- c(1:the$number_of_states)
the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                   the$number_of_states,
                                   the$number_of_states,
                                   byrow = TRUE)
colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)

#' Report my default settings
#' @export
memento_settings <- function() {
  return(list("number_of_states" = the$number_of_states,
              "all_possible_states" = the$all_possible_states,
              "state_transition_mtx" = the$state_transition_mtx,
              "state_transition_vec_out" = the$state_transition_vec_out))
}

#' Change my settings
#' @param new_number_of_states integer value for number of states: 2, 3 or 4
#' @export
set_memento_settings <- function(new_number_of_states) {
  old_num_states <- the$number_of_states

  the$number_of_states <- new_number_of_states

  the$all_possible_states <- c(1:the$number_of_states)
  the$state_transition_mtx <- matrix(1:the$number_of_states**2,
                                     the$number_of_states,
                                     the$number_of_states,
                                     byrow = TRUE)
  colnames(the$state_transition_mtx) <- rownames(the$state_transition_mtx) <- as.character(1:the$number_of_states)

  the$state_transition_vec_out <- rep(c(1:the$number_of_states), the$number_of_states)

  invisible(old_num_states)
}

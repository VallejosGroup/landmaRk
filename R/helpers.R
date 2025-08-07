.eval_error_str <- function(error_str) {
  if (length(error_str) == 1) {
    stop(error_str)
  } else if (length(error_str) > 1) {
    # If there are multiple errors, format them nicely
    msg <- paste0(
      error_str[1],
      "Additionally, the following errors occurred:\n",
      paste(error_str[-1], collapse = "")
    )
    stop(msg)
  }
}

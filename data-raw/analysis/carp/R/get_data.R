get_data <- function(file) {
  read_csv(file, col_types = cols())
}
l2_dist <- function(x, y) sqrt(sum((x-y)^2))

status_codes <- clarabel::solver_status_descriptions()

if ( requireNamespace("tinytest", quietly = TRUE) ){
  tinytest::test_package("clarabel")
}

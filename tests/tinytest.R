if (requireNamespace("tinytest", quietly = TRUE)) {
  home <- tinytest::at_home()
  tinytest::test_package("investr", at_home = home)
}

context("euclidean_linker_cpp")


test_that(
  "outputs are valid",
  {
    expect_equal({
      input <- as.matrix(data.frame(x = 1:10, y = 1:10))
      crit_dist <- 1.5
      .euclidean_linker_cpp(input, crit_dist)
    },
    rep(11, 10)
    )

    expect_equal({
      set.seed(10)
      input <- as.matrix(data.frame(x = rnorm(10), y = rnorm(10)))
      crit_dist <- 0.4
      .euclidean_linker_cpp(input, crit_dist)
    },
    c(0, 11, 2, 3, 4, 5, 6, 7, 8, 11)
    )
  }
)
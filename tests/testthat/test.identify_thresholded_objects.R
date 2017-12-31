context("identify_thresholded_objects")

test_that(
  "invalid inputs produce an error",
  {
    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7) > 0
      identify_thresholded_objects(mat)
    })

    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      mat[mat < 0.8] <- 0
      identify_thresholded_objects(mat[-(1:7),])
    })

    expect_error({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      mat[mat < 0.8] <- 0
      identify_thresholded_objects(mat[,-(1:10)])
    })
  }
)




test_that(
  "outputs are valid",
  {
    expect_equal({
      set.seed(10)
      mat <- matrix(runif(70), nrow = 7)
      mat[mat < 0.8] <- 0
      identify_thresholded_objects(mat)
    },
    structure(c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, NA, NA, 10, 10, NA, NA, NA, NA, NA, 10, NA, NA,
                NA, NA, NA, NA, 10, NA, NA, 11, 11, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, NA, 6, NA, 7, NA, NA, NA, NA, NA, NA, NA, NA, NA,
                NA, NA, NA, 8, NA, NA, NA, NA, NA), .Dim = c(7L, 10L))
    )
  }
)
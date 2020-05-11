stopifnot(require("testthat"),
          require("selfisher"))

data(haddock)
dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))

test_that("error messages for user-specified start", {
  expect_error(
    selfisher(prop~Lengths, total=tot, dat, psplit=TRUE,
            start=list(betar=c(2))),
    "parameter vector length mismatch.*length\\(betar\\)=1, should be 2")
  expect_error(selfisher(prop~Lengths, total=tot, dat, psplit=TRUE,
                       start=list(junk=5)),
               "unrecognized vector 'junk'")
})

stopifnot(require("testthat"),
          require("selfisher"))

data(haddock)
dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))

context("Saving and loading selfisher objects")

test_that("summary consistency", {
	fm1 <- selfisher(prop~Lengths, p=~1, total=tot, psplit = TRUE, dat)
	s1 <- capture.output(print(summary(fm1)))
    save(fm1, file="fm1.Rdata")
    load("fm1.Rdata")
    file.remove("fm1.Rdata")
    s2 <- capture.output(print(summary(fm1)))
    expect_identical(s1, s2)
})

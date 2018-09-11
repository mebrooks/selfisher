stopifnot(require("testthat"),
          require("selfisher"))

data(haddock)
dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))


context("Very basic selfisher fitting")
m0 <- selfisher(prop~Lengths, p=~0, total=tot, dat)
m1 <- selfisher(prop~Lengths, p=~1, total=tot, dat)
m2 <- selfisher(prop~Lengths, total=tot, dat, cover=FALSE)

test_that("Fixed psplit=0.5", {
	expect_equal(unname(fixef(m0)$p), 0)
	expect_equal(unname(fixef(m0)$r), c(-36.314353,1.233999), tol=1e-3)
	})

test_that("Default is cover=FALSE and estimate psplit", {
	expect_equal(fixef(m1), fixef(m2))
	expect_equal(L50SR(m1), L50SR(m2))
	expect_is(m0, "selfisher")
	expect_is(m1, "selfisher")
})

test_that("AICtab is working", {
  expect_equal(unname(summary(m0)$AICtab), c(119.94059, 122.29669, -57.97029, 36.02903,  34.68329, 22.00000), tol=1e-3)
})

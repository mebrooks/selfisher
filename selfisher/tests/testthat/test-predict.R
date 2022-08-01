stopifnot(require("testthat"),
          require("selfisher"))

data(haddock)
dat <- transform(haddock, tot=nfine+nwide, prop=nwide/(nfine+nwide))
m1 <- selfisher(prop~Lengths, pformula=~1, total=tot, psplit = TRUE, dat)
nd <- data.frame(Lengths=20:50, tot=100)

context("Predict function")

test_that("Different types work", {
	expect_equal(predict(m1, type="response"),
		c(0.005, 0.011, 0.028, 0.065, 0.138,
		0.254, 0.381, 0.477, 0.530, 0.555,
		0.566, 0.570, 0.572, 0.572, 0.573,
		0.573, 0.573, 0.573, 0.573, 0.573,
		0.573, 0.573, 0.573, 0.573), tol=1e-2)
	expect_equal(predict(m1, type="selection"),
		c(0.003, 0.009, 0.021, 0.052, 0.120, 0.254,
		0.460, 0.680, 0.842, 0.930, 0.971, 0.988,
		0.995, 0.998, 0.999, 1.0, 1.0, 1.0, 1.0,
		1.000, 1.000, 1.000, 1.000, 1.000), tol=1e-2)
	expect_equal(predict(m1, type="prob"),
		rep(0.573, nrow(dat)), tol=1e-2)
})

m0 <- selfisher(prop~Lengths, pformula=~0,  psplit = TRUE,total=tot, dat)

test_that("Prediction with fixed psplit works", {
	expect_equal(predict(m0, type="response"),
			predict(m0, type="selection")/(predict(m0, type="selection")+1), tol=1e-5)
	expect_equal(unique(predict(m0, type="prob")), 0.5)
})

test_that("Predict with simple newdata", {
	expect_equal(predict(m0, newdata=nd)[c(1, 15, 30)],
		c(0.0000, 0.4991, 0.5000), tol=1e-3)
	expect_equal(predict(m1, newdata=nd)[c(1, 15, 30)],
		c(0.0001, 0.5656, 0.5728), tol=1e-3)
})

test_that("Predict with splines", {
  library(splines)
  library(selfisher)
  library(ggplot2)
  library(plyr)
  data(comphaddock)

  comphaddock <<- transform(comphaddock,
                          total = TEST1 + TEST2,
                          prop = TEST1 / (TEST1+TEST2),
                          ratio = TEST1 / TEST2
  )
  comphaddock <<- subset(comphaddock, !is.na(prop))
  m5 = selfisher(prop~bs(LENGTH, df = 4), total = total, comphaddock, haul = HAUL)

  newdata1 = data.frame(LENGTH=(0:100)+0.5, total=1, HAUL=NA)
  newdata2 = data.frame(LENGTH=unique(comphaddock$LENGTH), total=1, HAUL=NA)
  newdata2=newdata2[order(newdata2$LENGTH),]

  expect_warning(newdata1$prop <- predict(m5, newdata = newdata1, type = "response"))
  newdata2$prop = predict(m5, newdata = newdata2, type = "response")
  y1 = subset(newdata1, LENGTH%in%newdata2$LENGTH)$prop

  expect_equal(y1, newdata2$prop)
})

test_that("predict works for I(x^2)", {
  mpoly=selfisher(prop~poly(LENGTH, 3), total=total, comphaddock)
  expect_warning(mI <- selfisher(prop~LENGTH + I(LENGTH^2)+ I(LENGTH^3), total=total, comphaddock))
  expect_equal(predict(mpoly), predict(mI), tolerance = 1e-6)
})

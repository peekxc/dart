

testthat::context("Testing dynamic connectivity structures")

library(phtools)

dg <- new(phtools:::DynamicGraph4, 10)


testthat::test_that("Link/cut works", {
	testthat::expect_equal(dg$components(), seq(0, 9))
	testthat::expect_true(dg$link(0,1))
	expect_false(dg$link(0,1))
})



testthat::test_that("Linking clique and removing single edge works", {
	testthat::expect_equal(dg$components(), seq(0, 9))
	dg$link(0,1)
	dg$components()
	dg$print_tours()
	dg$link(0,2)
	dg$components()
	dg$print_tours()
	dg$link(1,2)
	dg$components()
	dg$print_tours()
	dg$cut(1,2)
	dg$components()
	dg$print_tours()
	dg$cut(0,2)
	dg$components()
	dg$print_tours()
	
	
	dg <- new(phtools:::DynamicGraph4, 10)
	dg$components()
	dg$link(0,1)
	dg$link(0,1)
	dg$components()
	dg$print_tours()
	dg$cut(0,1)
	dg$components()
	
	dg$print_tours()
	dg$link(0,1)
	dg$link(0,2)
	dg$link(1,2)
	dg$components()
	dg$print_tours()
	dg$cut(0,2)
	dg$components()
	dg$print_tours()
	dg$cut(0,1)
	dg$components()
	dg$print_tours()
})
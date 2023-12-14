source('../../../../env_setup_tests.R')
source("R/library_imports.R")
library(testthat)

describe("Test utility methods for acceptance calculations", {
    source("R/SBART/MCMC/perturbations/acceptance_ratios.R")

    # mock data
    mock_tree_structure <- list(
        position = c(1, 2, 3, 4, 5, 6, 7, 14, 15),
        parent = c(NA, 1, 1, 2, 2, 3, 3, 7, 7),
        terminal = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE)
    )

    # expected results
    expected_terminal_nodes <- c(4, 5, 6, 8, 9)
    expected_gen2_nodes <- c(2, 7)
    expected_depths <- c(0, 1, 1, 2, 2, 2, 2, 3, 3)

    # actual results
    terminal_nodes <- get_terminal_nodes(mock_tree_structure)
    gen2_nodes <- get_gen2_internal_nodes(mock_tree_structure)
    
    depths <- c(0)
    for (i in 1:length(mock_tree_structure$position)) {
        depths[i] <- get_node_depth(mock_tree_structure, i)
    }

    it("should return the terminal nodes", {
        expect_equal(expected_terminal_nodes, terminal_nodes)
    })

    it("should return the second geenration internal nodes", {
        expect_equal(expected_gen2_nodes, gen2_nodes)
    })

    it("should return the depth of the nodes", {
        expect_equal(expected_depths, depths)
    })
})
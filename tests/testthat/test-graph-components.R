test_that(".tobs_graph_components labels a connected graph as one component", {
  chain <- function(n) {
    a <- matrix(0L, n, n)
    for (s in seq_len(n)) {
      if (s > 1L) a[s, s - 1L] <- 1L
      if (s < n) a[s, s + 1L] <- 1L
    }
    a
  }
  lab <- tulpaObs:::.tobs_graph_components(chain(30L))
  expect_identical(lab, rep(1L, 30L))
  expect_identical(max(lab), 1L)
})

test_that(".tobs_graph_components finds unequal, non-contiguous components", {
  chain <- function(n) {
    a <- matrix(0L, n, n)
    for (s in seq_len(n)) {
      if (s > 1L) a[s, s - 1L] <- 1L
      if (s < n) a[s, s + 1L] <- 1L
    }
    a
  }
  n1 <- 7L; n2 <- 3L
  adj <- matrix(0L, n1 + n2, n1 + n2)
  adj[1:n1, 1:n1] <- chain(n1)
  adj[n1 + 1:n2, n1 + 1:n2] <- chain(n2)

  lab <- tulpaObs:::.tobs_graph_components(adj)
  expect_identical(max(lab), 2L)
  expect_identical(sort(tabulate(lab)), sort(c(n1, n2)))

  # Membership must come from the edges, not from a contiguous equal split: the
  # same graph with its nodes permuted has to give the same partition of nodes.
  set.seed(3)
  perm <- sample(n1 + n2)
  lab_p <- tulpaObs:::.tobs_graph_components(adj[perm, perm])
  expect_identical(max(lab_p), 2L)
  expect_identical(sort(tabulate(lab_p)), sort(c(n1, n2)))
  # Two nodes share a component in the permuted graph exactly when they did
  # before, which pins membership rather than only the component count.
  same_before <- outer(lab[perm], lab[perm], "==")
  same_after  <- outer(lab_p, lab_p, "==")
  expect_identical(same_before, same_after)
})

test_that(".tobs_graph_components isolates a node with no neighbours", {
  adj <- matrix(0L, 5L, 5L)
  adj[1, 2] <- adj[2, 1] <- 1L
  adj[2, 3] <- adj[3, 2] <- 1L
  adj[4, 5] <- adj[5, 4] <- 1L
  # node ordering has the pair 4-5 after an untouched run; every node is placed
  expect_identical(max(tulpaObs:::.tobs_graph_components(adj)), 2L)

  iso <- matrix(0L, 4L, 4L)
  iso[1, 2] <- iso[2, 1] <- 1L
  iso[2, 3] <- iso[3, 2] <- 1L
  lab <- tulpaObs:::.tobs_graph_components(iso)
  expect_identical(max(lab), 2L)
  expect_identical(sort(tabulate(lab)), c(1L, 3L))
})

test_that("a connected graph reports nothing at term construction", {
  chain <- function(n) {
    a <- matrix(0L, n, n)
    for (s in seq_len(n)) {
      if (s > 1L) a[s, s - 1L] <- 1L
      if (s < n) a[s, s + 1L] <- 1L
    }
    a
  }
  adj <- chain(12L)
  # The overwhelmingly common case stays silent, and the term it builds is
  # unchanged by the component report.
  expect_silent(t1 <- tulpaObs:::.tobs_term_icar(graph = adj))
  expect_silent(tulpaObs:::.tobs_term_bym2(graph = adj))
  expect_silent(tulpaObs:::.tobs_term_car_proper(graph = adj))
  expect_identical(t1$n_units, 12L)
  expect_identical(tulpaObs:::.tobs_report_graph_components(adj, "icar"), NULL)
})

test_that("a disconnected graph reports its component count and sizes", {
  chain <- function(n) {
    a <- matrix(0L, n, n)
    for (s in seq_len(n)) {
      if (s > 1L) a[s, s - 1L] <- 1L
      if (s < n) a[s, s + 1L] <- 1L
    }
    a
  }
  adj <- matrix(0L, 10L, 10L)
  adj[1:6, 1:6] <- chain(6L)
  adj[7:10, 7:10] <- chain(4L)

  expect_message(tulpaObs:::.tobs_term_icar(graph = adj), "2 connected components")
  expect_message(tulpaObs:::.tobs_term_icar(graph = adj), "sizes 6, 4")
  expect_message(tulpaObs:::.tobs_term_icar(graph = adj), "one sum-to-zero per component")
  # Every areal entry point goes through the same check.
  expect_message(tulpaObs:::.tobs_term_bym2(graph = adj), "2 connected components")
  expect_message(tulpaObs:::.tobs_term_car(graph = adj), "2 connected components")
  expect_message(tulpaObs:::.tobs_term_car_proper(graph = adj), "2 connected components")

  # The report is the only difference: the term itself is what it always was.
  t_msg <- suppressMessages(tulpaObs:::.tobs_term_icar(graph = adj))
  expect_identical(t_msg$n_units, 10L)
  expect_identical(t_msg$type, "icar")
})

test_that("a single-node component shows in the sizes and is not pinned to zero", {
  adj <- matrix(0L, 6L, 6L)
  adj[1, 2] <- adj[2, 1] <- 1L
  adj[2, 3] <- adj[3, 2] <- 1L
  adj[4, 5] <- adj[5, 4] <- 1L
  # node 6 has no neighbours, so it is its own component and the sizes say so
  expect_message(tulpaObs:::.tobs_term_icar(graph = adj), "3 connected components")
  expect_message(tulpaObs:::.tobs_term_icar(graph = adj), "sizes 3, 2, 1")

  # The areal cover paths reject an isolated node outright; that guard stands and
  # is the more specific message, so the component report does not restate it.
  expect_error(tulpaObs:::.occu_cover_icar_Q(adj), "isolated node")
  expect_error(tulpaObs:::.occu_cover_adj_to_csr(adj), "isolated node")

  # Where such a node IS accepted it keeps a proper N(0, 1/tau) effect: the
  # identification augments the precision, so the charge for value v on a size-1
  # component is -0.5 * tau * v^2. A hard sum-to-zero would pin it to exactly 0.
  skip_if_not_installed("tulpa")
  csr <- tulpaObs:::adjacency_to_csr(adj)
  lp <- function(x, tau) {
    tulpa:::cpp_test_log_prior_icar(x, nrow(adj), tau, csr$row_ptr,
                                    csr$col_idx, csr$n_neighbors)
  }
  tau <- 1.3
  x0 <- c(0, 0, 0, 0, 0, 0)
  x2 <- c(0, 0, 0, 0, 0, 2)
  expect_equal(lp(x0, tau) - lp(x2, tau), 0.5 * tau * 4, tolerance = 1e-8)
})

test_that("the engine pins one sum-to-zero per component, not one globally", {
  skip_if_not_installed("tulpa")
  chain <- function(n) {
    a <- matrix(0L, n, n)
    for (s in seq_len(n)) {
      if (s > 1L) a[s, s - 1L] <- 1L
      if (s < n) a[s, s + 1L] <- 1L
    }
    a
  }
  n1 <- 20L; n2 <- 12L
  A1 <- chain(n1); A2 <- chain(n2)
  adj <- matrix(0L, n1 + n2, n1 + n2)
  adj[1:n1, 1:n1] <- A1
  adj[n1 + 1:n2, n1 + 1:n2] <- A2

  lp <- function(x, A, tau) {
    csr <- tulpaObs:::adjacency_to_csr(A)
    tulpa:::cpp_test_log_prior_icar(x, nrow(A), tau, csr$row_ptr,
                                    csr$col_idx, csr$n_neighbors)
  }
  tau <- 1.3
  set.seed(21)
  x1 <- rnorm(n1); x1 <- x1 - mean(x1)
  x2 <- rnorm(n2); x2 <- x2 - mean(x2)
  x  <- c(x1, x2)

  # The pooled log-prior is exactly the sum of the two independent ones.
  expect_equal(lp(x, adj, tau), lp(x1, A1, tau) + lp(x2, A2, tau),
               tolerance = 1e-10)

  # Raise component 1 and lower component 2, leaving the GLOBAL sum at zero. One
  # global constraint would not charge for this at all; one per component charges
  # tau * (n1 + n2) * c^2 / 2.
  cc <- 0.3
  d <- lp(c(x1 + cc, x2 - cc), adj, tau) - lp(x, adj, tau)
  expect_equal(d, -0.5 * tau * (n1 * cc^2 + n2 * cc^2), tolerance = 1e-8)
  expect_lt(d, -1e-6)
})

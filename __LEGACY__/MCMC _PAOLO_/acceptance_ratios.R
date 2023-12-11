getTerminalNodes <- function(tree) {
    terminal_nodes <- which(tree$terminal)
    return(terminal_nodes)
}

getGen2InternalNodes <- function(tree) {
    terminal_nodes <- getTerminalNodes(tree)
    gen2_internal_nodes <- which(table(tree$parent[terminal_nodes]) == 2)
    return(gen2_internal_nodes)
}

getNodeDepth <- function(tree, node_index) {
    position <- tree$position[node_index]
    depth <- floor(log2(position))
    return(depth)
}

logTransitionRatio.Grow <- function(
    p.grow, 
    p.prune,
    tree.current,
    tree.new
) {
    terminal_nodes.current.count <- length(getTerminalNodes(tree.current))    
    gen2_internal_nodes.new.count <- length(getGen2InternalNodes(tree.new))

    # missing the probability of picking a predictor, should simplify
    log_transition_ratio <- log(p.grow) - log(p.prune) + log(terminal_nodes.current.count) - log(gen2_internal_nodes.new.count)
    return(log_transition_ratio)
}

logLikelihoodRatio.Grow <- function(
    sigma2,
    sigma_mu2,
    residuals, 
    obs.left, 
    obs.right
) {
    n.left <- length(obs.left)
    n.right <- length(obs.right)
    n.all <- n.left + n.right

    variance.all <- sigma2 + n.all * sigma_mu2
    variance.left <- sigma2 + n.left * sigma_mu2
    variance.right <- sigma2 + n.right * sigma_mu2

    obs.all <- union(obs.left, obs.right)
    residuals_sum2.all <- sum(residuals[obs.all])^2
    residuals_sum2.left <- sum(residuals[obs.left])^2
    residuals_sum2.rigth <- sum(residuals[obs.right])^2

    sqrt_term <- 0.5 * (log(sigma2) + log(variance.all) - log(variance.left) - log(variance.right))
    
    exp_term <- 0.5 * sigma_mu2 / sigma2 * (
        residuals_sum2.left / variance.left +
        residuals_sum2.rigth / variance.right - 
        residuals_sum2.all / variance.all 
    )

    log_likelihood_ratio <- sqrt_term + exp_term
    return(log_likelihood_ratio)
}

logTreeStructureRatio.Grow <- function(
    alpha, 
    beta,
    node_to_grow.index
    tree.current
) {
    depth <- getNodeDepth(tree.current, node_to_grow.index)
    
    numerator <- 1 - alpha / (2 + depth)^beta
    denominator <-(1 + depth)^beta - alpha

    # missing the probability of picking a predictor, should simplify
    log_structure_ratio <- log(alpha) + 2 * log(numerator) - log(denominator) 
    return(log_structure_ratio)
}

logTransitionRatio.Prune <- function(
    p.grow, 
    p.prune,
    tree.current
) {
    terminal_nodes.current.count <- length(getTerminalNodes(tree.current))    
    gen2_internal_nodes.current.count <- length(getGen2InternalNodes(tree.current))

    # missing the probability of picking a predictor, should simplify
    log_transition_ratio <- log(p.prune) - log(p.grow) + log(gen2_internal_nodes.current.count) - log(terminal_nodes.current.count - 1) 
    return(log_transition_ratio)
}

# This is just the inverse of the likelihood ratio for the GROW proposal
logLikelihoodRatio.Prune <- function(
    sigma2,
    sigma_mu2,
    residuals, 
    obs.left, 
    obs.right
) {
    log_likelihood_ratio <- -logLikelihoodRatio.Grow(sigma2, sigma_mu2, residuals, obs.left, obs.right)
    return(log_likelihood_ratio)
}

# This is just the inverse of the structure ratio for the GROW proposal
logTreeStructureRatio.Prune <- function(
    alpha, 
    beta,
    node_to_grow.index
    tree.current
) {
    log_structure_ratio <- -logLikelihoodRatio.Grow(alpha, beta, node_to_grow.index, tree.current)
    return(log_structure_ratio)
}

# Transition and structure ratios cancel out
logAcceptanceRatio.Change <- function(
    sigma2, 
    sigma_mu2,
    obs.left,
    obs.right,
    obs_star.left,
    obs_star.right
) {
    n.left <- length(obs.left)
    n.right <- length(obs.right)
    n_star.left <- length(obs_star.left)
    n_star.right <- length(obs_star.right)

    variance.left <- sigma2 + n.left
    variance.right <- sigma2 + n.right
    variance_star.left <- sigma2 + n_star.left
    variance_star.right <- sigma2 + n_star.right

    residuals_sum2.left <- sum(residuals[obs.left])
    residuals_sum2.right <- sum(residuals[obs.right])
    residuals_star_sum2.left <- sum(residuals[obs_star.left])
    residuals_star_sum2.right <- sum(residuals[obs_star.right])

    sqrt_term <- 0.5 * (log(variance.left) + log(variance.right) - log(variance_star.left - log(variance_star.right)))

    exp_term <- 0.5 * sigma_mu2 / sigma2 * (
        residuals_star_sum2.left / variance_star.left + 
        residuals_star_sum2.right / variance_star.right - 
        residuals_sum2.left / variance.left - 
        residuals_sum2.right / variance.right
    )

    log_likelihood_ratio <- sqrt_term + exp_term
    return(log_likelihood_ratio)
}

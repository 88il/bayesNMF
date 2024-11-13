library(truncnorm)
library(fitdistrplus)

# 79 x 96 matrix LOAD
cosmic_signatures <- t(cosmic_matrix)

# Density function for truncated normal
dtnorm <- function(x, mean, sd, lower = 0, upper = 1) {
    dtruncnorm(x, a = lower, b = upper, mean = mean, sd = sd)
}

# Distribution function for truncated normal
ptnorm <- function(q, mean, sd, lower = 0, upper = 1) {
    ptruncnorm(q, a = lower, b = upper, mean = mean, sd = sd)
}

# Quantile function for truncated normal
qtnorm <- function(p, mean, sd, lower = 0, upper = 1) {
    qtruncnorm(p, a = lower, b = upper, mean = mean, sd = sd)
}


# Get NULL Dist by fitting MLE of mvntruncnorm where covariance matrix is
# diagonal (signatures are independent in H_0)

# Function to fit truncated normal to each mutation type
fit_truncnorm <- function(data) {
    fit <- fitdist(data, "tnorm",
                   start = list(mean = mean(data), sd = sd(data)),
                   lower = c(0, 0),  # Lower bounds for mean and sd
                   upper = c(1, Inf))  # Upper bound for mean, no upper bound for sd
    return(fit$estimate)
}

# Fit truncated normal to each mutation type
fitted_params <- apply(cosmic_signatures, 2, fit_truncnorm)

# Extract means and sds... MLE
mus <- fitted_params[1,]
sigmas <- fitted_params[2,]



mutation_types <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

contexts <- paste0(rep(c("A", "C", "G", "T"), each = 24),
                   "[", rep(mutation_types, each = 4), "]",
                   rep(c("A", "C", "G", "T"), times = 24))

# Ensure we have 96 contexts
stopifnot(length(contexts) == 96)

# Function to sample signatures assuming independent mutation types
sample_independent_signature <- function(n_mutations = 10000) {
    probs <- rtruncnorm(96, a = 0, b = 1, mean = mus, sd = sigmas)
    probs <- probs / sum(probs)  # Normalize to ensure sum is 1

    sampled_mutations <- sample(contexts, size = n_mutations, replace = TRUE, prob = probs)
    mutation_counts <- table(factor(sampled_mutations, levels = contexts))

    signature <- as.vector(mutation_counts) / sum(mutation_counts)
    names(signature) <- contexts
    return(signature)
}

# # Function to compute correlation matrix
# compute_correlation_matrix <- function(signatures) {
#     cor(t(signatures))
# }

# Observed correlation matrix
observed_cor <- cor(cosmic_signatures)

# Permutation test
n_permutations <- 1000
n_signatures <- 79
null_distributions <- array(0, dim = c(96, 96, n_permutations))

set.seed(42)  # For reproducibility

for (i in 1:n_permutations) {
    sampled_signatures <- replicate(n_signatures, sample_independent_signature(), simplify = FALSE)
    sampled_signatures <- do.call(cbind, sampled_signatures)
    null_distributions[,,i] <- compute_correlation_matrix(sampled_signatures)
}

# Compute permutation p-values
p_values <- matrix(0, nrow = 96, ncol = 96)

for (i in 1:96) {
    for (j in 1:96) {
        null_dist <- null_distributions[i,j,]
        p_values[i,j] <- mean(abs(null_dist) >= abs(observed_cor[i,j]))
    }
}

# Visualize results
library(pheatmap)

pheatmap(p_values,
         main = "Permutation p-values for COSMIC signature correlations",
         color = colorRampPalette(c("red", "white"))(100),
         breaks = seq(0, 1, length.out = 101,),
         cluster_rows=FALSE, cluster_cols=FALSE)

# Adjust for multiple testing
p_values_adj <- p.adjust(p_values, method = "BH")
significant_correlations <- which(p_values_adj < 0.05, arr.ind = TRUE)

pheatmap(p_values_adj,
         main = "Permutation adjusted p-values for COSMIC signature correlations",
         color = colorRampPalette(c("red", "white"))(100),
         breaks = seq(0, 1, length.out = 101,))#,
         #cluster_rows=FALSE, cluster_cols=FALSE)

print(paste("Number of significant correlations out of 96 x 96 = 9216:", length(significant_correlations)))


## WHICH PVALUES ARE THE MOST SIGNIFICANT?

# Flatten the matrix and get indices of the smallest p-values
flat_indices <- order(p_values, decreasing = FALSE)[1:50]
smallest_values <- p_values[flat_indices]
indices <- arrayInd(flat_indices, dim(p_values))

# Combine results
results <- data.frame(Value = smallest_values, Row = indices[,1], Column = indices[,2])
print(results)

# Print the contexts corresponding to these indices
for (i in 1:5) { # nrow(indices)
    cat("P-value:", smallest_values[i],
        "Row Context:", contexts[indices[i, 1]],
        "Column Context:", contexts[indices[i, 2]], "\n")
}



## ANALYZE CLUSTERS
# Generate the heatmap and capture the result
heatmap_result <- pheatmap(p_values,
                           main = "Permutation p-values for COSMIC signature correlations",
                           color = colorRampPalette(c("red", "white"))(100),
                           breaks = seq(0, 1, length.out = 101),
                           cluster_rows = TRUE, cluster_cols = TRUE)

# Extract row order from the heatmap result
ordered_row_indices <- heatmap_result$tree_row$order

# Print mutation types in clustered order
ordered_contexts <- contexts[ordered_row_indices]
print(ordered_contexts)

library('dplyr')
library('tidyr')
library('broom')
library('purrr')
library('ggplot2')
theme_set(theme_minimal())


simulateData <- function(alpha, beta, gamma, n) {
    epsilon <- rnorm(n)
    nu <- rnorm(n)
    x <- rnorm(n)

    d <- x*gamma + nu
    y <- d*alpha + x*beta + epsilon

    return(tibble(y, d, x))
}

estimateModels <- function(data) {
    long <- lm(y ~ d + x, data)
    short <- lm(y ~ d, data)
    corr <- lm(d ~ x, data)
    list(long = long, short = short, corr = corr) %>% map(tidy)
}

getRelevantParams <- function(models, n) {
    list(
        beta_p = filter(models$long, term == "x")$p.value,
        long = filter(models$long, term == "d")$estimate,
        long_p = filter(models$long, term == "d")$p.value,
        long_se = filter(models$long, term == "d")$std.error,
        short = filter(models$short, term == "d")$estimate,
        short_p = filter(models$short, term == "d")$p.value,
        short_se = filter(models$short, term == "d")$std.error,
        gamma_p = filter(models$corr, term == "x")$p.value,
        n = n
    )
}

simulateEstimation <- function(alpha, beta, gamma = 0.8, n = 100) {
    simulateData(alpha, beta, gamma, n) %>% 
    estimateModels() %>% 
    getRelevantParams(n)
}

calculateFinalEstimates <- function(results, method, ttest_p_cutoff = 0.05, eyeball_cutoff) {
    if (method == 'single') {
        condition <- quo(beta_p > ttest_p_cutoff)
    } else if (method == 'single-eyeball') {
        condition <- quo(beta_p > ttest_p_cutoff & abs(long - short) <= eyeball_cutoff)
    } else if (method == 'double') {
        condition <- quo(beta_p > ttest_p_cutoff & gamma_p > ttest_p_cutoff)
    } else {
        stop("Invalid value for method.")
    }
    
    mutate(results,
        final = ifelse(!! condition, short, long),
        final_p = ifelse(!! condition, short_p, long_p), 
        final_se = ifelse(!! condition, short_se, long_se),
        final_df = ifelse(!! condition, n - 2, n - 3)
    )
}

# SIMULATION - takes a while depending on the speed of your computer
set.seed(201706)
n_sim <- 10000
simulation_results <- map_df(seq(n_sim), ~{simulateEstimation(alpha = 0, beta = 0.2)})
methods <- list("single", "single-eyeball", "double")
simulation_results_methods <- map_df(methods, function(m) {
    simulation_results %>%
    calculateFinalEstimates(method = m, ttest_p_cutoff = 0.05, eyeball_cutoff = 0.05) %>%
    mutate(method = m)
}) %>% mutate(method = ordered(method, methods))


simulation_results_methods %>%
    select(method, long, final) %>%
    gather(model, estimate, -method) %>%
    ggplot(aes(estimate, fill = model)) + 
    geom_density(alpha = 0.5, linetype = 0) + 
    facet_wrap(~method)
ggsave("figure/coefficients.png", scale = 0.7)

simulation_results_methods %>%
    ggplot(aes(final_se, long_se)) + 
    geom_point(size = 2, alpha = 0.1) +
    geom_abline(slope = 1, linetype = 'dashed') +
    xlab("Final model") +
    ylab("Correct model") +
    coord_fixed(xlim = c(0.05, 0.15), ylim = c(0.05, 0.15)) +
    facet_wrap(~method)
ggsave("figure/standard-errors.png", scale = 0.7)


calculateERP <- function(df) {
    summarize(df, mean(final_p < 0.05) - 0.05) %>% 
    as.numeric()
}

ERPs <- map_df(
    seq(0, 0.2, 0.005), function(c) {
    calculateFinalEstimates(simulation_results, "single-eyeball", eyeball_cutoff = c) %>% 
        calculateERP() %>%
        tibble(single_selection_eyeball = ., tolerance = c)
    }) %>%
    mutate(
        double_selection = calculateFinalEstimates(simulation_results, "double") %>% calculateERP(),
        single_selection = max(single_selection_eyeball)
    ) %>%
    gather(method, ERP, -tolerance)

ERPs %>%
    ggplot(aes(x = tolerance, y = ERP, color = method)) + 
    geom_line(size = 1) +
    ylab("Error in rejection probability (ERP)")
ggsave("figure/erp.png", scale = 0.7)

calculateRejectionRate <- function(df, null_hyp) {
    mutate(df, final_tstat = abs((final - null_hyp)/final_se)) %>% 
    summarize(mean(final_tstat > abs(qt(0.025, final_df)))) %>% 
    as.numeric()
}

simulation_results_all <- simulation_results %>%
    mutate(effect = 'no') %>%
    rbind(
        map_df(
            seq(n_sim), 
            ~{simulateEstimation(alpha = 0.2, beta = 0)}
        ) %>% mutate(effect = 'yes')
    )

map_df(methods, function(m) {
    map_df(seq(0, 0.4, 0.005), function(null_hyp) {
        filter(simulation_results_all, effect == 'yes') %>%
        calculateFinalEstimates(method = m, ttest_p_cutoff = 0.05, eyeball_cutoff = 0.01) %>%
        calculateRejectionRate(., null_hyp) %>%
        tibble(rejection = ., null = null_hyp, method = m)    
    })
}) %>%
    ggplot(aes(x = null, y = rejection, col = method)) +
    geom_line(size = 1) +
    scale_y_continuous(lim = c(0, NA), label = scales::percent) +
    ylab("Rejection rate") +
    xlab("Null hypothesis")
ggsave("figure/power.png", scale = 0.7)

map_df(methods, function(m) {
    simulation_results_all %>% 
    group_by(effect) %>%
    do(
        calculateFinalEstimates(., m, eyeball_cutoff = 0.01) %>% 
        calculateRejectionRate(., 0) %>%
        tibble(rejection_rate = .)
    ) %>%
    mutate(method = m)
}) %>%
    spread(effect, rejection_rate) %>%
    ggplot(aes(no, yes, color = method)) +
    geom_point(size = 2) +
    geom_abline(slope = 1, linetype = 'dashed') +
    ylab("True Positive Rate") +
    xlab(paste("False Positive Rate")) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1))
ggsave("figure/roc.png", scale = 0.7)

source('~/Dropbox (Partners HealthCare)/github_repo/aou_gwas/R/constants.R')
packages = c('GGally', 'reshape2', 'data.table', 'scales', 'Matrix', 'ggplot2', 'extrafont', 'gridExtra', 'grDevices', 'grid',
             'RCurl', 'tidyverse',  'devtools', 'broom', 'plotly', 'slackr', 'magrittr', 'gapminder', 'readr',
             'purrr', 'skimr', 'gganimate', 'gghighlight', 'plotROC', 'naniar', 'BiocManager', 'cowplot', 'corrplot', 'corrr', 'ggridges', 'RColorBrewer',
             'ggpubr', 'pbapply', 'RMySQL', 'egg', 'ggwordcloud', 'patchwork', 'ggthemes', 'ggrepel', 'dplyr')

for(p in packages){
  if(!require(p, character.only = T)){
    install.packages(p)
  }
}

################################################################################################
# Github repo: https://github.com/mjpirinen/linemodels
# Code example: https://github.com/mjpirinen/linemodels/blob/main/linemodels_examples.R
# Documentation: https://github.com/mjpirinen/linemodels/blob/main/linemodels_vignette.pdf
## Line model example
# install.packages('devtools')
# library(devtools)
# install_github("mjpirinen/linemodels")
library(linemodels)

###################### Example ########################
# visualize.line.models(scales, slopes, cors,
#                       model.names = c("M1","M.5","M0"),
#                       model.cols = c("dodgerblue", "orange", "red"),
#                       legend.position = "topleft",
#                       xlim = c(0, 0.6), ylim = c(0, 0.6),
#                       xlab = expression(beta[1]),
#                       ylab = expression(beta[2]),
#                       emphasize.axes = FALSE)


## Load Data
get_pheno_label <- function(pheno){
  label <- if_else(pheno == 'spirometry', 'FEV[1]/FVC', str_to_sentence(pheno))
  label <- if_else(label == 'Cigday', 'Cig/Day', label)
  label <- if_else(label == 'Lung_cancer', 'Lung Cancer', label)
  label <- if_else(label == 'Copd', toupper(label), label)
  return(label)
}

load_clumped_data <- function(pheno1, pheno2, pop, root_path = '~/Downloads/'){
  label1 <- get_pheno_label(pheno1)
  label2 <- get_pheno_label(pheno2)
  pop1 = pop2 = pop
  if(pop == 'meta' & pheno1 %in% c('eosinophil', 'neutrophil', 'cigday', 'smkinit')) pop1 = 'eur'
  if(pop == 'meta' & pheno2 %in% c('eosinophil', 'neutrophil', 'cigday', 'smkinit')) pop2 = 'eur'

  print(paste0('-------Loading full sumstats comparison data (', pheno1, '-', pheno2, '-', pop, ')-------'))
  data <- read_delim(paste0(root_path, pheno1, '_', pop1, '_', pheno2, '_', pop2,'.txt.bgz'), delim = '\t', col_types = cols(locus = col_character())) %>%
    filter(!is.na(get('BETA_1')) & !is.na(get('BETA_2'))) %>%
    mutate(label = factor(if_else(get('P_1') < 5e-8 & get('P_2') < 5e-8, 'Both',
                                  if_else(get('P_1') < 5e-8, label1, label2)),
                          levels = c('Both', label1, label2)),
    )
  print(paste0('Number of variants:', nrow(data)))

  path1 <- paste0(root_path, pheno1, '_clumped_sumstats_GRCh37.txt.bgz')
  path2 <- paste0(root_path, pheno2, '_clumped_sumstats_GRCh37.txt.bgz')
  if(file.exists(path1) & file.exists(path2) & nrow(data) > 10000){
    # print(paste0('-------Loading clumped ', label1 ,' sumstats data-------'))
    pheno1_hit <- read_delim(path1, delim = '\t', col_types = cols(locus = col_character())) %>%
      select(locus, alleles) %>%
      mutate(pheno1_top_hit = T)

    # print(paste0('-------Loading clumped ', label2 ,' sumstats data-------'))
    pheno2_hit <- read_delim(path2, delim = '\t', col_types = cols(locus = col_character())) %>%
      select(locus, alleles) %>%
      mutate(pheno2_top_hit = T)

    data <- data %>%
      merge(., pheno1_hit, by = c('locus', 'alleles'), all = T) %>%
      merge(., pheno2_hit, by = c('locus', 'alleles'), all=T) %>%
      mutate(pheno1_top_hit = if_else(is.na(pheno1_top_hit), FALSE, pheno1_top_hit),
             pheno2_top_hit = if_else(is.na(pheno2_top_hit), FALSE, pheno2_top_hit),) %>%
      filter(pheno1_top_hit | pheno2_top_hit) %>%
      filter(complete.cases(.))
    print(paste0('Number of clumped variants:', nrow(data)))
  }
  return(data)
}

relabel_data <- function(data, pheno1, pheno2){
  data1 <- data %>% filter(label %in% c('Both', get_pheno_label(pheno1))) %>% mutate(label = get_pheno_label(pheno1))
  data2 <- data %>% filter(label %in% c('Both', get_pheno_label(pheno2))) %>% mutate(label = get_pheno_label(pheno2))
  data <- rbind(data1, data2)
  return(data)
}

## Estimate slope with EM
em_algorithm <- function(x, y, beta_init, tol = 1e-6, max_iter = 1000) {
  beta <- beta_init
  converged <- FALSE
  iter <- 0
  n <- length(y)

  while (!converged && iter < max_iter) {
    iter <- iter + 1

    # E-Step: Estimate missing y values
    y_hat <- ifelse(is.na(y), beta * x, y)

    # M-Step: Update beta
    beta_new <- sum(x * y_hat) / sum(x^2)

    # Check for convergence
    if (abs(beta_new - beta) < tol) {
      converged <- TRUE
    }

    beta <- beta_new
  }

  return(list(beta = beta, iterations = iter, converged = converged))
}

estimate_beta_from_em <- function(data, beta1, beta2, init_slopes, labels){
  slopes <- c()
  for(pheno in labels){

    data_to_fit <- data %>% filter(label %in% c(pheno))

    x <- unname(unlist(data_to_fit[, beta1]))
    y <- unname(unlist(data_to_fit[, beta2]))

    # Run the EM algorithm
    result <- em_algorithm(x, y, init_slopes[pheno])

    # Print the results
    if(! result$converged){
      print(paste("Estimated beta:", result$beta))
      print(paste("Number of iterations:", result$iterations))
      print(paste("Converged:", result$converged))
    }
    slopes <- c(slopes, result$beta)
  }
  print(paste0('Init slopes: ', init_slopes))
  print(paste0('Estimated slopes: ', slopes))
  return(slopes)
}

## Estimate empirical correlation
compute_group_corr <- function(data, pheno1, pheno2, labels){
  corrs <- data %>%
    # group_by(label) %>%
    mutate(label = factor(label, levels = labels)) %>%
    summarize(tidy(cor.test(get('BETA_1'), get('BETA_2')))) %$%
    estimate
  # names(corrs) <- labels
  return(corrs)
}

### ggplot version of plotting
prior.V <- function(scale = 0, slope = 1, cor = 0){

  # Returns prior covariance matrix for a line model defined by
  #  scale, slope and cor.
  # INPUT
  # scale > 0
  # slope in (-Inf,Inf]
  # cor >= 0, NOTE: negatively correlated effects are modeled by a negative slope,
  #                 not by a negative correlation

  if(cor > 1) stop("cor > 1")
  if(cor < 0) stop("cor < 0")
  if(scale < 0) stop("scale < 0")

  theta = atan(slope) - pi/4
  # R is rotation of angle theta
  R = matrix(c(cos(theta), -sin(theta),
               sin(theta), cos(theta)),
             ncol = 2, byrow = T)
  #Rotate diagonal case of correlation 'cor' by angle theta to S:
  S = R %*% matrix(c(1, cor, cor, 1), ncol = 2) %*% t(R)
  S = scale^2 * S / max(as.vector(S))
  return(S)
}


visualize_line_models <- function(scales, slopes, cors,
                                  model.names = NULL, model.cols = NULL,
                                  legend.position = "bottomright",
                                  xlim = NULL, ylim = NULL,
                                  xlab = "EFFECT1", ylab = "EFFECT2") {

  K <- length(slopes) # number of models
  if (length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if (length(cors) != K) stop("Length of 'scales' and 'cors' do not match.")
  if (any(cors > 1) | any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")

  lim <- 3 * max(scales) # Default: show models within 3 SDs
  if (is.null(xlim)) xlim <- c(-lim, lim)
  if (is.null(ylim)) ylim <- c(-lim, lim)
  if (is.null(model.names)) model.names <- paste0("M", 1:K)
  if (is.null(model.cols)) model.cols <- 1:K

  data_list <- list()

  for (ii in 1:K) {
    if (scales[ii] < 1e-16) {
      data_list[[ii]] <- data.frame(x = 0, y = 0, y.upper = 0, y.lower = 0, model = model.names[ii])
      next
    }
    if (!is.finite(slopes[ii])) {
      data_list[[ii]] <- data.frame(x = c(0, 0), y = ylim, y.upper = ylim, y.lower = ylim, model = model.names[ii])
    } else {
      x <- seq(xlim[1], xlim[2], length.out = 500)
      y <- slopes[ii] * x
      y.upper <- y
      y.lower <- y

      if (cors[ii] < 1) {
        Sigma <- prior.V(scale = scales[ii], slope = slopes[ii], cor = cors[ii])
        b <- qchisq(0.95, df = 2)
        a <- as.numeric(t(solve(Sigma)))
        y.upper <- (-(a[2] + a[3]) * x + sqrt((a[2] + a[3])^2 * x^2 - 4 * a[4] * (a[1] * x^2 - b))) / (2 * a[4])
        y.lower <- (-(a[2] + a[3]) * x - sqrt((a[2] + a[3])^2 * x^2 - 4 * a[4] * (a[1] * x^2 - b))) / (2 * a[4])
      }
      data_list[[ii]] <- data.frame(x = x, y = y, y.upper = y.upper, y.lower = y.lower, model = model.names[ii])
    }
  }

  plot_data <- do.call(rbind, data_list)

  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_line(aes(color=model)) +
    geom_ribbon(aes(ymin = y.lower, ymax = y.upper, fill = model, color = NULL), alpha = 0.1) +
    labs(x = xlab, y = ylab, color=NULL, fill=NULL) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    scale_fill_manual(breaks = model.names, values = model.cols) +
    scale_color_manual(breaks = model.names, values = model.cols)

  if (!is.null(legend.position)) {
    p <- p + theme(legend.position = legend.position)
  }

  print(p)
}

### SCENARIO [1]: PHENOTYPE(1) ONLY vs. PHENOTYPE(2) ONLY vs. BOTH
### SCENARIO [2]: PHENOTYPE(1)  vs. PHENOTYPE(2)
### SCENARIO [3]: PHENOTYPE(1) ONLY vs. PHENOTYPE(2) ONLY
run_analysis <- function(
    pheno1,
    pheno2,
    pop,
    scenario_idx,
    include_corr = FALSE,
    optimize_par = rbind(c(FALSE, FALSE, FALSE), # scale, slope, corr for model (1)
                         c(FALSE, TRUE, FALSE), # scale, slope, corr for model (2)
                         c(FALSE, FALSE, FALSE)), # scale, slope, corr for model (3)
    model.priors = c(1,1,1),
    scales = c(0.2, 0.2, 0.2),
    init_slopes = c(0.5, 0.5, 0.5),
    cors = c(0.99, 0.99, 0.99),
    colors = c("orange", "navy", 'red'),
    posterior_probs = c(0.95, 0.9, 0.8),
    root_path = '~/Downloads/',
    save = F
){
  r.lkhood <- 0
  label1 <- get_pheno_label(pheno1)
  label2 <- get_pheno_label(pheno2)
  print(paste0('--------', label1, ' vs. ', label2, '---------'))
  clumped_data <- load_clumped_data(pheno1=pheno1, pheno2=pheno2, pop=pop, root_path = root_path)

  if(scenario_idx == 1){
    labels <- c(label1, label2, 'Both')
  }else{
    model.priors <- model.priors[1:2]
    init_slopes <- init_slopes[1:2]
    scales <- scales[1:2]
    cors <- cors[1:2]
    colors <- colors[1:2]
    labels <- c(label1, label2)
    if(scenario_idx == 3){
      clumped_data <- relabel_data(clumped_data, pheno1, pheno2)
    }
    clumped_data <- clumped_data %>%
      mutate(label = factor(label, levels = labels))
  }

  names(colors) <- labels
  if(all(table(clumped_data$label)>0) & nrow(clumped_data) > 0){
    model.names <- labels
    names(init_slopes) <- labels
    print(init_slopes)
    slopes <- estimate_beta_from_em(
      data = clumped_data,
      beta1 = 'BETA_1',
      beta2 = 'BETA_2',
      init_slopes = init_slopes,
      labels = labels)

    cor_tag <- ''
    if(include_corr){
      cor_tag <- '_corr_included'
      r.lkhood <- compute_group_corr(clumped_data, pheno1, pheno2, labels)
    }

    opt_tag <- ''
    if(!is.null(optimize_par)){
      opt_tag <- '_par_optimized'
      par.include = optimize_par
      line.models.optimize(
        X = clumped_data[,c('BETA_1', 'BETA_2')],
        SE = clumped_data[,c('SE_1','SE_2')],
        par.include = par.include,
        init.scales = scales,
        init.slopes = slopes,
        init.cors = cors,
        model.priors = model.priors,
        model.names = model.names,
        r.lkhood = r.lkhood)
    }

    res = line.models(X = clumped_data[,c('BETA_1', 'BETA_2')],
                      SE = clumped_data[,c('SE_1','SE_2')],
                      scales = scales,
                      slopes = slopes,
                      cors = cors,
                      model.names = model.names,
                      model.priors = model.priors,
                      r.lkhood = r.lkhood,
    )
    clumped_data = cbind(clumped_data, res)
    clumped_data$cols = factor(labels[apply(res, 1, which.max)], levels = labels)
    # write_csv(clumped_data, paste0('~/Desktop/tmp_results/', pop, '_',label1,'_', if_else(pheno2 == 'spirometry', pheno2, label2), '_scenario_', scenario_idx, cor_tag, opt_tag, '.csv'))

    for(prob in posterior_probs){
      png_path = paste0('~/Desktop/tmp_figures/', pop, '_',label1,'_', if_else(pheno2 %in% c('spirometry', 'cigday'), pheno2, label2), '_prob_', prob, '_scenario_', scenario_idx, cor_tag, opt_tag, '.png')

      ylab = substitute(beta[sub], list(sub =  label2))
      if(pheno2 == 'spirometry') ylab = substitute(beta[sub], list(sub = substitute(FEV[sub]/FVC, list(sub = 1))))
      p <- visualize_line_models(
        scales=scales,
        slopes=slopes,
        cors=cors,
        model.names = model.names,
        model.cols = colors,
        legend.position = "bottomright",
        xlim = c(min(clumped_data[,'BETA_1']), max((clumped_data[,'BETA_1']))),
        ylim = c(min(clumped_data[,'BETA_2']), max((clumped_data[,'BETA_2']))),
        # xlab = expression(beta[1]~get(label1)),
        xlab =  substitute(beta[sub], list(sub = label1)),
        ylab = ylab
      )
      #color by max model
      ind = (apply(res, 1, max) > prob) #color only high probability cases
      p <- p +
        geom_point(data = clumped_data, aes(x = BETA_1, y = BETA_2, color=label), pch = 3, alpha = 0.15, size = 1.5) +
        geom_point(data = clumped_data[ind,], aes(x = BETA_1, y = BETA_2, color=cols), pch = 19, size = 1, alpha = 0.9) +
        themes
      if(!file.exists(png_path) | save){
        png(png_path,width=4, height=4, units = 'in', res = 300)
        print(p)
        dev.off()
      }
    }
    return(list(clumped_data, p))
  }
}

visualize_line_models_v2 <- function(scales, slopes, cors,
                                  model.names = NULL, model.cols = NULL,
                                  legend.position = "bottomright",
                                  xlim = NULL, ylim = NULL,
                                  xlab = "EFFECT1", ylab = "EFFECT2") {

  K <- length(slopes) # number of models
  if (length(scales) != K) stop("Length of 'scales' and 'slopes' do not match.")
  if (length(cors) != K) stop("Length of 'scales' and 'cors' do not match.")
  if (any(cors > 1) | any(cors < 0)) stop("Some value of 'cors' is outside [0,1].")

  lim <- 3 * max(scales) # Default: show models within 3 SDs
  if (is.null(xlim)) xlim <- c(-lim, lim)
  if (is.null(ylim)) ylim <- c(-lim, lim)
  if (is.null(model.names)) model.names <- paste0("M", 1:K)
  if (is.null(model.cols)) model.cols <- 1:K

  data_list <- list()

  for (ii in 1:K) {
    if (scales[ii] < 1e-16) {
      data_list[[ii]] <- data.frame(x = 0, y = 0, y.upper = 0, y.lower = 0, model = model.names[ii])
      next
    }
    if (!is.finite(slopes[ii])) {
      data_list[[ii]] <- data.frame(x = c(0, 0), y = ylim, y.upper = ylim, y.lower = ylim, model = model.names[ii])
    } else {
      x <- seq(xlim[1], xlim[2], length.out = 500)
      y <- slopes[ii] * x
      y.upper <- y
      y.lower <- y

      if (cors[ii] < 1) {
        Sigma <- prior.V(scale = scales[ii], slope = slopes[ii], cor = cors[ii])
        b <- qchisq(0.95, df = 2)
        a <- as.numeric(t(solve(Sigma)))
        y.upper <- (-(a[2] + a[3]) * x + sqrt((a[2] + a[3])^2 * x^2 - 4 * a[4] * (a[1] * x^2 - b))) / (2 * a[4])
        y.lower <- (-(a[2] + a[3]) * x - sqrt((a[2] + a[3])^2 * x^2 - 4 * a[4] * (a[1] * x^2 - b))) / (2 * a[4])
      }
      data_list[[ii]] <- data.frame(x = x, y = y, y.upper = y.upper, y.lower = y.lower, model = model.names[ii])
    }
  }

  plot_data <- do.call(rbind, data_list) %>%
    mutate(model = factor(model, levels = c(model.names, 'Both')))

  p <- ggplot(plot_data, aes(x = x, y = y)) +
    geom_line(aes(color=model)) +
    geom_ribbon(aes(ymin = y.lower, ymax = y.upper, fill = model, color = NULL), alpha = 0.1) +
    labs(x = xlab, y = ylab, color=NULL, fill=NULL) +
    coord_cartesian(ylim = ylim, xlim = xlim) +
    scale_fill_manual(labels = c(model.names, 'Both'),breaks = c(model.names, 'Both'), values = c(model.cols, 'Both' = 'red')) +
    scale_color_manual(labels = c(model.names, 'Both'),breaks = c(model.names, 'Both'), values = c(model.cols, 'Both' = 'red'))

  if (!is.null(legend.position)) {
    p <- p + theme(legend.position = legend.position)
  }

  print(p)
}

run_analysis_v2 <- function(
    pheno1,
    pheno2,
    pop,
    include_corr = FALSE,
    optimize_par = rbind(c(FALSE, FALSE, FALSE), # scale, slope, corr for model (1)
                         c(FALSE, TRUE, FALSE)), # scale, slope, corr for model (2)
    model.priors = c(1,1),
    scales = c(0.2, 0.2),
    init_slopes = c(0.5, 0.5),
    cors = c(0.99, 0.99),
    colors = c("orange", "navy", 'red'),
    posterior_probs = c(0.95, 0.9, 0.8),
    root_path = '~/Downloads/',
    save = F
){
  r.lkhood <- 0
  label1 <- get_pheno_label(pheno1)
  label2 <- get_pheno_label(pheno2)
  print(paste0('--------', label1, ' vs. ', label2, '---------'))
  clumped_data <- load_clumped_data(pheno1=pheno1, pheno2=pheno2, pop=pop, root_path = root_path)
  labels <- c(label1, label2, 'Both')

  names(colors) <- labels
  slope_data <- clumped_data %>% filter(label != 'Both') %>% mutate(label = factor(label, levels=c(label1, label2)))
  if(nrow(clumped_data) > 0 & all(table(slope_data$label)>0)){
    model.names <- labels[1:2]
    names(init_slopes) <- labels[1:2]
    print(init_slopes)
    print(paste0('N variant left for slope computation:', table(slope_data$label)))
    slopes <- estimate_beta_from_em(
      data = slope_data,
      beta1 = 'BETA_1',
      beta2 = 'BETA_2',
      init_slopes = init_slopes,
      labels = labels[1:2])

    cor_tag <- ''
    if(include_corr){
      cor_tag <- '_corr_included'
      r.lkhood <- compute_group_corr(clumped_data, pheno1, pheno2, labels)
    }

    opt_tag <- ''
    if(!is.null(optimize_par)){
      opt_tag <- '_par_optimized'
      par.include = optimize_par
      line.models.optimize(
        X = clumped_data[,c('BETA_1', 'BETA_2')],
        SE = clumped_data[,c('SE_1','SE_2')],
        par.include = par.include,
        init.scales = scales,
        init.slopes = slopes,
        init.cors = cors,
        model.priors = model.priors,
        model.names = model.names,
        r.lkhood = r.lkhood)
    }

    res = line.models(X = clumped_data[,c('BETA_1', 'BETA_2')],
                      SE = clumped_data[,c('SE_1','SE_2')],
                      scales = scales,
                      slopes = slopes,
                      cors = cors,
                      model.names = model.names,
                      model.priors = model.priors,
                      r.lkhood = r.lkhood,
    )
    clumped_data = cbind(clumped_data, res)
    clumped_data$cols = factor(labels[apply(res, 1, which.max)], levels = labels[1:2])
    result_output <- paste0('~/Desktop/tmp_results/', pop, '_',label1,'_', if_else(pheno2 %in% c('spirometry', 'cigday'), pheno2, label2), cor_tag, opt_tag, '_new.csv')
    if(!file.exists(result_output)){
      write_csv(clumped_data, result_output)
    }


    for(prob in posterior_probs){
      png_path = paste0('~/Desktop/tmp_figures/', pop, '_',label1,'_', if_else(pheno2 %in% c('spirometry', 'cigday'), pheno2, label2), '_prob_', prob , cor_tag, opt_tag, '.png')

      ylab = substitute(beta[sub], list(sub =  label2))
      if(pheno2 == 'spirometry') ylab = substitute(beta[sub], list(sub = substitute(FEV[sub]/FVC, list(sub = 1))))
      clumped_data <- clumped_data %>% mutate(label=factor(label, levels=labels),
                                              cols = factor(cols, levels=labels))
      p <- visualize_line_models_v2(
        scales=scales,
        slopes=slopes,
        cors=cors,
        model.names = model.names,
        model.cols = colors[1:2],
        legend.position = "bottomright",
        xlim = c(min(clumped_data[,'BETA_1']), max((clumped_data[,'BETA_1']))),
        ylim = c(min(clumped_data[,'BETA_2']), max((clumped_data[,'BETA_2']))),
        # xlab = expression(beta[1]~get(label1)),
        xlab =  substitute(beta[sub], list(sub = label1)),
        ylab = ylab
      )
      #color by max model
      ind = (apply(res, 1, max) > prob) #color only high probability cases
      p <- p +
        # scale_color_manual(values = colors, labels=labels) +
        # scale_fill_manual(values = colors, labels=labels) +
        geom_point(data = clumped_data, aes(x = BETA_1, y = BETA_2, color=label), pch = 3, alpha = 0.2, size = 2) +
        geom_point(data = clumped_data[ind,], aes(x = BETA_1, y = BETA_2, color=cols), pch = 19, size = 1) +
        themes +
        guides(
          fill = 'none',
        )
      if(!file.exists(png_path) | save){
        png(png_path,width=4, height=4, units = 'in', res = 300)
        print(p)
        dev.off()
      }
    }
    return(list(clumped_data, p))
  }
}


### Gene annotation
# install.packages("BiocManager")
# BiocManager::install("rtracklayer")
# BiocManager::install("GenomicRanges")
library(rtracklayer)
library(GenomicRanges)


annotate_genes <- function(data, ref_genome='GRCh37'){
  if(ref_genome == 'GRCh37'){
    gtf_path <- '~/Dropbox (Partners HealthCare)/analysis/lung_function_pleiotropy/paper/Homo_sapiens.GRCh37.75.gtf'
  }else{
    gtf_path <- '~/Dropbox (Partners HealthCare)/analysis/lung_function_pleiotropy/paper/Homo_sapiens.GRCh38.105.gtf'
  }

  gtf_data <- import(gtf_path, format = "gtf")
  # Extract gene information
  gtf_genes <- gtf_data[gtf_data$type == "gene"]
  gtf_gene_info <- data.frame(
    seqnames = seqnames(gtf_genes),
    start = start(gtf_genes),
    end = end(gtf_genes),
    strand = strand(gtf_genes),
    gene_id = mcols(gtf_genes)$gene_id,
    gene_name = mcols(gtf_genes)$gene_name
  )
  # Convert SNP data to a GRanges object
  snps_gr <- GRanges(
    seqnames = data$CHR,
    ranges = IRanges(start = data$POS, end = data$POS)
  )

  # Find the nearest genes for each SNP
  nearest_genes <- nearest(snps_gr, gtf_genes)

  # Extract the nearest gene information
  nearest_gene_info <- gtf_gene_info[nearest_genes,]

  # Combine the SNPs with the nearest gene information
  snps_with_genes <- cbind(data, nearest_gene_info)
  return(snps_with_genes)
}

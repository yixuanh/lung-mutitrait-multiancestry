source('~/Dropbox (Partners HealthCare)/github_repo/random_scripts/lung_multitrait_linemodel_utils.R')

### Update raw linemodel figure
root_path = '~/Dropbox (Partners HealthCare)/analysis/lung_function_pleiotropy/comparison_data/'
pheno1s <- c(rep('asthma', 7), rep('copd', 6), rep('lung_cancer', 5))
pheno2s <- c('copd', 'lung_cancer', 'eosinophil', 'neutrophil', 'smkinit', 'spirometry', 'cigday',
             'lung_cancer', 'eosinophil', 'neutrophil', 'smkinit', 'spirometry', 'cigday',
             'eosinophil', 'neutrophil', 'smkinit', 'spirometry', 'cigday')
pops <- c('meta', 'afr', 'amr', 'eas', 'eur')
for(pop in pops){
  for(i in 1:length(pheno1s)){
    pheno1 = pheno1s[i]
    pheno2 = pheno2s[i]
    label1 <- get_pheno_label(pheno1)
    label2 <- get_pheno_label(pheno2)
    if(pop == 'amr' & (pheno1 == 'lung_cancer' | pheno2 == 'lung_cancer')) next
    if(pop == 'meta' & (pheno2 %in% c('eosinophil', 'neutrophil', 'smkinit', 'cigday'))) next
    png_path =  paste0('~/Desktop/tmp_figures/', pop, '_',label1,'_', if_else(pheno2 %in% c('spirometry', 'cigday'), pheno2, label2), '_prob_0.99.png')
    if(file.exists(png_path)) next
    run_analysis_v2(
      pheno1=pheno1,
      pheno2=pheno2,
      pop = pop,
      include_corr = FALSE,
      optimize_par = NULL,
      model.priors = c(1,1),
      scales = c(0.2, 0.2),
      init_slopes = c(0, 0.5),
      cors = c(0.99, 0.99),
      colors = c("orange", "navy", 'red'),
      posterior_prob = c(0.99),
      root_path = root_path
    )
  }
}

### Main Figure 3b and 3c
root_path = '~/Dropbox (Partners HealthCare)/analysis/lung_function_pleiotropy/comparison_data/'
results1 <- run_analysis_v2(
  pheno1='copd',
  pheno2='cigday',
  pop = 'eur',
  include_corr = FALSE,
  optimize_par = NULL,
  model.priors = c(1,1),
  scales = c(0.2, 0.2),
  init_slopes = c(0, 0.5),
  cors = c(0.99, 0.99),
  colors = c("orange", "navy", 'red'),
  posterior_prob = c(0.99),
  root_path = root_path
)

table1 <- results1[[1]] %>%
  filter(COPD > 0.99 | `Cig/Day` > 0.99) %>%
  mutate(CHR = as.numeric(str_split(locus, ':') %>% map_chr(., 1)),
         POS = as.numeric(str_split(locus, ':') %>% map_chr(., 2)))
table1 <- annotate_genes(table1)
table1 <- table1 %>%
  select(CHR, POS, Alleles=alleles, rsID = rsid_1, Phenotype_1 = pheno_1, Population_1 = pop_1, BETA_1, P_1, SE_1, AF_1,N_case_1, N_ctrl_1,
         Phenotype_2 = pheno_2, Population_2 = pop_2, BETA_2, P_2, SE_2, AF_2, N_2, Significant_traits = label, COPD, `Cig/Day`,
         Confidently_associated = cols, nearest_gene_id = gene_id, nearest_gene_name = gene_name)
write_csv(table1, '~/Desktop/supp_table_7.csv')
table1 <- table1 %>%
  filter(COPD > 0.99) %>%
  group_by(gene_name) %>%
  mutate(max = max(COPD)) %>%
  filter(COPD == max)

results2 <- run_analysis_v2(
  pheno1='lung_cancer',
  pheno2='cigday',
  pop = 'eur',
  include_corr = FALSE,
  optimize_par = NULL,
  model.priors = c(1,1),
  scales = c(0.2, 0.2),
  init_slopes = c(0, 0.5),
  cors = c(0.99, 0.99),
  colors = c("orange", "navy", 'red'),
  posterior_prob = c(0.99),
  root_path = root_path
)

table2 <- results2[[1]] %>%
  filter(`Lung Cancer` > 0.99 | `Cig/Day` > 0.99) %>%
  mutate(CHR = as.numeric(str_split(locus, ':') %>% map_chr(., 1)),
         POS = as.numeric(str_split(locus, ':') %>% map_chr(., 2)),)
table2 <- annotate_genes(table2)
table2 <- table2 %>%
  select(CHR, POS, Alleles=alleles, rsID = RS_1, Phenotype_1 = pheno_1, Population_1 = pop_1, BETA_1, P_1, SE_1, AF_1,
         Phenotype_2 = pheno_2, Population_2 = pop_2, BETA_2, P_2, SE_2, N_2, Significant_traits = label, `Lung Cancer`, `Cig/Day`,
         Confidently_associated = cols, nearest_gene_id = gene_id, nearest_gene_name = gene_name)
write_csv(table2, '~/Desktop/supp_table_8.csv')
table2 <- table2 %>%
  filter(`Lung Cancer` > 0.99) %>%
  group_by(gene_name) %>%
  mutate(max = max((`Lung Cancer`))) %>%
  filter(`Lung Cancer` == max)


p <- ggpubr::ggarrange(results1[[2]] +
                         geom_text_repel(data = table1, aes(x = BETA_1, y = BETA_2, color = cols, label = gene_name), size = 2, vjust = 2, hjust =1, max.overlaps = 100) +
                         theme(plot.margin = unit(c(0.5,0,0,0), "cm")) ,
                       results2[[2]]+
                         geom_text_repel(data = table2, aes(x = BETA_1, y = BETA_2, color = cols, label = gene_name), size = 2, vjust = 2, hjust =1, max.overlaps = 100) +
                         theme(plot.margin = unit(c(0.5,0,0,0), "cm")),
                       ncol=2, widths = c(1, 1),
                       # labels = c(paste('(A) COPD (\u03b21) vs. Cigday (\u03b22) - EUR'), paste('(B) Lung cancer (\u03b21) vs. Cigday (\u03b22) -EUR')),
                       labels = NULL,
                       font.label = list(size = 9, color = "gray45", face = "bold", family = NULL), hjust =0, vjust = 2)

png(paste0('~/Desktop/main_figure_3_prob_0.99.png'), width=8, height=4, units = 'in', res = 300)
print(p)
dev.off()


empty_p <- ggplot(mtcars, aes()) + geom_blank() + theme(axis.line.x = element_blank(), axis.line.y = element_blank())
png_path <- function(pop, pheno1, pheno2){
  return(paste0('~/Desktop/tmp_figures/', pop, '_', pheno1, '_', pheno2, '_prob_0.99.png'))
}


## Supp figure 6 - Asthma (Two groups)
pheno1 <-'asthma'
pheno2s <- c('copd','lung_cancer', 'smkinit', 'cigday', 'spirometry')
j = 1
p_lst <- list()
for(pop in c('afr', 'eas', 'eur', 'meta')){
  for(i in 1:length(pheno2s)){
    print(pop)
    tmp_get_pheno_label <- function(pheno){
      label <- if_else(pheno == 'spirometry', "FEV[1]/FVC", str_to_sentence(pheno))
      label <- if_else(label == 'Copd', toupper(label), label)
      return(label)
    }
    pheno1 <- tmp_get_pheno_label(pheno1)
    pheno2 <- if_else(pheno2s[i] != 'spirometry', tmp_get_pheno_label(pheno2s[i]), 'spirometry')
    print(pheno2)
    raw_png_path <- png_path(pop=pop, pheno1=if_else(pheno1 == 'Lung_cancer', 'Lung Cancer', pheno1), pheno2=if_else(pheno2 == 'Lung_cancer', 'Lung Cancer', pheno2))
    print(raw_png_path)
    if(!file.exists(raw_png_path)){
      tmp_p <- empty_p
    }else{
      tmp_p <- run_analysis_v2(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1),
        scales = c(0.2, 0.2),
        init_slopes = c(0, 0.5),
        cors = c(0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.99),
        root_path = root_path
      )[[2]]
    }
    print(tmp_p)
    p_lst[[j]] <- tmp_p +
      theme(plot.margin = unit(c(0.5,0,0,0), "cm"),
            legend.position = 'None')
    j = j+1
  }
}

p <- egg::ggarrange(plots = p_lst[1:20],
                    ncol=5, widths = c(rep(1, 5)),
                    # labels = c(paste('(A) COPD (\u03b21) vs. Cigday (\u03b22) - EUR'), paste('(B) Lung cancer (\u03b21) vs. Cigday (\u03b22) -EUR')),
                    labels = c('a','', '', 'b', '', letters[3:6], '', letters[7:13], '', '', letters[14]),
                    font.label = list(size = 2, face = "bold"), hjust =0, vjust = 2)
top_label <-  paste('     ', sapply(pheno2s, function(x) get_pheno_label(x)), collapse = '               ')
p <- annotate_figure(p,
                     top =text_grob(bquote(~ ~ ~ ~ ~ ~ ~COPD~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Lung~ ~Cancer~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Smkinit~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Cig/Day~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~FEV[1]/FVC)),
                     left = paste(toupper(rev(c('afr', 'eas', 'eur', 'meta'))), collapse = '                         '))

png(paste0('~/Desktop/supp_figure_6_', pheno1,'_two_groups.png'), width=8, height=7, units = 'in', res = 300)
print(p)
dev.off()


#### Supp figure S7 - COPD (Two groups)
pheno1 <-'copd'
pheno2s <- c('lung_cancer', 'smkinit', 'cigday', 'spirometry')
j = 1
p_lst <- list()
for(pop in c('afr', 'eas', 'eur', 'meta')){
  for(i in 1:length(pheno2s)){
    print(pop)
    tmp_get_pheno_label <- function(pheno){
      label <- if_else(pheno == 'spirometry', "FEV[1]/FVC", str_to_sentence(pheno))
      label <- if_else(label == 'Copd', toupper(label), label)
      return(label)
    }
    pheno1 <- tmp_get_pheno_label(pheno1)
    pheno2 <- if_else(pheno2s[i] != 'spirometry', tmp_get_pheno_label(pheno2s[i]), 'spirometry')
    print(pheno2)
    raw_png_path <- png_path(pop=pop, pheno1=if_else(pheno1 == 'Lung_cancer', 'Lung Cancer', pheno1), pheno2=if_else(pheno2 == 'Lung_cancer', 'Lung Cancer', pheno2))
    print(raw_png_path)
    if(!file.exists(raw_png_path)){
      tmp_p <- empty_p
    }else{
      tmp_p <- run_analysis_v2(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1),
        scales = c(0.2, 0.2),
        init_slopes = c(0, 0.5),
        cors = c(0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.99),
        root_path = root_path
      )[[2]]
    }
    print(tmp_p)
    p_lst[[j]] <- tmp_p +
      theme(plot.margin = unit(c(0.5,0,0,0), "cm"),
            legend.position = 'None')
    j = j+1
  }
}

p <- egg::ggarrange(plots = p_lst[1:16],
                    ncol=4, widths = c(rep(1, 4)),
                    # labels = c(paste('(A) COPD (\u03b21) vs. Cigday (\u03b22) - EUR'), paste('(B) Lung cancer (\u03b21) vs. Cigday (\u03b22) -EUR')),
                    labels = c('a','', 'b', '', letters[3:5], '', letters[6:10], '', '', letters[11]),
                    # labels = letters[1:16],
                    font.label = list(size = 2, face = "bold"), hjust =0, vjust = 2)
top_label <- paste('', sapply(pheno2s, get_pheno_label), collapse = '                       ')
p <- annotate_figure(p,
                     top = text_grob(bquote(~ ~Lung~ ~Cancer~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Smkinit~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Cig/Day~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~FEV[1]/FVC)),
                     left = paste(toupper(rev(c('afr', 'eas', 'eur', 'meta'))), collapse = '                         '))

png(paste0('~/Desktop/supp_figure_', pheno1,'_two_groups.png'), width=8, height=7, units = 'in', res = 300)
print(p)
dev.off()

#### Supp figure S8 - Lung cancer (Two groups)
pheno1 <-'lung_cancer'
pheno2s <- c('smkinit', 'cigday', 'spirometry')
j = 1
p_lst <- list()
for(pop in c('afr', 'eas', 'eur', 'meta')){
  for(i in 1:length(pheno2s)){
    print(pop)
    tmp_get_pheno_label <- function(pheno){
      label <- if_else(pheno == 'spirometry', "FEV[1]/FVC", str_to_sentence(pheno))
      label <- if_else(label == 'Copd', toupper(label), label)
      return(label)
    }
    pheno1 <- tmp_get_pheno_label(pheno1)
    pheno2 <- if_else(pheno2s[i] != 'spirometry', tmp_get_pheno_label(pheno2s[i]), 'spirometry')
    print(pheno2)
    raw_png_path <- png_path(pop=pop, pheno1=if_else(pheno1 == 'Lung_cancer', 'Lung Cancer', pheno1), pheno2=if_else(pheno2 == 'Lung_cancer', 'Lung Cancer', pheno2))
    print(raw_png_path)
    if(!file.exists(raw_png_path)){
      tmp_p <- empty_p
    }else{
      tmp_p <- run_analysis_v2(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1),
        scales = c(0.2, 0.2),
        init_slopes = c(0, 0.5),
        cors = c(0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.99),
        root_path = root_path
      )[[2]]
    }
    print(tmp_p)
    p_lst[[j]] <- tmp_p +
      theme(plot.margin = unit(c(0.5,0,0,0), "cm"),
            legend.position = 'None')
    j = j+1
  }
}

p <- egg::ggarrange(plots = p_lst[4:12],
                    ncol=3, widths = c(rep(1, 3)),
                    # labels = c(paste('(A) COPD (\u03b21) vs. Cigday (\u03b22) - EUR'), paste('(B) Lung cancer (\u03b21) vs. Cigday (\u03b22) -EUR')),
                    labels = c('a', 'b', '', 'c', 'd', 'e', '', '', 'f'),
                    # labels = letters[1:12],
                    font.label = list(size = 2, face = "bold"), hjust =0, vjust = 2)
top_label <- paste('', sapply(pheno2s, get_pheno_label), collapse = '                       ')
p <- annotate_figure(p,
                     top = text_grob(bquote( ~Smkinit~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Cig/Day~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~FEV[1]/FVC)),
                     left = paste(toupper(rev(c('eas', 'eur', 'meta'))), collapse = '                         '))

png(paste0('~/Desktop/supp_figure_8_', pheno1,'_two_groups.png'), width=6, height=6, units = 'in', res = 300)
print(p)
dev.off()

#### Supp figure S9 - Asthma (Three groups)
# pheno1 <-'asthma'
# pheno2s <- c('copd','lung_cancer', 'smkinit', 'cigday', 'spirometry')
# pheno1 <-'copd'
# pheno2s <- c('lung_cancer', 'smkinit', 'cigday', 'spirometry')
pheno1 <-'lung_cancer'
pheno2s <- c('smkinit', 'cigday', 'spirometry')
j = 1
p_lst <- list()
for(pop in c('afr', 'eas', 'eur', 'meta')){
  for(i in 1:length(pheno2s)){
    print(pop)
    tmp_get_pheno_label <- function(pheno){
      label <- if_else(pheno == 'spirometry', "FEV[1]/FVC", str_to_sentence(pheno))
      label <- if_else(label == 'Copd', toupper(label), label)
      return(label)
    }
    pheno1 <- tmp_get_pheno_label(pheno1)
    pheno2 <- if_else(pheno2s[i] != 'spirometry', tmp_get_pheno_label(pheno2s[i]), 'spirometry')
    print(pheno2)
    png_path1 <- png_path(pop=pop, pheno1=pheno1, pheno2=pheno2, label=1)
    png_path2 <- png_path(pop=pop, pheno1=pheno1, pheno2=pheno2, label=2)
    print(png_path2)
    if(!file.exists(png_path1) & !file.exists(png_path2)){
      tmp_p <- empty_p
    }else if(file.exists(png_path1)){
      tmp_p <- run_analysis(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        scenario_idx=1,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1,1),
        scales = c(0.2, 0.2, 0.2),
        init_slopes = c(0, 0.5, 0.5),
        cors = c(0.99, 0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.95),
        root_path = root_path
      )[[2]]
    }else{

      tmp_p <- run_analysis(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        scenario_idx=2,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1,1),
        scales = c(0.2, 0.2, 0.2),
        init_slopes = c(0, 0.5, 0.5),
        cors = c(0.99, 0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.95),
        root_path = root_path
      )[[2]]
    }
    print(class(tmp_p))
    p_lst[[j]] <- tmp_p +
      theme(plot.margin = unit(c(0.5,0,0,0), "cm"),
            legend.position = 'None')
    j = j+1
  }
}


p <- egg::ggarrange(plots = p_lst[1:20],
                    ncol=5, widths = c(rep(1, 5)),
                    # labels = c(paste('(A) COPD (\u03b21) vs. Cigday (\u03b22) - EUR'), paste('(B) Lung cancer (\u03b21) vs. Cigday (\u03b22) -EUR')),
                    labels = c('a','', '', 'b', '', letters[3:6], '', letters[7:13], '', '', letters[14]),
                    font.label = list(size = 2, face = "bold"), hjust =0, vjust = 2)
top_label <-  paste('     ', sapply(pheno2s, function(x) get_pheno_label(x)), collapse = '               ')
p <- annotate_figure(p,
                     top =text_grob(bquote(~ ~ ~ ~ ~ ~ ~COPD~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Lung~ ~Cancer~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Smkinit~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Cig/Day~  ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~FEV[1]/FVC)),
                     left = paste(toupper(rev(c('afr', 'eas', 'eur', 'meta'))), collapse = '                         '))

png(paste0('~/Desktop/supp_figure_9_', pheno1,'_three_groups.png'), width=8, height=7, units = 'in', res = 300)
print(p)
dev.off()

#### Supp figure S10 - COPD (Three group)
pheno1 <-'copd'
pheno2s <- c('lung_cancer', 'smkinit', 'cigday', 'spirometry')
j = 1
p_lst <- list()
for(pop in c('afr', 'eas', 'eur', 'meta')){
  for(i in 1:length(pheno2s)){
    print(pop)
    tmp_get_pheno_label <- function(pheno){
      label <- if_else(pheno == 'spirometry', "FEV[1]/FVC", str_to_sentence(pheno))
      label <- if_else(label == 'Copd', toupper(label), label)
      return(label)
    }
    pheno1 <- tmp_get_pheno_label(pheno1)
    pheno2 <- if_else(pheno2s[i] != 'spirometry', tmp_get_pheno_label(pheno2s[i]), 'spirometry')
    print(pheno2)
    png_path1 <- png_path(pop=pop, pheno1=pheno1, pheno2=pheno2, label=1)
    png_path2 <- png_path(pop=pop, pheno1=pheno1, pheno2=pheno2, label=2)
    print(png_path2)
    if(!file.exists(png_path1) & !file.exists(png_path2)){
      tmp_p <- empty_p
    }else if(file.exists(png_path1)){
      tmp_p <- run_analysis(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        scenario_idx=1,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1,1),
        scales = c(0.2, 0.2, 0.2),
        init_slopes = c(0, 0.5, 0.5),
        cors = c(0.99, 0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.95),
        root_path = root_path
      )[[2]]
    }else{

      tmp_p <- run_analysis(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        scenario_idx=2,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1,1),
        scales = c(0.2, 0.2, 0.2),
        init_slopes = c(0, 0.5, 0.5),
        cors = c(0.99, 0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.95),
        root_path = root_path
      )[[2]]
    }
    print(class(tmp_p))
    p_lst[[j]] <- tmp_p +
      theme(plot.margin = unit(c(0.5,0,0,0), "cm"),
            legend.position = 'None')
    j = j+1
  }
}


p <- egg::ggarrange(plots = p_lst[1:16],
                    ncol=4, widths = c(rep(1, 4)),
                    # labels = c(paste('(A) COPD (\u03b21) vs. Cigday (\u03b22) - EUR'), paste('(B) Lung cancer (\u03b21) vs. Cigday (\u03b22) -EUR')),
                    labels = c('a','', 'b', '', letters[3:5], '', letters[6:10], '', '', letters[11]),
                    # labels = letters[1:16],
                    font.label = list(size = 2, face = "bold"), hjust =0, vjust = 2)
top_label <- paste('', sapply(pheno2s, get_pheno_label), collapse = '                       ')
p <- annotate_figure(p,
                     top = text_grob(bquote(~ ~Lung~ ~Cancer~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Smkinit~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Cig/Day~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~FEV[1]/FVC)),
                     left = paste(toupper(rev(c('afr', 'eas', 'eur', 'meta'))), collapse = '                         '))

png(paste0('~/Desktop/supp_figure_10_', pheno1,'_three_groups.png'), width=8, height=7, units = 'in', res = 300)
print(p)
dev.off()

#### Supp figure S11 - Lung cancer (Three groups)
pheno1 <-'lung_cancer'
pheno2s <- c('smkinit', 'cigday', 'spirometry')
j = 1
p_lst <- list()
for(pop in c('afr', 'eas', 'eur', 'meta')){
  for(i in 1:length(pheno2s)){
    print(pop)
    tmp_get_pheno_label <- function(pheno){
      label <- if_else(pheno == 'spirometry', "FEV[1]/FVC", str_to_sentence(pheno))
      label <- if_else(label == 'Copd', toupper(label), label)
      return(label)
    }
    pheno1 <- tmp_get_pheno_label(pheno1)
    pheno2 <- if_else(pheno2s[i] != 'spirometry', tmp_get_pheno_label(pheno2s[i]), 'spirometry')
    print(pheno2)
    png_path1 <- png_path(pop=pop, pheno1=pheno1, pheno2=pheno2, label=1)
    png_path2 <- png_path(pop=pop, pheno1=pheno1, pheno2=pheno2, label=2)
    print(png_path2)
    if(!file.exists(png_path1) & !file.exists(png_path2)){
      tmp_p <- empty_p
    }else if(file.exists(png_path1)){
      tmp_p <- run_analysis(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        scenario_idx=1,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1,1),
        scales = c(0.2, 0.2, 0.2),
        init_slopes = c(0, 0.5, 0.5),
        cors = c(0.99, 0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.95),
        root_path = root_path
      )[[2]]
    }else{

      tmp_p <- run_analysis(
        pheno1=pheno1,
        pheno2=pheno2s[i],
        pop = pop,
        scenario_idx=2,
        include_corr = FALSE,
        optimize_par = NULL,
        model.priors = c(1,1,1),
        scales = c(0.2, 0.2, 0.2),
        init_slopes = c(0, 0.5, 0.5),
        cors = c(0.99, 0.99, 0.99),
        colors = c("orange", "navy", 'red'),
        posterior_prob = c(0.95),
        root_path = root_path
      )[[2]]
    }
    print(class(tmp_p))
    p_lst[[j]] <- tmp_p +
      theme(plot.margin = unit(c(0.5,0,0,0), "cm"),
            legend.position = 'None')
    j = j+1
  }
}

p <- egg::ggarrange(plots = p_lst[4:12],
                    ncol=3, widths = c(rep(1, 3)),
                    # labels = c(paste('(A) COPD (\u03b21) vs. Cigday (\u03b22) - EUR'), paste('(B) Lung cancer (\u03b21) vs. Cigday (\u03b22) -EUR')),
                    labels = c('a', 'b', '', '', 'c', 'd', '', '', 'e'),
                    # labels = letters[1:12],
                    font.label = list(size = 2, face = "bold"), hjust =0, vjust = 2)
top_label <- paste('', sapply(pheno2s, get_pheno_label), collapse = '                       ')
p <- annotate_figure(p,
                     top = text_grob(bquote( ~Smkinit~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~Cig/Day~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~FEV[1]/FVC)),
                     left = paste(toupper(rev(c('eas', 'eur', 'meta'))), collapse = '                         '))

png(paste0('~/Desktop/supp_figure_11_', pheno1,'_three_group.png'), width=6, height=6, units = 'in', res = 300)
print(p)
dev.off()


#### Supp figure S12 - S13
get_snp_lst <- function(pop, pheno1='Lung_cancer', pheno2='Cigday'){
  data <- read_csv(paste0('~/Desktop/tmp_results/', pop,'_', pheno1,'_', pheno2, '_new.csv'), col_types = cols(locus = col_character())) %>%
    filter(get(pheno1) > 0.99) %>%
    mutate(x  = paste0(locus, ':', alleles)) %$%
    x %>%
    unlist(.)
  return(data)
}

get_snp_lst_df <- function(pop, pheno1='Lung_cancer', pheno2='Cigday'){
  data <- read_csv(paste0('~/Desktop/tmp_results/', pop,'_', pheno1,'_', pheno2, '_new.csv'), col_types = cols(locus = col_character())) %>%
    filter(get(pheno1) > 0.99)%>%
    select(locus, alleles, rsid=RSID_2)
  data[, paste0('in_', toupper(pop))] <- T
  return(data)
}
pheno1s <- c('copd','lung_cancer')
pheno2s <- c('cigday', 'cigday')

library(VennDiagram)
for(i in 1:length(pheno1s)){
  label1 <- get_pheno_label(pheno1s[i])
  label2 <- get_pheno_label(pheno2s[i])
  pheno_label1 = label1
  pheno_label2 = if_else(pheno2s[i] %in% c('spirometry', 'cigday'), pheno2s[i], label2)
  print(label1)
  print(label2)
  comparison_data <- merge(get_snp_lst_df('eas', pheno1 = pheno_label1, pheno2 = pheno_label2),
                           get_snp_lst_df('eur', pheno1 = pheno_label1, pheno2 = pheno_label2),
                           by = c('locus', 'alleles'), all=T)
  comparison_data <- merge(comparison_data, get_snp_lst_df('meta', pheno1 = pheno_label1, pheno2 = pheno_label2), by = c('locus', 'alleles'), all=T)
  comparison_data <- comparison_data %>%
    mutate(in_EAS  = if_else(is.na(in_EAS), FALSE, in_EAS),
           in_EUR  = if_else(is.na(in_EUR), FALSE, in_EUR),
           in_META  = if_else(is.na(in_META), FALSE, in_META),
    )
  write_csv(comparison_data, paste0('~/Desktop/tmp_results/', pheno_label1, '_', pheno_label2, '_new.csv'))
  # Chart
  venn.diagram(
    x = list(get_snp_lst('eas', pheno1 = pheno_label1, pheno2 = pheno_label2),
             get_snp_lst('eur', pheno1 = pheno_label1, pheno2 = pheno_label2),
             get_snp_lst('meta', pheno1 = pheno_label1, pheno2 = pheno_label2)),
    category.names = c("EAS" , "EUR" , "META"),
    filename = paste0('~/Desktop/tmp_figures/', pheno_label1, '_', pheno_label2, '_venn_diagramm.png'),
    output=T,
    # Output features
    imagetype="png" ,
    # height = 480 ,
    # width = 480 ,
    resolution = 300,
    compression = "lzw",
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c(color_eas, color_eur, color_oth),
    # Numbers
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",

    # Set names
    cat.cex = 2,
    cat.fontface = "bold",
    # cat.default.pos = "outer",
    cat.pos = c(0, 0, 0),
    cat.dist = c(0.025, 0.025, 0.005),
    cat.fontfamily = "sans",
    cat.col = c(color_eas, color_eur, color_oth),
    rotation = 1
  )
}


venn_table_copd <- read_csv('~/Desktop/tmp_results/COPD_cigday_new.csv', col_types = cols(locus = col_character()))
venn_table_copd <- venn_table_copd %>%
  mutate(CHR = as.numeric(str_split(locus, ':') %>% map_chr(., 1)),
         POS = as.numeric(str_split(locus, ':') %>% map_chr(., 2)))
venn_table_copd <- annotate_genes(venn_table_copd)
write_csv(venn_table_copd, '~/Desktop/tmp_results/COPD_cigday_new.csv')

venn_table_lung <- read_csv('~/Desktop/tmp_results/Lung Cancer_cigday_new.csv', col_types = cols(locus = col_character()))
venn_table_lung <- venn_table_lung %>%
  mutate(CHR = as.numeric(str_split(locus, ':') %>% map_chr(., 1)),
         POS = as.numeric(str_split(locus, ':') %>% map_chr(., 2)))
venn_table_lung <- annotate_genes(venn_table_lung)
write_csv(venn_table_lung, '~/Desktop/tmp_results/Lung Cancer_cigday_new.csv')

#### Check shared variants confidently identified in both Lung cancer vs. Cig/Day and COPD vs. Cig/Day
table1 <- read_csv('~/Desktop/supp_table_7.csv') %>% select(CHR, POS, Alleles, nearest_gene_name, COPD, `Cig/Day`)
table2 <- read_csv('~/Desktop/supp_table_8.csv') %>% select(CHR, POS, Alleles, nearest_gene_name, `Lung Cancer`, `Cig/Day`)
table <- merge(table1, table2, by = c('CHR', 'POS', 'Alleles'), all=T)

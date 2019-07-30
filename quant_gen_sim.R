#!/usr/bin/Rscript

# visualizations of simple quantitative / evolutionary genetics simulations
# (for teaching)
# Max Press July 28, 2019

sim_trait = function(num_loci = 10, h2 = 0.8, gen_model = "additive", num_inds = 100, ploidy = 1, 
                     rand_seed=NULL, allele_freq = 0.5, effect_size_dist = "normal", plot_method="circles") {
  # num_loci: number of loci contributing to trait
  # rand_seed: a seed in case you want to repeat the computation
  # h2: heritability, e.g. the proportion of genetic variance out of total variance of trait
  # gen_model: how loci contribute to the trait. currently only supports additive.
  # num_inds: size of population. currently only supports 20, 50, 100
  # ploidy: how many copies it's possible to have of each allele.
  # TODO: should have params controlling:
  # -allele frequency distribution (NOT JUST AVERAGE FREQUENCY OF "1" ALLELE)
  # -effect size distribution (assumes normal)
  # -simulate/handle loci NOT contributing to the trait
  if (!(is.null(rand_seed))) {
    set.seed(rand_seed)
  }
  
  if (length(allele_freq) == 1) {
    # would be handled differently (more realistically) with different allele frequencies across loci
    genotypes = rbinom(n = (num_loci * num_inds), size = ploidy, prob = allele_freq)
  } else if (length(allele_freq) == num_loci) {
    
  } else {
    stop(paste("length of allele_freq parameter must be either 1 or equal to the number of loci.", 
               "Instead it is of length", length(allele_freq)))
  }
  
  # individuals are rows and genotypes are columns
  geno_matrix = matrix(data = genotypes, ncol = num_loci)
  
  ### if split out, geno sim method returns here.
  
  ### phenotype compute method goes here
  # sim effect sizes- a little rough. would be good to nail down distributions a little more.
  effect_sizes = sim_effect_sizes(num_loci = num_loci, effect_size_dist = effect_size_dist, mean_sd_ratio = 10)
  genetic_pheno = (geno_matrix %*% effect_sizes)[,]
  print(length(genetic_pheno))
  V_g = var(genetic_pheno)
  V_e = ((1 - h2) / h2) * V_g
  print(V_g)
  print(V_e)
  # might be nice to tune env effect distribution too
  env_dist = rnorm(n = num_inds, sd = sqrt(V_e))
  #env_var = env_dist * V_e
  phenotypes = genetic_pheno + env_dist
  
  # plotting method
  if (plot_method == "circles") {
    circle_plots(phenotypes)
  } else {
    hist(phenotypes, 20, xlab = "Phenotype value")
  }
  
}

sim_effect_sizes = function(num_loci, effect_size_dist = "normal", mean_sd_ratio = 5, constant=4) {
  # mean_sd_ratio is only as stated for normal effect size dist. should tailor more to handle edge cases,
  # e.g. negative values.
  # for exp rate is inverse of mean_sd_ratio
  if (effect_size_dist == "normal") {
    effect_sizes = rnorm(num_loci, mean = mean_sd_ratio)
  } else if (effect_size_dist == "exponential") {
    effect_sizes = rexp(n = num_loci, rate = (1 / mean_sd_ratio))
  }
  return(effect_sizes + constant) #
}

circle_plots = function(data) {
  root = sqrt(length(data))
  x_axes = c()
  y_axes = c()
  for (i in 1:root) {
    x_axes = append(x_axes, 1:root)
    y_axes = append(y_axes, rep(i, times = root))
  }
  data[data < 0] = 0
  #print(data)
  scalar = mean(data) / 2
  plot(x_axes, y_axes, cex = data / scalar, 
       pch = 19, ylim = c(0.5, root+.5), xlim = c(0.5, root+.5),
       xlab = "", ylab = "", main = "Population phenotypes")
  grid = (0:root)+.5
  abline(h = grid, v = grid)
}
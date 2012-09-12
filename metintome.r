# A variety of methods for testing (and visualizing) interactions
# between species in a microbiome.

library(gplots)
library(entropy)
library(stringr)


namespecies = function(x, replace = F){
  if (! is.null(names(x)) && ! replace){
    # Do nothing
    return(x)
  }
  else {
    numbers = 1:length(x)
    len = floor(log10(length(x))) + 1
    num_strings = str_pad(numbers, len, "left", pad = 0)
    names(x) = paste("SP", num_strings, sep = '')
    return(x)
  }
}

# rxs stands for replicate-by-species (matrix)
rabund = function(rxs_mat){
  col_totals = apply(rxs_mat, c(2), sum)
  grand_total = sum(col_totals)
  rabunds = col_totals / grand_total
  rabunds = namespecies(rabunds)
  return(rabunds)
}

sim_neut_1rep = function(sample_size, rabunds){
  rabunds = namespecies(rabunds)
  species = names(rabunds)
  out = array(0, dim = length(rabunds),
              dimnames = list(species = species))
  samp = factor(sample(species, size = sample_size,
                       replace = T, prob = rabunds))
  counts = table(samp)
  out[names(counts)] = counts
  return(out)
}

sim_neut_1trial = function(reps, sample_size, rabunds){
  rabunds = namespecies(rabunds)
  num_species = length(rabunds)
  species = names(rabunds)
  if (length(sample_size) == 1){
    sample_size = array(sample_size, dim = c(reps))
  }
  out_mat = array(NA, dim = c(reps, num_species),
                  dimnames = list(rep = 1:reps,
                                  species = species))
  for (rep in 1:reps){
    out_mat[rep,] = sim_neut_1rep(sample_size[rep], rabunds)
  }
  return(out_mat)
}

sim_neut = function(trials, reps, sample_size, rabunds){
  rabunds = namespecies(rabunds)
  species = names(rabunds)
  num_species = length(species)
  out_mat = array(NA, dim = c(reps, num_species, trials),
                  dimnames = list(rep = 1:reps,
                                  species = species,
                                  trial = 1:trials))
  for (trial in 1:trials){
    out_mat[,,trial] = sim_neut_1trial(reps, sample_size, rabunds)
  }
  return(out_mat)
}

inter_mat_1trial = function(rxs_mat, method = "spearman"){
  # Return the covariance, correlation, or other interaction metric
  # as a matrix.
  if (method %in% c("pearson", "spearman", "kendall")){
    return(cor(rxs_mat, method = method))
  }
  if (method %in% c("covar", "effint")){
    covar = cov(rxs_mat)
    ifelse(method == "effint", return(-solve(covar)), return(covar))
  }
  if (method == "mutinfo"){
    return(mi.Dirichlet(rxs_mat, a = 0))
  }
  # Implement additional measures here by adding an
  # if (method %in% ...) statement.
  stop(simpleError(paste("Unrecognized method passed to measure:",
                         method)))
}

inter_mat = function(rxs_mats, method = "spearman"){
  num_reps = dim(rxs_mats)[1]
  num_species = dim(rxs_mats)[2]
  num_trials = dim(rxs_mats)[3]
  species = dimnames(rxs_mats)$species
  
  if (is.null(species)){
    species = paste("SP", 1:num_species, sep = '')
  }
  
  out_mat = array(NA, dim = c(num_species, num_species, num_trials),
                  dimnames = list(species = species,
                                  species = species,
                                  trial = 1:num_trials))
  for (trial in 1:num_trials){
    out_mat[,,trial] = inter_mat_1trial(rxs_mats[,,trial],
                                        method = method)
  }
  return(out_mat)
}

percentile = function(obs, pop){
  species = dimnames(pop)$species
  num_species = length(species)
  out_mat = array(NA, dim = c(num_species, num_species),
                  dimname = list(species = species,
                                 species = species))
  for (i in species){
    for (j in species){
      sample = pop[i, j,]
      if (all(is.na(sample))){
        out_mat[i, j] = NA
      }
      else{
        out_mat[i, j] = 
          ecdf(pop[i, j,])(obs[i, j])
      }
    }
  }
  diag(out_mat) = NA
  return(out_mat)
}

signless_dist = function(x, method = 'euclidian'){
  unidirectional = abs(x - 0.5) + 0.5
  return(dist(unidirectional, method = method))
}

binary_dist = function(x, cutoff = 0.01){
  mat = array(0, dim = dim(x))
  mat[x < cutoff | x > 1 - cutoff] = 1
  return(dist(mat, method = 'binary'))
}

percentile_heatmap = function(mat, cutoff1 = 0.01, cutoff2 = 0.001, cluster = T, ...){
  note = array(NA, dim = dim(mat))
  note[signif(mat, cutoff1)] = '.'
  note[signif(mat, cutoff2)] = '*'
  breaks = c(seq(0, 0.01, 0.0002),
             seq(0.02, 0.98, 0.01),
             seq(0.99, 1, 0.0002))
  color_scheme = colorpanel(length(breaks) - 1,
                            'red', 'black', 'green')
  if (cluster == T){
    Rowv = T
    dendrogram = 'row'
  }
  else {
    Rowv = F
    dendrogram = 'none'
  }
  heatmap.2(mat,
            Rowv = Rowv, dendrogram = dendrogram, symm = T,
            distfun = dist,
            breaks = breaks, col = color_scheme,
            trace = 'none', density.info = 'density',
            cellnote = note, notecol = 'black',
            keysize = 1.5,
            ...)
}

lower_left = function(mat){
  lower_left = array(NA, dim = dim(mat), dimnames = dimnames(mat))
  lower_left[t(upper.tri(mat))] = mat[t(upper.tri(mat))]
  return(lower_left)
}

perc_signif = function(mat, cutoff = 0.01){
  total_comps = length(which(! is.na(mat)))
  signif_comps =
    length(which(significant(mat, cutoff) == T))
  return(signif_comps / total_comps)
}

signif = function(mat, cutoff = 0.01){
  mat < cutoff | mat > 1 - cutoff
}

subsample = function(rxs_mat, size){
  reps = dim(rxs_mat)[1]
  out_mat = array(NA, dim = dim(rxs_mat),
                  dimnames = dimnames(rxs_mat))
  for (rep in 1:reps){
    counts = rxs_mat[rep,]
    rabunds = counts / sum(counts)
    out_mat[rep,] = sim_neut_1rep(sample_size = size,
                                  rabunds = rabunds)
  }
  return(out_mat)
}
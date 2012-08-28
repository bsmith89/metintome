library(gplots)

namespecies = function(x, replace = F){
  if (! is.null(names(x)) && ! replace){
    # Do nothing
    return(x)
  }
  else {
    names(x) = paste("SP", 1:length(x), sep = '')
    return(x)
  }
}

# rxs stands for replicate-by-species (matrix)
rel_abund = function(rxs_mat){
  col_totals = apply(rxs_mat, c(2), sum)
  grand_total = sum(col_totals)
  rel_abunds = col_totals / grand_total
  rel_abunds = namespecies(rel_abunds)
  return(rel_abunds)
}

sim_neut_1rep = function(sample_size, rel_abunds){
  rel_abunds = namespecies(rel_abunds)
  species = names(rel_abunds)
  samp = factor(sample(names(rel_abunds), size = sample_size,
                replace = T, prob = rel_abunds))
  return(tabulate(samp, nbins = length(rel_abunds)))
}

sim_neut_1trial = function(reps, sample_size, rel_abunds){
  rel_abunds = namespecies(rel_abunds)
  num_species = length(rel_abunds)
  species = names(rel_abunds)
  if (length(sample_size) == 1){
    sample_size = array(sample_size, dim = c(reps))
  }
  out_mat = array(NA, dim = c(reps, num_species),
                  dimnames = list(rep = 1:reps,
                                  species = species))
  for (rep in 1:reps){
    out_mat[rep,] = sim_neut_1rep(sample_size[rep], rel_abunds)
  }
  return(out_mat)
}

sim_neut = function(trials, reps, sample_size, rel_abunds){
  rel_abunds = namespecies(rel_abunds)
  species = names(rel_abunds)
  num_species = length(species)
  out_mat = array(NA, dim = c(reps, num_species, trials),
                  dimnames = list(rep = 1:reps,
                                  species = species,
                                  trial = 1:trials))
  for (trial in 1:trials){
    out_mat[,,trial] = sim_neut_1trial(reps, sample_size, rel_abunds)
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
    ifelse(method == "effint", return(-1/covar), return(covar))
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

percentile_heatmap = function(mat, cutoff1 = 0.01, cutoff2 = 0.001){
  note = array(NA, dim = dim(mat))
  note[signif(mat, cutoff1)] = '.'
  note[signif(mat, cutoff2)] = '*'
  breaks = c(seq(0, 0.01, 0.0002),
             seq(0.01, 0.99, 0.01),
             seq(0.99, 1, 0.0002))
  color_scheme = colorpanel(length(breaks) - 1,
                            'red', 'black', 'green')
  heatmap.2(mat, Rowv = T, dendrogram = 'none', symm = T,
            breaks = breaks,
            col = color_scheme,
            trace = 'none', density.info = 'none',
            cellnote = note, notecol = 'black',
            keysize = 1.5)
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
    rel_abunds = counts / sum(counts)
    out_mat[rep,] = sim_neut_1rep(sample_size = size,
                                  rel_abunds = rel_abunds)
  }
  return(out_mat)
}
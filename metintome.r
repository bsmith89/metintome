library(gplots)

sim_neutral = function(reps, sample_size, rel_abunds)
{
  names = names(rel_abunds)
  species = length(names)
  if (length(sample_size) == 1)
  {
    sample_size = array(sample_size, reps)
  }
  data_out = array(0, dim = c(reps, length(rel_abunds)),
                   dimnames = list(1:reps, names))
  for (i in 1:reps)
  {
    counts = tabulate(sample(1:species,
                             prob = rel_abunds,
                             size = sample_size,
                             replace = T),
                      nbins = species)
    data_out[i,] = counts
  }
  data_out = data.frame(data_out)
  return(data_out)
}


inter = function(y)
{
  # Return a matrix of effective interaction terms between
  # species.
  return(-(1 / cov(y)))
  # return(cov(y))
}

mean_rel_abunds = function(y)
{
  # Return the mean relative abundance where replicates are given
  # equal weight regardless of sample-size.
  totals = matrix(apply(y, c(1), sum),
                  nrow = nrow(y), ncol = ncol(y), byrow = F)
  rel_abunds = y / totals
  return(apply(rel_abunds, c(2), mean))
}

rel_abunds = function(y)
{
  col_totals = apply(y, c(2), sum)
  grand_total = sum(col_totals)
  return(col_totals / grand_total)
}

sample_neutral = function(n, reps, sample_size, rel_abunds)
{
  # Return a list of 3D matrices of stacked
  # interaction/covariance/correlation, matrices from a zero
  # interaction neutral model of community assembly
  # (see sim_neutral).  *n* is the number of interaction matrices
  # to sample, *reps* is the number of independent replicates to
  # pull in calculating each neutral interaction matrix,
  # *sample_size* is the number of individuals to count in each
  # replicate, and *rel_abunds* is the expected relative abundance
  # of each species.
  if (length(sample_size) == 1)
  {
    sample_size = array(sample_size, dim=c(reps))
  }
  species = length(rel_abunds)
  names = names(rel_abunds)
  dimnames = list(names, names)
  interaction = array(0, dim = c(species, species, n), dimnames = dimnames)
  spearman = array(0, dim = c(species, species, n), dimnames = dimnames)
  kendall = array(0, dim = c(species, species, n), dimnames = dimnames)
  covar = array(0, dim = c(species, species, n), dimnames = dimnames)
  for (i in 1:n)
  {
    simulation = sim_neutral(reps, sample_size, rel_abunds)
    interaction[,,i] = inter(simulation)
    spearman[,,i] = cor(simulation, method = "spearman")
    kendall[,,i] = cor(simulation, method = "kendall")
    covar[,,i] = cov(simulation)
  }
  return(list(inter = interaction, spearman = spearman,
              kendall = kendall, cov = covar))
}

analyze = function(y, sim_data = NULL,
                        n = 300, assign_intra = T)
{
  # Takes an abundance matrix *y* and returns percentiles on all
  # interaction metrics.
  reps = nrow(y)
  sample_size = apply(y, c(1), sum)
  obs_rel_abunds = rel_abunds(y)
  species = length(obs_rel_abunds)
  obs_spearman = cor(y, method = "spearman")
  obs_kendall = cor(y, method = "kendall")
  obs_inter = inter(y)
  obs_cov = cov(y)
  if (is.null(sim_data))
  {
    sim_data = sample_neutral(n, reps, sample_size,
                              obs_rel_abunds)
  }
  dimnames = list(names(y), names(y))
  percentile_spearman = array(0, dim=c(species, species),
                              dimnames = dimnames)
  percentile_kendall = array(0, dim=c(species, species),
                             dimnames = dimnames)
  percentile_inter = array(0, dim=c(species, species),
                              dimnames = dimnames)
  percentile_cov = array(0, dim=c(species, species),
                              dimnames = dimnames)
  for (i in 1:species)
  {
    for (j in 1:species)
    {
      if (i == j)
      {
          if (assign_intra == T)
          {
            percentile_spearman[i,j] = 0.5
            percentile_kendall[i,j] = 0.5
            percentile_inter[i,j] = 0.5
            percentile_cov[i,j] = 0.5
          }
          else
          {
            percentile_spearman[i,j] = NA
            percentile_kendall[i,j] = NA
            percentile_inter[i,j] = NA
            percentile_cov[i,j] = NA
          }
          next
      }
      percentile_spearman[i,j] = 
        ecdf(sim_data$spearman[i,j,])(obs_spearman[i,j])
      percentile_kendall[i,j] = 
        ecdf(sim_data$kendall[i,j,])(obs_kendall[i,j])
      percentile_inter[i,j] = 
        ecdf(sim_data$inter[i,j,])(obs_inter[i,j])
      percentile_cov[i,j] = 
        ecdf(sim_data$cov[i,j,])(obs_cov[i,j])
    }
  }
  return(list(sim_data = sim_data,
              observed = list(spearman = obs_spearman,
                              kendall = obs_kendall,
                              cov = obs_cov,
                              inter = obs_inter),
              percentile = list(spearman = percentile_spearman,
                                 kendall = percentile_kendall,
                                 inter = percentile_inter,
                                 cov = percentile_cov)))
}

percentile_heatmap = function(y, cutoff = 0.025)
{
  sign = array(NA, dim = c(nrow(y), ncol(y)))
  sign[which(y < cutoff | y > 1 - cutoff)] = '*'
  heatmap.2(y, Rowv = T, dendrogram = 'none', symm = T,
            breaks = c(0, 0.01, 0.025,
                       seq(0.05, 0.95, 0.01),
                       0.975, 0.99, 1),
            col = redgreen, trace = 'none', density.info = 'none',
            cellnote = sign, notecol = 'black')
}

# n = 300
# reps = 100
# sample_size = 1000
# species = 20
# true_rel_abunds = array(1 / species, dim = c(species))
# obs = sim_neutral(reps, sample_size = sample_size,
#                   rel_abunds = true_rel_abunds)
# interaction_results = interactions(obs, n = n, assign_intra = F)
# percentile_heatmap(interaction_results$cov)
# percentile_heatmap(interaction_results$kendall)
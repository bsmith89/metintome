sample2 = function(n, rel_abunds)
{
  # Returns a list of counts in a random sample of a
  # population with relative abundances for species i given
  # in rel_abunds[i]
  sample = NULL
  rand = runif(n)
  for (i in 1:n)
  {
    r = rand[i]
    pointer = 0
    for (j in 1:length(rel_abunds))
    {
      pointer = pointer + rel_abunds[j]
      if (r < pointer)
      {
        sample = c(sample, j)
        break
      }
    }
  }
  return(sample)
}

sim_neutral = function(reps, sample_size, rel_abunds)
{
  if (length(sample_size) == 1)
  {
    sample_size = array(sample_size, reps)
  }
  data_out = matrix(0, nrow = reps,
                    ncol = length(rel_abunds))
  for (i in 1:reps)
  {
    data_out[i,] = tabulate(sample2(sample_size[i],
                                    rel_abunds))
  }
  data_out = data.frame(data_out)
  return(data_out)
}

inter = function(y)
{
  # Return a matrix of effective interaction terms between
  # species.
  return(-(1 / cor(y)))
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
  interaction = array(0, dim = c(species, species, n))
  spearman = array(0, dim = c(species, species, n))
  kendall = array(0, dim = c(species, species, n))
  for (i in 1:n)
  {
    simulation = sim_neutral(reps, sample_size, rel_abunds)
    interaction[,,i] = inter(simulation)
    spearman[,,i] = cor(simulation, method = "spearman")
    kendall[,,i] = cor(simulation, method = "kendall")
  }
  return(list(interaction = interaction, spearman = spearman,
              kendall = kendall))
}

interactions = function(y, simulated_data = NULL,
                        n = 300, assign_intra = T)
{
  # Takes an abundance matrix *y* and returns percentiles on all
  # interaction metrics.
  if (is.null(simulated_data))
  {
    reps = nrow(y)
    sample_size = apply(y, c(1), sum)
    obs_rel_abunds = rel_abunds(y)
    species = length(obs_rel_abunds)
    obs_spearman = cor(y, method = "spearman")
    obs_kendall = cor(y, method = "kendall")
    obs_interaction = inter(y)
    simulated_data = sample_neutral(n, reps, sample_size,
                                    obs_rel_abunds)
  }
  percentile_spearman = array(0, dim=c(species, species))
  percentile_kendall = array(0, dim=c(species, species))
  percentile_interaction = array(0, dim=c(species, species))
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
            percentile_interaction[i,j] = 0.5
          }
          else
          {
            percentile_spearman[i,j] = NA
            percentile_kendall[i,j] = NA
            percentile_interaction[i,j] = NA
          }
          next
      }
      percentile_spearman[i,j] = 
        ecdf(simulated_data$spearman[i,j,])(obs_spearman[i,j])
      percentile_kendall[i,j] = 
        ecdf(simulated_data$kendall[i,j,])(obs_kendall[i,j])
      percentile_interaction[i,j] = 
        ecdf(simulated_data$interaction[i,j,])(obs_interaction[i,j])
    }
  }
  return(list(spearman = percentile_spearman,
              kendall = percentile_kendall,
              interaction = percentile_interaction,
              simulated_data = simulated_data))
}

obs = sim_neutral(100, sample_size = 1000,
                  rel_abunds = c(0.1, 0.1, 0.1, 0.1, 
                                 0.25, 0.25, 0.05, 0.05))

interaction_results = interactions(obs, n = 300)

# sample50_1000 = sample_neutral(300, 50, 1000,
#                           c(0.1, 0.1, 0.1, 0.1,
#                             0.25, 0.25, 0.05, 0.05))
# sample100_1000 = sample_neutral(300, 100, 1000,
#                            c(0.1, 0.1, 0.1, 0.1,
#                              0.25, 0.25, 0.05, 0.05))
# sample50_2000 = sample_neutral(300, 50, 2000,
#                                c(0.1, 0.1, 0.1, 0.1, 
#                                  0.25, 0.25, 0.05, 0.05))
# 
# par(mfrow=c(3,1))
# hist(sample50_1000$spearman[1,2,], xlim = c(-0.6, 0.3))
# hist(sample50_2000$spearman[1,2,], xlim = c(-0.6, 0.3))
# hist(sample100_1000$spearman[1,2,], xlim = c(-0.6, 0.3))

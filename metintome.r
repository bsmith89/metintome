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

sim_neutral = function(n, sample_size, rel_abunds)
{
  if (length(sample_size) == 1)
  {
    sample_size = array(sample_size, n)
  }
  data_out = matrix(0, nrow = n,
                    ncol = length(rel_abunds))
  for (i in 1:n)
  {
    data_out[i,] = tabulate(sample2(sample_size[i],
                                    rel_abunds))
  }
  data_out = data.frame(data_out)
  return(data_out)
}

inter_mat = function(y)
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












ab_corr_obs50_1000 = NULL
for (rep in 1:300)
{
  ab_corr_obs50_1000 = c(ab_corr_obs50_1000,
                  cor(sim_neutral(50, 1000, 
                                  c(0.05, 0.3, 0.2,
                                    0.05, 0.05, 0.35)),
                      method = "spearman")[1,2])
}
hist(ab_corr_obs50_1000)
qqnorm(ab_corr_obs50_1000)
qqline(ab_corr_obs50_1000)

ab_corr_obs10_1000 = NULL
for (rep in 1:300)
{
  ab_corr_obs10_1000 = c(ab_corr_obs10_1000,
                  cor(sim_neutral(10, 1000, 
                                  c(0.05, 0.3, 0.2,
                                    0.05, 0.05, 0.35)),
                      method = "spearman")[1,2])
}
hist(ab_corr_obs10_1000)
qqnorm(ab_corr_obs10_1000)
qqline(ab_corr_obs10_1000)

ab_corr_obs10_500 = NULL
for (rep in 1:300)
{
  ab_corr_obs10_500 = c(ab_corr_obs10_500,
                  cor(sim_neutral(10, 500, 
                                  c(0.05, 0.3, 0.2,
                                    0.05, 0.05, 0.35)),
                      method = "spearman")[1,2])
}
hist(ab_corr_obs10_500)
qqnorm(ab_corr_obs10_500)
qqline(ab_corr_obs10_500)
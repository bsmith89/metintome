# Interaction analysis of microbiome from non-smokers, both lung and
# oral samples. Values in data set are relative abundance of the given
# OTU from a sample size of 1000.

source("metintome.r")
library(VGAM) # for rpareto()

# # In this analysis I'm combining all 86 samples, both lung and oral,
# # in order to test the ability of the analysis to identify OTUs which
# # co-occur in the same environment.
# 
# data = read.csv("data/both_subsample.csv", header = T, sep = ",")
# data = data[, -c(1:2)]
# data = data * 1000
# totals = apply(data, 2, sum)
# sorted_names = names(sort(totals, decreasing = T))
# data = data[, sorted_names]
# totals = apply(data, 2, sum)
# data_all = data[, which(totals != 0)]
# 
# obs_rxs = data_all
# method = "spearman"
# reps = dim(obs_rxs)[1]
# sample_sizes = apply(obs_rxs, 1, sum)
# obs_rel_abunds = rel_abund(obs_rxs)
# neut_rxs_sample = sim_neut(1000, reps, sample_sizes, obs_rel_abunds)
# obs_inter_mat = inter_mat_1trial(obs_rxs, method = method)
# neut_inter_mat_sample = inter_mat(neut_rxs_sample, method = method)
# comparison = percentile(obs_inter_mat, neut_inter_mat_sample)
# percentile_heatmap(comparison[1:100, 1:100],
#                    cutoff1 = 0.005, cutoff2 = 0.001)




# Here I will create two simultad communities with 480 OTUs a piece. 
# The abundances of each OTU will be pulled from a pareto
# distribution, and then normalized to the total, so that the relative
# abundance sums to 1 I will then simulate 20 replicates from each of
# these meta-communities, and I will run the resulting 40 replicate,
# 500 OTU trial through the previously run analysis. My goal is to see
# if I can identify the OTUs which are both partial to one of the
# meta-communities through clustering over interactions. Basically,
# this is a positive environmental interaction control.

location = 1
shape = 1
num_species = 500

abund1 = rpareto(num_species, location, shape)
abund2 = rpareto(num_species, location, shape)
abund3 = rpareto(num_species, location, shape)
abund4 = rpareto(num_species, location, shape)

abund1 = namespecies(abund1)
abund2 = namespecies(abund2)
abund3 = namespecies(abund3)
abund4 = namespecies(abund4)

rabund1 = abund1/sum(abund1)
rabund2 = abund2/sum(abund2)
rabund3 = abund3/sum(abund3)
rabund4 = abund4/sum(abund4)

obs_reps1 = sim_neut_1trial(10, 1000, rabund1)
obs_reps2 = sim_neut_1trial(10, 1000, rabund2)
obs_reps3 = sim_neut_1trial(10, 1000, rabund3)
obs_reps4 = sim_neut_1trial(10, 1000, rabund4)
data = rbind(obs_reps1, obs_reps2, obs_reps3, obs_reps4)

totals = apply(data, 2, sum)
data_all = data[, which(totals != 0)]
sort(totals)[c(1, seq(50, 500, 50))]





obs_rxs = data_all
method = "spearman"
reps = dim(obs_rxs)[1]
sample_sizes = apply(obs_rxs, 1, sum)
obs_rel_abunds = rel_abund(obs_rxs)
neut_rxs_sample = sim_neut(1000, reps, sample_sizes, obs_rel_abunds)
obs_inter_mat = inter_mat_1trial(obs_rxs, method = method)
neut_inter_mat_sample = inter_mat(neut_rxs_sample, method = method)
comparison = percentile(obs_inter_mat, neut_inter_mat_sample)
percentile_heatmap(comparison[1:100, 1:100],
                   cutoff1 = 0.005, cutoff2 = 0.001)

# I'm now getting much closer to the expected clustering!!!
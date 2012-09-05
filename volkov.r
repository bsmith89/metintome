# An attempt to use Volkov's et al. (2009)'s interaction metric,
# which is based on maximum entropy, to relative abundance data from
# the mouth/lung microbiome.


source("metintome.r")

data = read.csv("data/bal.csv", header = T, sep = ",")
data = data[, -c(1:2)]
data = data * 1000
totals = apply(data, 2, sum)
sorted_names = names(sort(totals, decreasing = T))
data = data[, sorted_names]
totals = apply(data, 2, sum)
data_all = data[, which(totals != 0)]



species_range = 1:25
data_abundant = data[,species_range]
obs_rxs = data_all
eim_obs = inter_mat_1trial(data_abundant, method = 'effint')
reps = dim(obs_rxs)[1]
sample_sizes = apply(obs_rxs, 1, sum)
obs_rel_abunds = rel_abund(obs_rxs)
# neut_rxs_sample = sim_neut(1000, reps, sample_sizes, obs_rel_abunds)
neut_rxs_sample_abundant = neut_rxs_sample[,species_range,]


eim_sample = inter_mat(neut_rxs_sample_abundant,
                       method = 'effint')
comparison = percentile(eim_obs, eim_sample)
percentile_heatmap(comparison, 0.005, 0.001, cluster = F)
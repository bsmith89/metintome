# Interaction analysis of ling microbiome from non-smokers.
# Values in data set are relative abundance of the given OTU from
# a sample size of 1000.
source("metintome.r")

data = read.csv("bal.csv", header = T, sep = ",")
data = data[, -c(1:2)]
data = data * 1000
totals = apply(data, 2, sum)
sorted_names = names(sort(totals, decreasing = T))
data = data[, sorted_names]

# # If I cut out all but the top 20 OTUs
# data_20 = data[,1:20]
# obs_rxs20 = data_20
# method20 = "spearman"
# reps20 = dim(obs_rxs)[1]
# sample_sizes20 = apply(obs_rxs20, 1, sum)
# obs_rel_abunds20 = rel_abund(obs_rxs20)
# neut_rxs_sample20 = sim_neut(1000, reps20, sample_sizes20, obs_rel_abunds20)
# obs_inter_mat20 = inter_mat_1trial(obs_rxs20, method = method20)
# neut_inter_mat_sample20 = inter_mat(neut_rxs_sample20, method = method20)
# comparison20 = percentile(obs_inter_mat20, neut_inter_mat_sample20)
# percentile_heatmap(comparison20[1:40, 1:40],
#                    cutoff1 = 0.005, cutoff2 = 0.001)



data_all = data[, which(totals != 0)]

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

# Pretty cool heatmap, eh?  If there's clustering turned on, you'll 
# see that there exist ~2.5 clusters of OTUs which interact strongly 
# with each other.  These groups show up as reciprocal positive or
# negative interactions.  For instance, OTU A may interact 
# significantly with B (which means vice-versa as well).  A cycle is 
# formed if both A and B interact with a third OTU, C.  The existance 
# of an OTU C which fulfills this description is not trivial.  Even 
# less trivial would be a four OTU cycle, in which A and B are 
# connected to each other at the 1º, 2º (by both C and D), and 3º (by
# C through D) levels.  These graph structures are not expected from
# the existance of false positives in our analysis. However, we see a
# ton of them in our actual data.  In fact I see a cluster of 5 OTUs
# which are linked by every possible connection of 1º through 4º 
# degrees.  Groups of OTUs which are linked by a *large* fraction of
# all possible chains are also obvious.  These quantitative
# definitions of OTU clusters should be simple enough to evaluate from
# a network matrix.  You could also potentially weight the edges of 
# the (undirected) graph by the percentile (relative to simulated, 
# neutral data) or actual strength of the interaction.

# Some follow up questions that are worth asking: (1) Are OTUs which
# are more phylogenetically related (based on 16S) more likely to be
# clustered together as described above?  (2) Are phylogetically
# similar OTUs more likely to interact than not?  Positively or
# negatively.  The answers to many of these questions could be
# visualized by clustering OTUs based on their genetic distance in the
# heatmap.2 function, rather than by the default criteria. I also want
# to consider the prodominance of positive interactions relative to
# the neutral assumptions rather than negative interactions. It's
# interesting that positive interactions could mean a few different
# things.  There are direct relationships between species A and B
# which could result in positive interaction scores: (1) What we're
# most interested in, are mutalistic relationships, synergies, in
# which the presence of species A is beneficial to species B (and
# maybe vice-versa).  (2) Another possibility is that species A and B
# share environmental optima, so where we find A we expect to find B. 
# Indirectly: (3) a third species, C, could have a mutualistic
# relationship with both A and B, thereby giving A and B an
# interaction which is not directly based on their physiology.  (4+)
# Just about any combination of environmental and mutualistic
# relationships chaining A to B by one or more additional species
# could give A and B a significant interaction score. This means that
# given the existance of true interactions, blocks of >2 OTUs all of
# which interact are not surprising, given that higher order, indirect
# interactions, link all species to all other species.  This is not to
# be expected from coincidental (false positive) correlations.

# Another questions, however, is why we don't see more negative
# interactions.  If species A and B both occupy the same niche,
# shouldn't we expect their interaction score to be negative, or is
# this canelled out by the positive interaction inherent in their
# similar environmental optima? Interestingly, we do seem to have some
# tenuous clusters of negatively interacting OTUs.  It may be
# significant that we see these negative interactions between members
# of two clusters of highly (positively) connected OTUs. Perhaps we're
# seeing environmentally driven negative correlation, or maybe we're
# seeing negative correlations between A and C, and B and C resulting
# in positive correlations between A and B. I'd like to see a 
# reclustering of the heat map where positive and negative 
# interactions are clustered together.

comparison_signless = abs(comparison - 0.5) + 0.5
percentile_heatmap(comparison_signless[1:100, 1:100],
                   cutoff1 = 0.005, cutoff2 = 0.001)
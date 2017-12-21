#setwd("/afs/cad.njit.edu/u/g/r/gr93/public_html/uploads")
library(GeoDE)
data(expression_data)
data(sampleclass)
data(gammas)
#data(GeneOntology_BP.gmt)
chdir_analysis_example <- chdirAnalysis(example_expression_data,example_sampleclass,example_gammas,CalculateSig=TRUE,nnull=10)
#PAEAtest <- PAEAAnalysis(chdir_analysis_example$chdirprops, gmt[1:100], example_gammas)


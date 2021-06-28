### Calculate reads remaining at each step of DADA2
### Full .RData stored in login node, deleted off local machine after this code was run (too large)
### 6.23.2021

reads_remaining <- data.frame(NULL)
dataFiles <- list.files("/Users/kbja10/Documents/Cornell/Research/Round goby/Populations genetics: eDNA/Field_experiment_2019/data_analysis/out_dada2_6.23", pattern="Nmel*", full.names="TRUE")
dataFiles <- paste(dataFiles,"/Dada2.RData",sep="")

for (file in dataFiles) {
  e = local({load(file); environment()})
  tools:::makeLazyLoadDB(e, "New")
  reads_remaining <- rbind(reads_remaining,as.data.frame(e$track))
}
# write.csv(reads_remaining, "/Users/kbja10/Github/eDNA_goby_popgen/datasets/dada2_reads_remaining.csv")

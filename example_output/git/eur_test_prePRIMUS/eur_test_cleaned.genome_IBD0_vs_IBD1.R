data <-read.table("../example_output/git/eur_test_prePRIMUS/eur_test_cleaned.genome",header=T)
jpeg("../example_output/git/eur_test_prePRIMUS/eur_test_cleaned.genome_IBD0_vs_IBD1.jpeg", height=480, width=480)
plot(data[,'Z0'],data[,'Z1'],xlab="IBD0",xlim=c(0,1),ylim=c(0,1),ylab="IBD1",main="IBD0 vs IBD1 for eur_test")
dev.off()

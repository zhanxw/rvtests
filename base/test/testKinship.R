library(kinship2)

pedAll <- pedigree(id=sample.ped$id,
                   dadid=sample.ped$father, momid=sample.ped$mother,
                   sex=sample.ped$sex, famid=sample.ped$ped)

print(pedAll)
ped2basic <- pedAll['2']
kin2 <- kinship(ped2basic)

#write.table(sample.ped[,seq(5)], file = "testKinship2.ped", row.names = FALSE, col.names= FALSE)
write.table(subset(sample.ped[,seq(5)], ped == 2), file = "testKinship2.ped", row.names = FALSE, col.names= FALSE)
write.table(file = "R.kinship2.correct", kin2, col.names = FALSE, quote = FALSE)

pdf(file = "R.kinship2.pdf")
plot(ped2basic)
dev.off()


## a <- read.table("testKinship.ped", header = FALSE, stringsAsFactors = FALSE)
## a <- read.table("testKinship.ped", header = FALSE)

## colnames(a) <- c("famid", "id", "dadid", "momid")
## head(a)
## pedAll <- pedigree(id = as.numeric(a$id),
##                    famid = as.numeric(a$famid),
##                    dadid = as.numeric(a$dadid),
##                    momid = as.numeric(a$momid),
##                    sex = rep(1, nrow(a)))

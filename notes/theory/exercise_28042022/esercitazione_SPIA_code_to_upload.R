rm(list=ls());verbose<-TRUE; options(max.print=1000);getOption("max.print")

# install.packages("SPIAssay")
## set working dir
setwd("C:/Users/Maurizio/Desktop/HumanGenomics/notes/exercise_28042022"); getwd();

library(SPIAssay)

## set param for high MAF SNPs
spia_parameters <- list(Pmm = 0.1, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.7)

## load genotype data
genotype_data <- read.table("SPIA_input_genotype.tsv", header = TRUE, sep = "\t", as.is = TRUE ); 
dim(genotype_data); head(genotype_data);

## load high MAF SNPs unique IDs
selected_snps <- read.delim("SPIA_selected_SNPs.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,1]
dim(genotype_data); head(genotype_data);
length(selected_snps); head(selected_snps);

##TASK 1
## select genotype data based on SNP IDs
genotype_data_selected <- genotype_data[selected_snps,]
dim(genotype_data_selected);
genotype_data_selected

## run test
SPIA_selected <-  SPIATest(x = genotype_data_selected, row.names = FALSE, test.param = spia_parameters) 
## visualize distances
SPIAPlot(SPIA_selected)
## check output matrix
head(SPIA_selected["SPIAresult"]);
##look at distances histogram
hist(as.numeric(SPIA_selected$SPIAresult[,3]),breaks = 100,xlim=c(0,1))
## save output matrix
write.table(SPIA_selected["SPIAresult"],paste("output_","selected",".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)

##subsampling of high MAF SNPs: 100
c100<-sample(c(1:nrow(genotype_data_selected)),100,replace = FALSE);


SPIA_selected <-  SPIATest(x = genotype_data_selected[c100,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected)
## ...

##subsampling of high MAF SNPs, 100, 80, ...
c80 <- sample(c(1:nrow(genotype_data_selected)),80,replace = FALSE);
SPIA_selected <-  SPIATest(x = genotype_data_selected[c80,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected)

c60 <- sample(c(1:nrow(genotype_data_selected)),80,replace = FALSE);
SPIA_selected <-  SPIATest(x = genotype_data_selected[c60,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected)

c40 <- sample(c(1:nrow(genotype_data_selected)),80,replace = FALSE);
SPIA_selected <-  SPIATest(x = genotype_data_selected[c40,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected)

c20 <- sample(c(1:nrow(genotype_data_selected)),80,replace = FALSE);
SPIA_selected <-  SPIATest(x = genotype_data_selected[c20,], row.names = FALSE, test.param = spia_parameters) 
SPIAPlot(SPIA_selected)


##TASK 2
SPIA_all <-  SPIATest(x = genotype_data, row.names = F, test.param = spia_parameters) #all the SNPs
SPIAPlot(SPIA_all)
c100<-sample(c(1:nrow(genotype_data)),100,replace = FALSE);
genotype_data_random <- genotype_data[c100,]
SPIA_random <-  SPIATest(x = genotype_data_random, row.names = F, test.param = spia_parameters) 
SPIAPlot(SPIA_random)
hist(as.numeric(SPIA_random$SPIAresult[,3]),breaks = 100,xlim=c(0,1))


##TASK 3
## load somatic data
somatic_data <- read.table("SPIA_somatic_genotype.tsv", header = TRUE, sep = "\t", as.is = TRUE ); 
dim(somatic_data); head(somatic_data);
SPIA_somatic <-  SPIATest(x = somatic_data, row.names = F, test.param = spia_parameters) 
SPIAPlot(SPIA_somatic)

# look at SPIAresult
write.table(SPIA_somatic["SPIAresult"], file = paste("output_","somatic",".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE)
table_somatic <- read.table("output_somatic.txt", header = TRUE)
View(table_somatic)

#RESULTS OBTAINED 
#TODO add description

##TASK 4
## change parameters
## set param for high MAF SNPs, exploring different values for Pmm, Pmm_nonM, nsigma, nsigma_nonM (one at the time)
spia_parameters <- list(Pmm = 0.2, nsigma = 2, Pmm_nonM = 0.6, nsigma_nonM = 3, PercValidCall=0.7)
# ... #TODO

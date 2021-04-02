library(GEOquery)
library(limma)
library(dplyr)
library(ggplot2)

# Deleting all objects in the workspace
rm(list = ls()) 

###########################################################
# Loading in the GSE146446  data.
###########################################################
GSE146446 <- getGEO("GSE146446")

###########################################################
# Pulling out gene expression data and pheno type data. 
###########################################################
GSE146446.expr <- exprs(GSE146446[[1]])
GSE146446.p <- pData(GSE146446[[1]])

###########################################################
# Generating a boxplot of the first 10 samples, checking
# if data is alreay in log2 scale
###########################################################
boxplot(GSE146446.expr[,1:10], ylab = "Expression Data in log2", main = "Processed Data",
        xlab = "Samples")
# Based on boxplot samples the data is on the log2 scale and therefore
# normailzed
###########################################################

# Numebr of samples
ncol(GSE146446.expr)

# Number of probes
nrow(GSE146446.expr)

###########################################################
# Extracting the columns that contains the categories that
# are going to be compared. The data has to be processed first
###########################################################

# Creating Vector to combine the time and the responder status
response <- paste(GSE146446.p$source_name_ch1, GSE146446.p$characteristics_ch1.5, sep ="_")

# Constructing design Matrix
design <- model.matrix(~0+response)

# Changing Column Names
colnames(design) <- c("T0DLX0", "T0DLX1", "T8DLX0","T8DLX1","T0PLB0","T0PLB1", "T8PLB0",
                                            "T8PLB1")
# Number of samples in each group
table(response)[1:2]

# fitting a linear model to each row of the matrix
fit <- lmFit(GSE146446.expr, design)

# specifying the contrast
contrast.matrix <- makeContrasts(T0DLX0 - T0DLX1, levels = design )

# fiting model based on contrasts (T0DLX0 - T0DLX1)
fit2 <- contrasts.fit(fit, contrast.matrix)

# calculatign the moderated t-statistics by moderating standard errors
# toward a common value
fit2 = eBayes(fit2) 

# get all probes with FDR < 1, sorted by p-value
tt.1 <- topTable(fit2,sort.by = "p", p.value = 1, number = 50)

# number of diffrerntialy expressed probes
nrow(tt.1)

# Getting overall FDR
overallFDR <- tt.1[50,5];overallFDR

## get probes that are differentially expressed
m <- match(rownames(tt.1), rownames(GSE146446.expr))
X <- GSE146446.expr[m,]

# probe with the lowest adjusted p-value
# creating data frame with expression values and response
probe <- rownames(tt.1)[1]
m <- match(probe, rownames(GSE146446.expr))
df <- data.frame(expr = GSE146446.expr[m,], response = response)

# group data frame by response type, then summarize by expression value (new approach)
means <- df %>% group_by(response) %>% summarize(mean = mean(expr))

# convert to FC #
logFC <- tt.1$logFC[1]

## visualize ##
FC <- paste0("FC = ", round(2**logFC, 2))
main <- paste0("Expression of ", probe, ", ", FC, " FDR = ", overallFDR )

# Constucting box plot 
ggplot(df, aes(x = response, y = expr, fill = response)) + geom_boxplot() +
  scale_x_discrete(limits=c("DLX_T0_0_response: 0","DLX_T0_1_response: 1")) +
  ylab("log2 expression") + ggtitle(main) +
  theme_classic() + theme(legend.position = "none")

###############################################################
# Finding the genes associated with all probes
###############################################################
platform <- annotation(GSE146446[[1]])   
platform

# downloading the platform data 
pl <- getGEO(platform)
pl <- Table(pl)

# Finding the gene for the desired probe 
probes <- rownames(tt.1)
m <- match(probes, pl$ID)
genes <- pl$Symbol[m]

# data frame containing the probe names, logFC and adjusted p-values
dfInfo <- tt.1 %>% select(`logFC`, `adj.P.Val`)

# Adding gene names to data frame 
geneDF <- mutate(dfInfo, Names = genes); geneDF

# out putting the top 5 probes
geneDF[1:5,]

# removing all genes that are empty strings from the data frame
geneDF <- geneDF %>% na_if("") %>% na.omit; geneDF

# Getting the gene names specifically
geneList <- geneDF %>% select(Names)

# get a unique set of genes 
genes <- unique(geneList); genes

# Writing genes in a file to upload on DAVID
write.table(genes, row.names = FALSE, quote = FALSE, file="genes.txt")

# Brief Summary
# The GSE series analyzed was GSE146446. The samples analyzed are the response
# of antidepressed patients on the treament they recevied, whether it was effective
# or not. The number of samples are 406, the number of probes that were profiled
# are 47323, the number of samples in each group are 21 (DLX reponse 0) and 75
# (DLX reponse 1), there were 50 differentially expressed probes with a FDR of .99992, 
# the names of the top 3 genes are FABP1, C4orf6, RAET1E, and the top go terms
# assocaited with the phenotype cellular resposne to hypoxia, negative regulation
# of fatty acid oxidation, positive regulation of proteolysis, 
# negative regulation of apoptotic process, response to testosterone
# positive regulation of apoptotic process, response to calcium ion, 
# IRE1-mediated unfolded protein response. 

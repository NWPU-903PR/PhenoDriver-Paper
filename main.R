library(PhenoDriver)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Inpute data
expN <- read_csv("./data/BRCA_NT_Rawcounts.csv")
expT <- read_csv('./data/BRCA_TP_Rawcounts.csv')
mutationData <- read.delim('./data/MC3_BRCA.maf')
reactome <- read.delim('./data/NCBI2Reactome.txt',header = F)
STN <- read.delim('./data/SignalTransductionNet.txt')

# Data per-processing
expressionData <- list(Normal = expN, Tumor = expT)
cancerData <- dataPerProcessing(expressionData = expressionData,
                                mutationData = mutationData,
                                annotationcol = 2,
                                removeDuplicates = 2,
                                removeSilenceSNP = F)

# Differential expression analysis
diffexpgene <- DEA(exp = cancerData$expressionData,
                   parallelworker = 40,
                   genenamecol = 2,
                   annotationcol = 2)

reactome <- reactome[reactome$V6 == 'Homo sapiens', 1:4]
reactome <- unique(reactome)
reactome$V1 <- AnnotationDbi::select(org.Hs.eg.db, keys = reactome$V1, columns = 'SYMBOL')$SYMBOL
reactome <- unique(reactome)
reactome <- reactome[!(is.na(reactome$V1)),]
reactome <- reactome[reactome$V4 %in% names(table(reactome$V4)[table(reactome$V4) > 10 &
                                                                 table(reactome$V4) < 200]),]

# Get personalized abnormal pathways
enrichreactomeRes <- pAbnormalPathway(diff = diffexpgene,
                                      reactome = reactome,
                                      parallelworker = 40)

# Calculating driving force matrix for each patient
genescore_individual <- calDrivingForce(STN,
                                        diffexpgene,
                                        enrichreactomeRes,
                                        reactome,
                                        15)

# Get driver genes from driving force matrix
drivergenes <- getDriverGenes(drivingforcelist = genescore_individual,
                              cancerData = cancerData)


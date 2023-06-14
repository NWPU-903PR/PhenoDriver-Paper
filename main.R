library(PhenoDriverR)
library(readr)
library(AnnotationDbi)
library(org.Hs.eg.db)

# Inpute data
expN <- read_csv("./data/BRCA_NT_Rawcounts.csv")
expT <- read_csv('./data/BRCA_TP_Rawcounts.csv')
mutationData <- read_tsv('./data/MC3_BRCA.maf.gz', progress = TRUE, col_types = cols())
reactome <- read.delim('./data/NCBI2Reactome.txt', header = F)
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
                   parallelworker = 10,
                   genenamecol = 2,
                   annotationcol = 2)

reactome <- reactome[reactome$V6 == 'Homo sapiens', 1:4]
reactome <- unique(reactome)
reactome$V1 <- select(org.Hs.eg.db, keys = reactome$V1, columns = 'SYMBOL')$SYMBOL
reactome <- unique(reactome)
reactome <- reactome[!(is.na(reactome$V1)),]
reactome <- reactome[reactome$V4 %in% names(table(reactome$V4)[table(reactome$V4) > 10 &
                                                                 table(reactome$V4) < 200]),]

# Get personalized abnormal pathways
enrichreactomeRes <- pAbnormalPathway(diff = diffexpgene,
                                      reactome = reactome,
                                      parallelworker = 10)

# Calculating driving force matrix for each patient
genescore_individual <- calDrivingForce(network = STN,
                                        diffexpgene = diffexpgene,
                                        enrichreactomeRes = enrichreactomeRes,
                                        reactome = reactome,
                                        parallelworker = 10)

# Get driver genes from driving force matrix
drivergenes <- getDriverGenes(drivingforcelist = genescore_individual$genescore,
                              cancerData = cancerData)

# Get driver-associated abnormal pathways (TP53 as an example)
tp53relatedpathway <- getDriverRelatedPathway(drivingforcelist = genescore_individual$genescore,
                                              enrichreactomeRes = enrichreactomeRes,
                                              drivergenes = drivergenes,
                                              quary = 'TP53')

# Get edgelist of TP53 regulatory paths to some cell cycle associated abnormal pathways
regulation_path_edgelist <- getRegulatoryPath(drivergenes = drivergenes,
                                              zscore_ind = genescore_individual$zscore,
                                              network = STN,
                                              quary = 'TP53',
                                              quary_pathway = c('EML4 and NUDC in mitotic spindle formation',
                                                                'Resolution of Sister Chromatid Cohesion',
                                                                'Mitotic Prometaphase',
                                                                'RHO GTPases Activate Formins',
                                                                'Separation of Sister Chromatids',
                                                                'Amplification of signal from unattached kinetochores via a MAD2 inhibitory signal'))

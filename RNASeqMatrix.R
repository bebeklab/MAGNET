#get the ID of run job
hh <- commandArgs(trailingOnly =TRUE)
args <- unlist(hh)

FindID <- as.numeric(args[1])
randomID <- as.numeric(args[2])

print(args)

print("FindID")
print(FindID)
print("randomID")
print(randomID)


logName <- paste("/var/www/magnet/logs/",FindID, ".txt", sep= "")
if (file.exists(logName)) {
  file.remove(logName)
}

cat(paste(toString(Sys.time()), " - INFO - ", "Starting MAGNET ", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = FALSE, sep = "\n")
#install package to connect to SQL server
library(RMySQL)
library(sendmailR)
#library(psych) 

#Vishal add Check 1
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 1.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")


tryCatch({
#to round decimal
specifydecimal <- function(x, k) format(round(x, k), nsmall=k)

#connect to the database
con <- dbConnect(MySQL(), user="*******", password="********", dbname ="********" , host="localhost")
#get all info in table for that specific run
table <- dbGetQuery(con, paste("SELECT * FROM process_processrequest WHERE concat('',ID * 1) = " , FindID,  " LIMIT 1", sep = ""))
#table <- dbGetQuery(con, "SELECT * FROM process_interaction WHERE request_id = 1113 AND genename_a = 'EP300' AND genename_a = 'EP301'")
#to get the files. don't understand syntax so useless for now.
files <- table[20]
stringOfFiles <- toString(files)
sFiles <- strsplit(stringOfFiles, "\n")
RNASeqOrTCGA <- FALSE
for(i in 1:length(sFiles[[1]])) {
  if(sFiles[[1]][i] == "S'var/www/magnet/TCGA/dummyFile.txt'" | sFiles[[1]][i] == "S'/var/www/magnet/TCGA/dummyFile.txt'") {
    RNASeqOrTCGA = TRUE;
  }
}
#end of useless lines

#Vishal add Check 2
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 2.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")



#get GSE file
cat(paste(toString(Sys.time()), " - INFO - ", "Reading GPL file.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
#GPLFile  <- read.delim("bim_GPL1261_sample.txt" , header=F, stringsAsFactors=FALSE)
GPLFile <- read.delim(paste("/var/www/magnet/uploads/gpl_", FindID, ".txt", sep = ""), header=F, stringsAsFactors=FALSE)
#if GPL doesn;t exist, it is a TCGA File, so we can just import it regulary. no dummy gene list.
if(exists("GPLFile") == FALSE | RNASeqOrTCGA == TRUE) {
  #RNASeqExpressionTable <- read.delim("Documents/College/BioInformatics_Research/MAGNET_ALL/Downloading_TCGA/Final_Expression/acc_tcgaRNASeqExpression.txt", header=F, stringsAsFactors=FALSE)
  cat(paste(toString(Sys.time()), " - INFO - ", "No GPL File. Reading GSE file now.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  #RNASeqExpressionTable <- read.delim("gbm_tcgaRNASeqExpression.txt", header=F, stringsAsFactors=FALSE)
  RNASeqExpressionTable <- read.delim(paste("/var/www/magnet/uploads/gse0_", FindID, ".txt", sep = ""), header=F, stringsAsFactors=FALSE)
} else { #if it is not we have to combine the table. link probe IDs to gene names
  #Find the start of genes in GPL
  for(i in 1:nrow(GPLFile)) {
    if(GPLFile[i,1] == "ID") {
      IDLocation = i;
      break;
    }
  }
  GPLFile <- read.delim(paste("/var/www/magnet/uploads/gpl_", FindID, ".txt", sep = ""), header=F, stringsAsFactors=FALSE, skip = IDLocation-1)
  #reallocate GPLFile
  #GPLFile  <- read.delim("bim_GPL1261_sample.txt" , header=F, stringsAsFactors=FALSE, skip = IDLocation-1)
  cat(paste(toString(Sys.time()), " - INFO - ", "Reading GSE file.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  #GSEFile <- read.delim("bim_GSE_19338_sample.txt", header=F, stringsAsFactors=FALSE)
  GSEFile <- read.delim(paste("/var/www/magnet/uploads/gse0_", FindID, ".txt", sep = ""), header=F, stringsAsFactors=FALSE)
  #Find start of the table in the GSE
  for(i in 1:nrow(GSEFile)) {
    if(GSEFile[i,1] == "!Sample_title") {
      TableStartLocation = i;
      break;
    }
  }
  #reallocate GSEFile
  #GSEFile <- read.delim("bim_GSE_19338_sample.txt", header=F, stringsAsFactors=FALSE, skip = TableStartLocation-1)
  GSEFile <- read.delim(paste("/var/www/magnet/uploads/gse0_", FindID, ".txt", sep = ""), header=F, stringsAsFactors=FALSE, skip = TableStartLocation-1)
  #find the column that contains gene symbol is in the GPL File
  i = 2
  foundGeneHeader = FALSE
  while(ncol(GPLFile) > 2) {
    if(foundGeneHeader == FALSE & (GPLFile[1,i] == "Symbol" || GPLFile[1,i] == "Gene" || GPLFile[1,i] == "Gene Symbol" || GPLFile[1,i] == "GENE_SYMBOL" || GPLFile[1,i] == "Genename")) {
      i = i +1;
      foundGeneHeader = TRUE
    } else {
      GPLFile <- GPLFile[,-i]
    }
  }
  #GPL Files only has probes and gene names now
  #need to get characteristics out of GSE
  for(i in 1:nrow(GSEFile)) {
    if(GSEFile[i,1] == "ID_REF") {
      startLocation = i;
      break;
    }
  }
  count = 1;
  #removing characteristics (to add back in later)
  while(count < startLocation) {
    if(count == 1) {
      AllChars <- GSEFile[1,]
    } else {
      AllChars <-rbind(AllChars, GSEFile[1,])
    }
    GSEFile <- GSEFile[-1,]
    count = count + 1
  }
  colnames(GPLFile)[1] <- "ID_REF"
  colnames(GPLFile)[2] <- "Gene Symbol"
  colnames(GSEFile)[1] <- "ID_REF"
  #merge by IDREF
  GPLandGSE <- merge(GPLFile, GSEFile)
  #removes probe column.gene duplicates will now be merged.
  GPLandGSE <- GPLandGSE[,-1]
  IDS <- GSEFile[1,]
  colnames(IDS) [1] <- "V1"
  #to change into numeric, have to split first column from all others
  #change others to numeric and combine again
  GPLandGSE2 <- data.frame(lapply(GPLandGSE[,2:ncol(GPLandGSE)], as.numeric), stringsAsFactors=FALSE)
  GPLandGSE3 <- cbind(GPLandGSE[,1], GPLandGSE2)
  colnames(GPLandGSE3)[1] <- "V1"
  index <- which(duplicated(GPLandGSE3[,1]))
  genesAggregated <- GPLandGSE3[-index, ]
  #time to aggregate
  #genesAggregated <- aggregate(GPLandGSE3, by = list(GPLandGSE3[,1]), FUN = mean)
  #genesAggregated <-  genesAggregated[,-2]
  #to allow the rows to be bound
  #colnames(genesAggregated)[1] <- "V1"
  RNASeqExpressionTable <- rbind(AllChars, IDS, genesAggregated)
}
#get GeneListFile file (if any)
#geneList <- read.delim("Documents/College/BioInformatics_Research/RNASeqPipeline/FakeGeneList.txt", header=F, stringsAsFactors=FALSE)
#geneList <- read.delim("bim_gene_list_1_sample.txt", header=F, stringsAsFactors=FALSE)
if(file.exists(paste("/var/www/magnet/uploads/gene_list1_0_", FindID, ".txt", sep = ""))) {
  geneList <- read.delim(paste("/var/www/magnet/uploads/gene_list1_0_", FindID, ".txt", sep = ""), header=F, stringsAsFactors=FALSE)
  cat(paste(toString(Sys.time()), " - INFO - ", "Reading gene list.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  geneList[,1] <- toupper(geneList[,1])
  
}
#for testing charLists
#charList <- read.delim("Documents/College/BioInformatics_Research/RNASeqPipeline/FakeCharacteristics.txt", header=F, stringsAsFactors=FALSE)

#get Char columns
cat(paste(toString(Sys.time()), " - INFO - ", "Removing undesired characteristics.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
columns <- table[4]
str <- columns[1,]
charList_notMatrix <- strsplit(str, ",")
charList <- matrix(unlist(charList_notMatrix), byrow = TRUE)
if(length(charList) == 0) {
  charList <- NA
}

#Vishal add Check 3
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 3.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")


#get all desired columns. (desired columns are in SQL Table)
MinusCharRNASeqExpressionTable <- RNASeqExpressionTable[,1] 
if(is.na(charList) != TRUE) {
  for(i in 1:nrow(charList)) {
    MinusCharRNASeqExpressionTable = cbind(MinusCharRNASeqExpressionTable, RNASeqExpressionTable[,(as.numeric(charList[i])+1)]);
  }
}else{
  for(i in 2:ncol(RNASeqExpressionTable)) {
    MinusCharRNASeqExpressionTable = cbind(MinusCharRNASeqExpressionTable, RNASeqExpressionTable[,i]);
  }
}

#looking for start of the table
for(i in 1:nrow(MinusCharRNASeqExpressionTable)) {
  if(MinusCharRNASeqExpressionTable[i,1] == "ID_REF" || MinusCharRNASeqExpressionTable[i,1] == "IDREF"){
    startIndex = i+1;
    break;
  }
}

#makes a new table with only the desired genes
FinalRNASeqExpressionTable <- MinusCharRNASeqExpressionTable[1:startIndex-1,]

#goes through and adds all rows of genes that are needed
if(exists("geneList") == TRUE) {
  cat(paste(toString(Sys.time()), " - INFO - ", "Removing undesired genes.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  for(i in startIndex:nrow(MinusCharRNASeqExpressionTable)) {
    for (j in 1:nrow(geneList)) {
      if(toupper(MinusCharRNASeqExpressionTable[i,1]) == geneList[j,1]) {
        FinalRNASeqExpressionTable <- rbind(FinalRNASeqExpressionTable, MinusCharRNASeqExpressionTable[i,])
      }
    }
  }
  #if gene llist doesnt exist, just use all the genes
} else {
  FinalRNASeqExpressionTable <- MinusCharRNASeqExpressionTable
}

#changing matrix to work with cor function
FinalRNASeqExpressionTable <- FinalRNASeqExpressionTable[startIndex:nrow(FinalRNASeqExpressionTable),]
FinalRNASeqExpressionTable[,1] <- toupper(FinalRNASeqExpressionTable[,1])

FinalRNASeqExpressionTable2 <- FinalRNASeqExpressionTable[,-1]
FinalRNASeqExpressionTable2 <- matrix(as.numeric(FinalRNASeqExpressionTable2), nrow = nrow(FinalRNASeqExpressionTable2))
rownames(FinalRNASeqExpressionTable2) <- FinalRNASeqExpressionTable[,1]
#use later if needed
#test<- FinalRNASeqExpressionTable2[complete.cases(FinalRNASeqExpressionTable2),]
#if(ncol(test) != 0) {
#  FinalRNASeqExpressionTable2 <- test
#}
FinalRNASeqExpressionTable2 <- t(FinalRNASeqExpressionTable2)

localization <- table[16]
coclustering <- table[17]
literature <- table[18]
#check if we need to calculate correlation values for all of them
PPIN = FALSE
if(localization != 0 || literature != 0 || coclustering != 0) { 
  PPIN = TRUE
}

if(exists("geneList") == TRUE & PPIN == TRUE) {
  #also need to calculate correlation values for all genes. if there is no gene list,
  #FinalRNASeqExpressionTable2 will already have all genes. no need to calculate twice
  AllGeneCorrelations <- MinusCharRNASeqExpressionTable
  AllGeneCorrelations <- AllGeneCorrelations[startIndex:nrow(AllGeneCorrelations),]
  AllGeneCorrelations2 <- AllGeneCorrelations[,-1]
  
  AllGeneCorrelations2 <- matrix(as.numeric(AllGeneCorrelations2), nrow = nrow(AllGeneCorrelations2))
  rownames(AllGeneCorrelations2) <- toupper(AllGeneCorrelations[,1])
  AllGeneCorrelations2 <- t(AllGeneCorrelations2)
}


#corpij <- function(i,j,data) {
  #print(i) 
  #round(cor(data[,i],data[,j], method="pearson",use="na.or.complete"),digits=2)
#}
#corpij2 <- function(i,j,data) {round(cor(data[,i],data[,j], method="spearman",use="na.or.complete"),digits=2)}

spearman <- table[21]
cat(paste(toString(Sys.time()), " - INFO - ", "Calculating correlations.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

#Vishal add Check 4
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 4.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")


if(spearman == 0) {
#   FinalCorrelation <- cor(FinalRNASeqExpressionTable2[,1:(ncol(FinalRNASeqExpressionTable2)/2)], method="pearson",use="na.or.complete")
#   FinalCorrelation2 <- cor(FinalRNASeqExpressionTable2[,(ncol(FinalRNASeqExpressionTable2)/2):ncol(FinalRNASeqExpressionTable2)], method="pearson",use="na.or.complete")
#   FinalCorrelation <- cbind(FinalCorrelation, FinalCorrelation2)
  #avoid calculating twice. cor function takes a while on large dataset.
  FinalCorrelation <- cor(FinalRNASeqExpressionTable2[,1:ncol(FinalRNASeqExpressionTable2)], method="pearson",use="pairwise.complete.obs")
  #corp <- Vectorize(corpij, vectorize.args=list("i","j"))
  #FinalCorrelation<-outer(1:ncol(FinalRNASeqExpressionTable2),1:ncol(FinalRNASeqExpressionTable2),corp,data=FinalRNASeqExpressionTable2)
  #rownames(FinalCorrelation) <- FinalRNASeqExpressionTable[,1]
  #colnames(FinalCorrelation) <- FinalRNASeqExpressionTable[,1]
  if(exists("geneList") == TRUE & PPIN == TRUE) {
    FinalAllGenesCorrelation <- cor(AllGeneCorrelations2[,1:ncol(AllGeneCorrelations2)], method="pearson",use="pairwise.complete.obs")
    #corp <- Vectorize(corpij, vectorize.args=list("i","j"))
    #FinalAllGenesCorrelation<-outer(1:ncol(AllGeneCorrelations2),1:ncol(AllGeneCorrelations2),corp,data=AllGeneCorrelations2)
    #rownames(FinalAllGenesCorrelation) <- AllGeneCorrelations[,1]
    #colnames(FinalAllGenesCorrelation) <- AllGeneCorrelations[,1]
  } else {
    FinalAllGenesCorrelation <- FinalCorrelation
  }
} else {
  FinalCorrelation <- cor(FinalRNASeqExpressionTable2[,1:ncol(FinalRNASeqExpressionTable2)], method="spearman",use="pairwise.complete.obs")
  #avoid calculating twice. cor function takes a while on large dataset.
  #corp <- Vectorize(corpij2, vectorize.args=list("i","j"))
  #FinalCorrelation<-outer(1:ncol(FinalRNASeqExpressionTable2),1:ncol(FinalRNASeqExpressionTable2),corp,data=FinalRNASeqExpressionTable2)
  #rownames(FinalCorrelation) <- FinalRNASeqExpressionTable[,1]
  #colnames(FinalCorrelation) <- FinalRNASeqExpressionTable[,1]
  if(exists("geneList") == TRUE & PPIN == TRUE) {
    FinalAllGenesCorrelation <- cor(AllGeneCorrelations2[,1:ncol(AllGeneCorrelations2)], method="spearman",use="pairwise.complete.obs")
    #corp <- Vectorize(corpij2, vectorize.args=list("i","j"))
    #FinalAllGenesCorrelation<-outer(1:ncol(AllGeneCorrelations2),1:ncol(AllGeneCorrelations2),corp,data=AllGeneCorrelations2)
    #rownames(FinalAllGenesCorrelation) <- AllGeneCorrelations[,1]
    #colnames(FinalAllGenesCorrelation) <- AllGeneCorrelations[,1]
    } else {
    FinalAllGenesCorrelation <- FinalCorrelation
  }
}

#Vishal add Check 5
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 5.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

#PPIN selected.
if(PPIN == TRUE) {
  cat(paste(toString(Sys.time()), " - INFO - ", "Querying interactions database.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  taxID <- table[7]
  query <- "SELECT GENENAME_A, GENENAME_B, GOLDEN"
  if(coclustering == 1) {
    query <- paste(query, ", CLUSTERING", sep  = "")
  }
  if(literature == 1) {
    query <- paste(query, ", EXPERIMENT", sep  = "")
  }
  if(localization == 1) {
    query <- paste(query, ", LOCATION", sep  = "")
  }
  #depending on different tax IDs, have to select different protein protein interactions
  if(taxID == "" || taxID == "null") {
    query <- paste(query, " FROM IntAct", sep  = "")
    res10 <- dbSendQuery(con, query)
  } else if(taxID == "human") {
    query <- paste(query, " FROM IntAct WHERE ORGANISM_A_ID LIKE 'taxid:9606(Human)%' AND ORGANISM_B_ID LIKE 'taxid:9606(Human)%'", sep  = "")
    res10 <- dbSendQuery(con, query)
  } else if(taxID == "mammalian") {
    query <- paste(query, " FROM IntAct WHERE ((ORGANISM_A_ID LIKE 'taxid:9606(Human)%') OR (ORGANISM_A_ID LIKE 'taxid:10090(Mouse)%') OR (ORGANISM_A_ID LIKE 'taxid:10116(Rat)%')) AND ((ORGANISM_B_ID LIKE 'taxid:9606(Human)%') OR (ORGANISM_B_ID LIKE 'taxid:10090(Mouse)%') OR (ORGANISM_B_ID LIKE 'taxid:10116(Rat)%')", sep = "")
    res10 <- dbSendQuery(con, query)
  } else if(taxID == "yeast") {
    query <- paste(query, "FROM IntAct WHERE ((ORGANISM_A_ID LIKE 'taxid:559292(Saccharomyces cerevisiae (strain ATCC 204508 / S288c))%') OR (ORGANISM_A_ID LIKE 'taxid:4932(Baker''s yeast)%')) AND ((ORGANISM_B_ID LIKE 'taxid:559292(Saccharomyces cerevisiae (strain ATCC 204508 / S288c))%') OR (ORGANISM_B_ID LIKE 'taxid:4932(Baker''s yeast)%'))", sep  = "")
    res10 <- dbSendQuery(con, query)
  } else { #assume all to avoid errors
    query <- paste(query, " FROM IntAct", sep  = "")
    res10 <- dbSendQuery(con, query)
  }
  tableAll <- fetch(res10, n = -1)
  index <- which(duplicated(tableAll))
  tableAll <- tableAll[-index, ]
  
  #get number of iterations  
  numIts <- table[22]
  if(is.na(numIts) || numIts == "") {
    numIts = 100
  }
  numIts <-  as.numeric(unlist(numIts))
  
  #add row for expression to table all, put in all expression values
  cat(paste(toString(Sys.time()), " - INFO - ", "Populating table with calculated correlation values. This may take a while.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  tableAll[,(ncol(tableAll) + 1)] <- NA
  colnames(tableAll)[ncol(tableAll)] <- "EXPRESSION"
  tableAll[,1] <- toupper(tableAll[,1])
  tableAll[,2] <- toupper(tableAll[,2])
  for(j in 1:nrow(tableAll)) { 
    if((tableAll[j,1] %in% row.names(FinalAllGenesCorrelation)) & (tableAll[j,2] %in% row.names(FinalAllGenesCorrelation))) {
      tryCatch({
        print(j)
        tableAll[j,ncol(tableAll)] <- FinalAllGenesCorrelation[tableAll[j,1], tableAll[j,2]]  
      }, error=function(e){})
    }
  }
  tableAll <- tableAll[complete.cases(tableAll),]
  
  cat(paste(toString(Sys.time()), " - INFO - ", "Generating beta values.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  for(i in 1:numIts) {
    #500 random goldens
    table500Golden <- tableAll[sample(which(tableAll$GOLDEN == 1),500),]
    #table500Golden <- data.frame(matrix(unlist(table500Golden), nrow = nrow(table500Golden)))
    #500 normals
    table500Normal <- tableAll[sample(which(tableAll$GOLDEN == 0),500),]
    #table500Normal <- data.frame(matrix(unlist(table500Normal), nrow = nrow(table500Normal)))
    #change column names and then make one whole table
    tableGoldenAndNormal<-rbind(table500Golden, table500Normal)
    #colnames(tableGoldenAndNormal)[3] <- "GOLDEN"
    #ugliest code I think I have ever written. Only doing it bc in a time crunch. If editing - FIX THIS FIRST.
    if(coclustering == 1 & literature == 0 & localization == 0) {
      colnames(tableGoldenAndNormal)[4] <- "CLUSTERING"
      temp <- glm( "GOLDEN ~ as.numeric(CLUSTERING)+as.numeric(EXPRESSION)" , data = tableGoldenAndNormal, family = binomial(logit))
    }else if(coclustering == 0 & literature == 1 & localization == 0) {
      colnames(tableGoldenAndNormal)[4] <- "EXPERIMENT"
      temp <- glm( "GOLDEN ~ as.numeric(EXPERIMENT)+as.numeric(EXPRESSION)" , data = tableGoldenAndNormal, family = binomial(logit))
    }else if(coclustering == 0 & literature == 0 & localization == 1) {
      colnames(tableGoldenAndNormal)[4] <- "LOCATION"
      temp <- glm( "GOLDEN ~ as.numeric(LOCATION)+as.numeric(EXPRESSION)" , data = tableGoldenAndNormal, family = binomial(logit))
    }else if(coclustering == 1 & literature == 1 & localization == 0) {
      colnames(tableGoldenAndNormal)[4] <- "CLUSTERING"
      colnames(tableGoldenAndNormal)[5] <- "EXPERIMENT"
      temp <- glm( "GOLDEN ~ as.numeric(CLUSTERING)+as.numeric(EXPERIMENT)+as.numeric(EXPRESSION)" , data = tableGoldenAndNormal, family = binomial(logit))
    }else if(coclustering == 1 & literature == 0 & localization == 1) {
      colnames(tableGoldenAndNormal)[4] <- "CLUSTERING"
      colnames(tableGoldenAndNormal)[5] <- "LOCATION"
      temp <- glm( "GOLDEN ~ as.numeric(CLUSTERING)+as.numeric(LOCATION)+as.numeric(EXPRESSION)" , data = tableGoldenAndNormal, family = binomial(logit))
    }else if(coclustering == 0 & literature == 1 & localization == 1) {
      colnames(tableGoldenAndNormal)[4] <- "EXPERIMENT"
      colnames(tableGoldenAndNormal)[5] <- "LOCATION"
      temp <- glm( "GOLDEN ~ as.numeric(EXPERIMENT)+as.numeric(LOCATION)+as.numeric(EXPRESSION)" , data = tableGoldenAndNormal, family = binomial(logit))
    }else {#ALL
      colnames(tableGoldenAndNormal)[4] <- "CLUSTERING"
      colnames(tableGoldenAndNormal)[5] <- "EXPERIMENT"
      colnames(tableGoldenAndNormal)[6] <- "LOCATION"
      temp <- glm( "GOLDEN ~ as.numeric(CLUSTERING)+as.numeric(EXPERIMENT)+as.numeric(LOCATION)+as.numeric(EXPRESSION)" , data = tableGoldenAndNormal, family = binomial(logit))
    }
    if(i==1){
      coef<-temp$coefficients
    }else{
      coef<-temp$coefficients+coef
    }
  }
  coef <- coef/numIts
  tableAll$PROBABILITY <- -1 # or 0
  temp$coefficients<- coef
  cat(paste(toString(Sys.time()), " - INFO - ", "Calculating probabilities", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  tableAll$PROBABILITY <- predict(temp, tableAll,type = "response")
  
  #need to upload betas
  dbSendQuery(con, paste("DELETE FROM process_beta WHERE request_id =", FindID))
  hold <- summary(temp)$coefficients
  size = nrow(hold)
  for(i in 1:size) {
    toAdd <- data.frame(hold[i,])
    toAdd <- t(toAdd)
    p_value = toAdd[4]
    p_value <- as.character(p_value)
    e_loc = gregexpr("e", p_value)
    e_loc = e_loc[[1]][1]
    if(e_loc!= -1) {
      p_value <- paste(substring(p_value, 1,5), substring(p_value, e_loc,nchar(p_value)), sep = "")
      p_value <- gsub("-", "", p_value)
    }
    else{
      p_value <- as.numeric(p_value)
      p_value <- signif(p_value, digits = 4)
    }
    toAdd <- toAdd[-4]
    toAdd <- as.numeric(formatC(toAdd, format="f", digits=3))
    toAdd <- t(toAdd)
    name = rownames(hold)[i]
    name =  strsplit(name, "(", fixed = TRUE)[[1]][2]
    name = substring(name, first = 1, last = nchar(name)-1)
    name <- paste(substring(name, first = 1, last = 1), tolower(substring(name, first = 2, last = nchar(name))), sep = "")
    toAdd <- cbind(FindID, toAdd, p_value , name)
    toAdd <- data.frame(toAdd, stringsAsFactors = FALSE)
    dbWriteTable(con, value = toAdd, field.types=list("request_id" = toAdd[,1], "estimate" = toAdd[,2], "std_error" = toAdd[,3], "z_value" = toAdd[,4], "p_value" = toAdd[,5], "row_name" = toAdd[,6]), name = "process_beta", append = TRUE, row.names=FALSE) 
  }
}
#Vishal add Check 6
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 6.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
#read in max/min
min <- as.numeric(table[10])
max <- as.numeric(table[9])

#read in ID, make FileNames
label <- table[3]
if(PPIN == FALSE) {
  fileName = paste(FindID, "_correlation", sep= "")
} else if(PPIN == TRUE) {
  fileName = paste(FindID, "_interaction", sep= "")
}

fileNameEDA = paste("/var/www/magnet/download/", fileName, ".eda", sep = "")
fileNameSIF = paste("/var/www/magnet/download/", fileName, ".sif", sep = "")
fileNameTXT = paste("/var/www/magnet/download/", fileName, ".txt", sep = "")
fileNameHIST = paste("/var/www/magnet/download/", fileName,"_hist",".pdf", sep = "")
fileNameFORMATTED = paste("/var/www/magnet/download/", fileName,"_formatted",".txt", sep = "")
if (file.exists(fileNameEDA)) {
  file.remove(fileNameEDA)
}
if (file.exists(fileNameSIF)) {
  file.remove(fileNameSIF)
}
if (file.exists(fileNameTXT)) {
  file.remove(fileNameTXT)
}
if (file.exists(fileNameHIST)) {
  file.remove(fileNameHIST)
}
if (file.exists(fileNameFORMATTED)) {
  file.remove(fileNameFORMATTED)
}

cat(paste(toString(Sys.time()), " - INFO - ", "Creating files.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}
#Vishal add Check 7
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 7.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

#write .txt file, SIF and EDA files
if(PPIN == FALSE) {
   #write .txt file
   write.table(format(FinalCorrelation, digits=3), file = fileNameTXT, sep = "\t", quote = FALSE, na = "")#, append = TRUE)
}

#Vishal add Check 7
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 7-2.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

if(PPIN == FALSE) {
  require(reshape2)
  FinalCorrelation[lower.tri(FinalCorrelation,diag=TRUE)] <- NA
  FinalCorrelation[ (FinalCorrelation < max) &  (FinalCorrelation > min)] <-NA

  siffile<- melt(FinalCorrelation,na.rm = TRUE)

  if(nrow(siffile)==0){
  cat(paste(toString(Sys.time()), " - INFO - ", "There are no correaltions in this range (Max= ",max,", Min=",min,"). Please re-submit with lower cut off values.",sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")  
 stop(structure(simpleError("There are no correlations."),class=c("specialError","error","condition")) ) 		
  }



  siffile$pp<-"pp"    
    siffile$pppp <-"(pp)"    
    siffile$eq   <-"="
  write.table((siffile[,c(1,4,2)]), file=fileNameSIF,sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE )
  write.table(format(siffile[,c(1,5,2,6,3)],digits=3), file=fileNameEDA,sep="\t" ,quote=FALSE,row.names = FALSE,col.names = FALSE )
  
 # counter<-1
 # n_size<-(nrow(FinalCorrelation)*ncol(FinalCorrelation))
 # # Use preallocated vectors
 # value <- numeric(n_size)
 # genename_a <- character(n_size)
 # 2D genename_b <- character(n_size)
 # # dbSendQuery(con, paste("DELETE FROM process_correlation WHERE request_id =", FindID))
 # for(j in 1:(nrow(FinalCorrelation) - 1)) {
 #   for(k in (j+1):ncol(FinalCorrelation)) {
 #     if(is.na(FinalCorrelation[j,k]) == FALSE 
 #        & (FinalCorrelation[j,k] >=max | FinalCorrelation[j,k] <= min)) {
 #       genename_a[counter]<-colnames(FinalCorrelation)[k]
 #       genename_b[counter]<-row.names(FinalCorrelation)[j]
 #       value[counter]<-specifydecimal(FinalCorrelation[j,k],3)
 #       counter<- counter+1;
 #     }
 #   }
 # }
 # siffile<-data.frame(genename_a, "pp", genename_b,stringsAsFactors=FALSE)
 # write.table(siffile[1:(counter-1),],file=fileNameSIF,sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE )
 # siffile<-data.frame(genename_a, "(pp)", genename_b,"=",value,stringsAsFactors=FALSE)
 # write.table(siffile[1:(counter-1),],file=fileNameEDA,sep="\t" ,quote=FALSE,row.names = FALSE,col.names = FALSE )

} else { #PPIN
  dbSendQuery(con, paste("DELETE FROM process_interaction WHERE request_id =", FindID))
  for(i in 1:nrow(tableAll)) {
    if((tableAll[i,1] %in% row.names(FinalCorrelation)) & (tableAll[i,2] %in% row.names(FinalCorrelation))) {
      checkDuplicate <- dbGetQuery(con, paste("SELECT * FROM process_interaction WHERE request_id = ", FindID, " AND genename_a = '", tableAll[i,2], "'AND genename_b = '", tableAll[i,1], "'", sep = ""))
      if(nrow(checkDuplicate) == 0) {
        #get all values to upload to sql table
        pub_ids <- dbGetQuery(con, paste("SELECT PUB_ID, INTERACTOR_1_HPRD, INT_ID, DB, DB, PID_A, PID_B FROM IntAct WHERE (genename_a = '", tableAll[i,1], "'AND genename_b = '", tableAll[i,2], "')", "OR", "(genename_a = '", tableAll[i,2], "'AND genename_b = '", tableAll[i,1], "')" ,sep = ""))
        #pub_ids <- dbGetQuery(con, paste("SELECT PUB_ID, INTERACTOR_1_HPRD, INT_ID, DB, DB, PID_A, PID_B FROM IntAct WHERE (genename_a = '", "CSNK2B", "'AND genename_b = '", "CSNK2A1", "')", "OR", "(genename_a = '", "CSNK2A1", "'AND genename_b = '", "CSNK2B", "')" ,sep = ""))
        most_common <- mode(pub_ids[,7])
        index = which(pub_ids$PID_A != most_common & pub_ids$PID_B != most_common)
        if(is.na(index[1]) == FALSE) {
          pub_ids <- pub_ids[-index,]
        }
        pubMedString <- ""
        intActString <- ""
        IntAct = 0
        HPRD = 0
        HPRD_ID = -1
        for(j in 1:nrow(pub_ids)) {
          pubMedList <- strsplit(pub_ids[j,1], "pubmed:", fixed = TRUE)[[1]][2]
          if(is.na(pubMedList) == FALSE) {
            pubMedList2 <- strsplit(pubMedList, "|", fixed = TRUE)[[1]][1]
            if( pubMedString == "") {
              pubMedString <- paste(pubMedString, pubMedList2)
            } else {
              pubMedString <- paste(pubMedString, ", ", pubMedList2)
            }
          } else{
            if( pubMedString == "") {
              pubMedString <- paste(pubMedString, pub_ids[j,1])
            } else {
              pubMedString <- paste(pubMedString, ",",pub_ids[j,1])
            }
          }
          intActList <- strsplit(pub_ids[j,3], "intact:", fixed = TRUE)[[1]][2]
          if(is.na(intActList) == FALSE) {
            intActList2 <- strsplit(intActList, "|", fixed = TRUE)[[1]][1]
            if( intActString == "") {
              intActString <- paste(intActString, intActList2)
            } else {
              intActString <- paste(intActString, ", ", intActList2)
            }
          } else{
            if(intActString == "") {
              intActString <- paste(intActString, pub_ids[j,3])
            } else {
              intActString <- paste(intActString, ",",pub_ids[j,3])
            }
          }
          if(pub_ids[j,4] == "IntAct") {
            IntAct = 1
          }else{
            HPRD = 1
            HPRD_ID = pub_ids[j,2]
            HPRD_ID = formatC(HPRD_ID, width = 5, format = "d", flag = "0")
          }
        }
        rowToWrite <- data.frame(a = FindID, b= tableAll[i,1], c= tableAll[i,2], d=pubMedString, e=HPRD_ID, f=intActString, g = 0, h = 0,i = specifydecimal(tableAll[i,]$EXPRESSION,3), j = 0, k = specifydecimal(tableAll[i,]$PROBABILITY,3), l = IntAct, m = HPRD, n = pub_ids[1,6], o = pub_ids[1,7])
        if(localization == 1) {
          rowToWrite$g <- tableAll[i,]$LOCATION
        }
        if(coclustering == 1) {
          rowToWrite$h <- specifydecimal(tableAll[i,]$CLUSTERING,3)
        }
        if(literature == 1) {
          rowToWrite$j <- tableAll[i,]$EXPERIMENT
        }
        dbWriteTable(con, value = rowToWrite, field.types=list("request_id" = rowToWrite$a,"genename_a"=rowToWrite$b, "genename_b" =rowToWrite$c, "pubmed_id" = rowToWrite$d, "hprd_id" = rowToWrite$e, "intact_id" = rowToWrite$f,"location" = rowToWrite$g,"coclustering" = rowToWrite$h,"coexpression" = rowToWrite$i,"observations" = rowToWrite$j,"probability" = rowToWrite$k,"intact" = rowToWrite$l,"hprd" = rowToWrite$m,"pid_a" = rowToWrite$n,"pid_b" = rowToWrite$o), name = "process_interaction", append = TRUE, row.names=FALSE)    
        #for sif
        write(paste(tableAll[i,1],"pp",tableAll[i,2]), file=fileNameSIF, sep="\t", append = TRUE)
        #for eds
        write(paste(tableAll[i,1],"(pp)", tableAll[i,2], "=", specifydecimal(tableAll[i,ncol(tableAll)], 3)),file=fileNameEDA,sep="\t", append = TRUE)
      }
    }
  }
}
#Vishal add Check 8
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 8.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

#write Histogram
if(PPIN == FALSE) {
  for(i in 1:nrow(FinalCorrelation)) {
    FinalCorrelation[i,i] <- NA;
  }
  pdf(fileNameHIST)
  par(oma = c(2,2,2,2))
  hist(FinalCorrelation[1:(nrow(FinalCorrelation)),1:(ncol(FinalCorrelation))], main = "Histogram of Coexpression Correlations", xlim = c(-1,1),xlab = "Correlations", col = "lightblue")
  mtext(paste("Uploaded file: " , fileName, sep=""), outer= TRUE, 1 , col= "blue")
  dev.off()
} else {
  pdf(fileNameHIST)
  par(oma = c(2,2,2,2))
  hist(tableAll$PROBABILITY, main = "Histogram of All PPI Probabilites", xlab = "PPI Probabilites", col = "lightblue")
  mtext(paste("Uploaded file: " , fileName, sep=""), outer= TRUE, 1 , col= "blue")
  dev.off()
}

cat(paste(toString(Sys.time()), " - INFO - ", "Histogram is generated.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

#write .txt file
#if(PPIN == FALSE) {
#  for(i in 1:nrow(FinalCorrelation)) {
#    for(j in 1:ncol(FinalCorrelation)) {
#      if(i ==j) {
#        FinalCorrelation[i,j] <- 1.000
#      }# else if(is.na(FinalCorrelation)[i,j] == TRUE){
#       # FinalCorrelation[i,j] <- NA
#      #} else if(i > j) {
#      #  FinalCorrelation[i,j] <- NA
#      #} else if(FinalCorrelation[i,j] < max & FinalCorrelation[i,j] > min) {
#      #  FinalCorrelation[i,j] <- NA
#      #}
#    }
#  }
#  write.table(format(FinalCorrelation, digits=3), file = fileNameTXT, sep = "\t", quote = FALSE, na = "")#, append = TRUE)
#}

formattedMatrix = table[13]
if(formattedMatrix == 1) {
  write.table(format(RNASeqExpressionTable, digits=3), file = fileNameFORMATTED, sep = "\t", quote = FALSE, na = "")
}

cat(paste(toString(Sys.time()), " - INFO - ", "Finished job.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

email = table[2]
print(email)
#Vishal add Check 9
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 9.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

cat(paste(toString(Sys.time()), " - INFO - ", paste("View your results ", " <a href=\"http://magnet.case.edu/process/results?uid=", email,"$$", FindID, "\">here</a>", sep = "") , sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
dbSendQuery(con, paste("UPDATE process_processrequest SET status = 'processed' WHERE id = ", FindID, sep = ""))

#Vishal add Check 10
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 10.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

if(email != "") {
  sendmail(from = "magnet@case.edu", to = email, subject = paste("MAGNET v.2 Job Results (",label," - ",FindID, " [",randomID ,"]" , ")",sep=""), msg= paste("View your results at ", "http://magnet.case.edu/process/results?uid=", email,"$$", FindID, sep = ""),control=list(smtpServer="smtp.case.edu"))
}

#Vishal add Check 11
#cat(paste(toString(Sys.time()), " - INFO - ", "CHECK 11.", sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")

}, error=function(e){
  cat(paste(toString(Sys.time()), " - INFO - ", "Something went wrong, check the formatting of your files and your parameters. Contact gurkan {at} case {dot} edu for more help with your job" , sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  cat(paste(toString(Sys.time()), " - INFO - ", toString(e) , sep = ""), file = paste("/var/www/magnet/logs/",FindID, ".txt", sep= ""), append = TRUE, sep = "\n")
  if(email != "") {
    sendmail(from = "magnet@case.edu", to = email, subject = "MAGNET Job Results", msg= "Something went wrong, check the formatting of your files and your parameters. Contact gurkan {at} case {dot} edu for more help with your job",control=list(smtpServer="smtp.case.edu"))
  }
})
#END
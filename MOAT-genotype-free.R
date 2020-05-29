
#clean-up
rm(list=ls())

MOATversion <- "2-3"
Updated <- "28-05-2020"


#############################################################
#############################################################
##                                                           ##
###### MOAT DATA INPUT FROM BESPOKE PIPELINE CSV - GENOTYPE FREE ####################################################
##                                                           ##
#############################################################
#############################################################

#search criteria for run names
search <- "Test"

#directory for run data
scriptPath <- function() {
  getSrcDirectory(scriptPath);
}
wd <- paste0(scriptPath(),"/")
setwd (wd)

#MOATpanel used for filtering of SNPs in MOAT assay
MOATpanel <- read.csv(paste0(wd,"panel_template_MOAT.csv"))

###Find run data based on search. Results File has to begin with "Results"
files <- list.files(wd, pattern="^Results", recursive=T)
files <- files[grep(search,files)]
show(files)

#add all run data to list
SNPcsvALL <- lapply(files, function(x){
  read.csv(x, header=F)
})
lapply(SNPcsvALL, head)

#Create result folder for each RUN
for(run in 1:length(SNPcsvALL)){
  SNPcsv1 <- SNPcsvALL[[run]]
  resfol <- paste0(dirname(files[run]),"/")
  Worksheet1 <- as.character(SNPcsv1[2,2])
  setwd(paste0(wd,resfol))
  dir.create(paste0("MOAT_" , Worksheet1))
  setwd(paste0("MOAT_" , Worksheet1))
  cols <- NCOL(SNPcsv1)
  
  #save patient raw data as individual csv files and folders
  for (patient in 2:cols){
    ptdata <- cbind(as.matrix(SNPcsv1$V1))
    ptdata <- cbind(ptdata, as.matrix(SNPcsv1[,patient]))
    ID <- (SNPcsv1[1,patient])
    Worksheet <- SNPcsv1[2,patient]
    Episode <- SNPcsv1[3,patient]
    PtName <- SNPcsv1[4,patient]
    PtName <- gsub(" ","_",PtName)
    Type <- SNPcsv1[5,patient]
    dir.create(paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Results"))
    setwd(paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Results"))
    write.table(ptdata, file = ((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Raw_Data.csv"))), sep=",",  col.names=FALSE, row.names=FALSE)
    SNPcsv <- read.csv(paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Raw_Data.csv"), header = FALSE)

    #return to run directory for next sample in loop
    setwd(paste0(wd, resfol,"MOAT_" , Worksheet1))
  }
}


     #############################################################
     #############################################################
    ##                                                           ##
###### MOAT DATA INPUT FROM BESPOKE PIPELINE CSV - GENOTYPE FREE #########
    ##                                                           ##
     #############################################################
     #############################################################

setwd (wd)


###Find and select individual sample data
files <- list.files(wd, pattern="Raw_Data.csv", recursive=T)
files <- files[grep(search,files)]
show(files)

##Get all patient data files
SNPcsvALL <- lapply(files, function(x){
  read.csv(x, header=F)
})
lapply(SNPcsvALL, head)

#Master loop for all runs and samples
for (patient in 1:length(files)){
  SNPcsv <- SNPcsvALL[[patient]]
  resfol <- paste0(dirname(files[patient]),"/")
  setwd(paste0(wd,resfol))
  
SNPcsvmat <- as.matrix(SNPcsv)
MOATpanelmat <- as.matrix(MOATpanel)

#Pulling sample/run details for folder and file names
ID <- SNPcsv[1,2]
Worksheet <- SNPcsv[2,2]
Episode <- SNPcsv[3,2]
PtName <- SNPcsv[4,2] #replace spaces with underscores
PtName <- gsub(" ","_",PtName)
Type <- SNPcsv[5,2]

resultsfolder <- ((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_MOAT_Results"))) #used in filenames


#Assigning SNP names from MOAT panel to list SNPnames
MOATpanellist <- as.list(MOATpanel)
MOATpanellist <- MOATpanellist[["Group"]]
SNPnames <- MOATpanellist[4:(length(MOATpanellist))]
SNPnames <- as.matrix(SNPnames)
SNPsearch <- SNPnames
SNPs <- length(SNPsearch)

SNPoutput <- matrix(0, nrow=SNPs, ncol=4, byrow=TRUE)
dimnames(SNPoutput) = list (SNPnames, c("A", "B", "Genotype", "MAF"))

#######################################Genotype SNPs######################################################################
# Genotype assigned AA or BB if frequency is over 86%, otherwise will be AB
# MOAT ANALYSIS of cfDNA will not be able to detect SNPs greater than 14% (i.e 100-86)
##########################################################################################################################

###start of loop#############################################################################################################

for(i in 1:SNPs){
  #search by list of SNPnames from MOAT panel
  SNPsearch <- SNPnames[i]
  #output to temp
  tmp <- SNPcsv[grep(SNPsearch, SNPcsv$V1), ]
  if (is.na(tmp[1,1])) {
    next
  }
  tmplst <- as.list(tmp[["V1"]])
  
  #removing "c" SNPs i.e. anything with a row 3
  if (length(tmplst) == 3) {
    tmp <- tmp[-3,]
  }
  #name as 1st column i.e. SNP name
  tmp2 <- data.frame (tmp, row.names = 1)
  #x <-> Y
  tmp3 <- t(tmp2)
  rownames(tmp3) <- c(SNPsearch)
  
  #Calculating allele frequencies
  Afrq <- as.numeric(tmp3[1,1])/(as.numeric(tmp3[1,1])+as.numeric(tmp3[1,2]))
  Bfrq <- as.numeric(tmp3[1,2])/(as.numeric(tmp3[1,1])+as.numeric(tmp3[1,2]))
  Fail <- is.nan(Afrq) 
  Afrq[is.nan(Afrq)] <- 0 # ?not actually working, still appearing as NaN in Processed_Data
  Bfrq[is.nan(Bfrq)] <- 0 
  
  #Assigning homozygous genotype assuming 86%
  if (Afrq > 0.86) {
    genotype <- "AA"
  } else if (Bfrq > 0.86) {
    genotype <- "BB"
  } else
    genotype <- "AB"
  
  #add genotype to table
  tmp4 <- cbind(tmp3, genotype )
  tmp4df <- data.frame(tmp4)
  
  # calculate MAF
  if (("AA" %in% tmp4) == FALSE) {
    MAF <- as.numeric(tmp4[1,1])/(as.numeric(tmp4[1,1])+as.numeric(tmp4[1,2]))*100
  } else 
    MAF <- as.numeric(tmp4[1,2])/(as.numeric(tmp4[1,1])+as.numeric(tmp4[1,2]))*100
  # add MAF to table
  tmp5 <- cbind(tmp4, MAF)
  
  #Final assign to SNPname data
  (assign(SNPsearch, tmp5))
  SNPoutput[i,] <- tmp5
}


#GENOTYPE_DATA: output processed data - all SNPs, genotype, MAFs
write.csv(SNPoutput, file = ((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Genotype_Data_allSNPs_","v",MOATversion,".csv"))))

#remove poor SNPs with count less than 6337
poorSNPsrow <- 0
for(i in 1:SNPs){
  if ((as.numeric(SNPoutput[i,1])+as.numeric(SNPoutput[i,2])) < 6255) {
    poorSNPsrow <- c(poorSNPsrow,i) 
  }
}
poorSNPsrow <- poorSNPsrow[-1]
poorSNPs <- SNPoutput[poorSNPsrow,]
lenpoorSNPs <- length(poorSNPs[,1])
write.csv(poorSNPs, file = ((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_lowcount_SNPs_","v",MOATversion,".csv"))))

#add good SNPs to file
goodSNPsrow <- 0
for(i in 1:SNPs){
  if ((as.numeric(SNPoutput[i,1])+as.numeric(SNPoutput[i,2])) >= 6255) {
    goodSNPsrow <- c(goodSNPsrow,i) 
  }
}
goodSNPsrow <- goodSNPsrow[-1]
goodSNPs <- SNPoutput[goodSNPsrow,]

#error catching if only 1 row fail else
if (length(goodSNPs) == 4) {
  lengoodSNPs <- 1
  write.csv(as.array(goodSNPs), file = ((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Genotype_Data_1SNP_","v",MOATversion,".csv"))))
  next
} else if (length(goodSNPs) == 0) {
  next
} else {
lengoodSNPs <- length(goodSNPs[,1])
}

#output all SNPs with high counts
write.csv(goodSNPs, file = ((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Genotype_Data_filtered_","v",MOATversion,".csv"))))

#output recipient informative SNPs (AA or BB)
#INFORMATIVE SNPS for analysis: Recipient AA or BB
AASNPs <- ""
BBSNPs <- ""
informativeSNPs <- (goodSNPs[- grep("AB", goodSNPs[,3]),])

write.csv(informativeSNPs, file = ((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Informative_SNPs_","v",MOATversion,".csv"))))


     #############################################################
     #############################################################
    ##                                                           ##
######            MOAT GENOTYPE FREE ANALYSIS                    #################################################
    ##                                                           ##
     #############################################################
     #############################################################

library(mixtools)

#################################  Y counting and fraction  ########################################################

#assign sex chromosome counts
SRYcounts <- SNPcsvmat[23,2]
ZFXcounts <- SNPcsvmat[21,2]
ZFYcounts <- SNPcsvmat[22,2]

#If ZFY present calculate Yfraction
#if Y fraction about 70% assuming male recipient otherwise assume female recipient
#If ZFY not present assume female recipient and female donor
#WARNING in low fraction patients there can be dropout of SRY and/or ZFY
if (as.numeric(ZFYcounts) > 0) { 
  Yfraction <- ((as.numeric(ZFYcounts))*2*100) / (as.numeric(ZFXcounts)+as.numeric(ZFYcounts))
  if (Yfraction > 70) {
    sex <- "XY/X?"
  } else if (Yfraction > 0) {
  sex <- "XX/XY"
  }
} else {
  sex <- "XX/XX"
  Yfraction <- 0
}

################################  SNP ANALYSIS  #################################################################
#using normalmixEM to model for expected SNP patterns based on raw SNP type counts
#################################################################################################################

#import recipient homozyous SNPs which can be used for analysis
datainput <- read.csv(((paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_Informative_SNPs_","v",MOATversion,".csv"))))
infSNPs <- datainput[,5]

#count informative SNPs
SNPcount <- length(infSNPs)

#used to estimate mu in modelling highest MAF estimate hom/diff SNP, half of this for hom/het SNPs
maxSNP <- max(infSNPs)

#repeat normalmixEM modelling
#set matrix for 500 loop of normalmixEM
iterations <- 500
lambdarep <- matrix(0, nrow=iterations, ncol=3, byrow=TRUE)
  dimnames(lambdarep) = list (c(), c("Background", "Het", "Hom"))
murep <- matrix(0, nrow=iterations, ncol=3, byrow=TRUE)
  dimnames(murep) = list (c(), c("Background", "Het", "Hom"))
sigmarep <- matrix(0, nrow=iterations, ncol=3, byrow=TRUE)
  dimnames(sigmarep) = list (c(), c("Background", "Het", "Hom"))
restartsrep <- matrix(0, nrow=iterations, ncol=1, byrow=TRUE)

#1st pass as 3 component model ie HOM/same, HOM/HET, HOM/diff#################################################################################################
#loop normalmixEM by the number of iterations to allow for averaging to reduce variation caused by random modelling
try(for(i in 1:iterations){
  results <- normalmixEM(infSNPs, k=3, mu = c(0.01, maxSNP/2, maxSNP), sigma= c(0.01,0.5,1), maxrestart = 1500, maxit = 1000) #mu (seq error rate based on Q30, het, maxSNP value assume HOM) #sigma start values improves analysis set to start around baseline level but not restrictive
  lambdarep[i,] <- (results[["lambda"]]) # proportions of each group
  murep[i,] <- (results[["mu"]]) # average of each group
  sigmarep[i,] <- (results[["sigma"]]) # SD of each group
  restartsrep[i,] <- (results[["restarts"]]) # record restarts as a measure of difficulty to model
}
)
#delete empty rows in cases of termination errors if the data cannot be modelled
#this may be gDNA samples, large differences in SNP MAFs, too few SNPs, trying to assing 1 SNP to a group
#some 1st Pass struggle to group and hit max retries therefore need to delete empty rows
(for (i in 1:iterations) {
  if  (lambdarep[i,2] == "0") {
    lambdarep <- lambdarep[-(i:iterations),]
    sigmarep <- sigmarep[-(i:iterations),]
    murep <- murep[-(i:iterations),]
    break
  }
  
}
)

#number of iterations for 1st pass
Nopass1 <- length(lambdarep)/3

#error catching if only 1 iteration
if (Nopass1 == 1) {
  lambdalist <- lambdarep
  mulist <- murep
  sigmalist <- sigmarep
} else {
  #averaging data from multiple iterations and add to new lists
  lambdalist <- c(mean(lambdarep[,1]), mean(lambdarep[,2]), mean(lambdarep[,3]))
  mulist <- c(mean(murep[,1]), mean(murep[,2]), mean(murep[,3]))
  sigmalist <- c(mean(sigmarep[,1]), mean(sigmarep[,2]), mean(sigmarep[,3]))
}

#error catching for failed normalMIX produce plot of all SNPs by A frequency
if (is.nan(lambdalist)) {
  AfreqallSNPs <- SNPoutput
  AfreqallSNPs <- cbind(SNPoutput,Afreq = 0)
  for (i in 1:(SNPs)) {
    if (AfreqallSNPs[i,3] == "AB") {
      AfreqallSNPs[i,5] <- AfreqallSNPs[i,4]
    } else 
      AfreqallSNPs[i,5] <- (as.numeric(AfreqallSNPs[i,1]) / (as.numeric(AfreqallSNPs[i,1]) + as.numeric(AfreqallSNPs[i,2]))*100)
  }
  setx1 <- rep(1,SNPcount)
  ymax <- (max(infSNPs))*1.05
  #All SNPs
  pdf(paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_plots_","v",MOATversion,".pdf"))
  par(mfrow=c(2,1), mar=c(0.8,3.5,3,1.5), mgp=c(2.3,1,0), oma=c(0,0,2,0))
  plot(jitter(rep(1,(SNPs)),4),as.numeric(AfreqallSNPs[,5]), col = rgb(0,0,0,0.6), xlim = c(0.1,2.1), ylim = c(0,100), cex.lab=0.8, xaxt = 'n', xlab="", ylab= "A frequency", pch=20, title(main="All Genotyped SNPs", cex.main=0.8, line=0.5))
  plot(jitter(setx1,1),infSNPs, col = rgb(0,0,0,0.7), pch=20, xaxt = 'n', xlim = c(0.1,2.2), ylim = c(0,ymax), ylab= "Minor allele %", cex.lab=0.8, title(main="Recipient Homozygous SNPs", cex.main=0.8, line=0.5))
  points(2,Yfraction, pch=20, col="blue", cex=0.8)
  text(1, mean(infSNPs), cex = 0.7, pos = 4, labels = paste0("n=", length(infSNPs)))
  abline(0.1168,0, col="red") #based on actual background SNPs from pooled spiked data (average+2SD)=99% - ***REVISIT WITH MORE DATA***
  abline(0.052,0, col="green") #based on modelled background SNPs from pooled spiked data (average+2SD)=99%
  legend("topleft", legend =c("Unallocated SNP","Background SNP", "Het SNP", "Hom-diff SNP", "ZFXY","LOQ", "LOQ - modelled") ,pch=c(20,20,20,20,20,3,3), col=c("black","red","orange","purple", "blue", "red", "green"), cex = 0.7)
  dev.off()
  next
}

#assigning SNPs
#prortion of SNPs * total number of SNPs
backSNPsa <- round(lambdalist[1]*SNPcount, 0) #lower mixing proportion = HOM/same = background
hetSNPs <- round(lambdalist[2]*SNPcount, 0) #middle mixing proportion = HOM/HET = informative for cell free fraction
homSNPs <- round(lambdalist[3]*SNPcount, 0) #upper mixing proportion = HOM/diff - informative for cell free fraction

sortSNPs <- sort(infSNPs, decreasing = FALSE)

#find the middle of HOM/same and HOM/HET to separate out HOM/same background SNPs
#works well with discrete groups of data, SDs are plotted and can be reviewed on pdf results
midlower <- (mulist[1]+mulist[2])/2

#anything below midlower classed as background SNPs
#count and assign SNPs to background and remainder are HOM/HET and HOM/diff
backSNPsb <- sum(infSNPs<midlower)
hethomSNPlist <- 0
hethomSNPlist <- sortSNPs[-(1:backSNPsb)] #remove background SNPs and send to variable
backSNPlist <- sortSNPs[-(backSNPsb+1:SNPcount)] # remove HOM/HET and HOM/diff SNPs and send to variable
sortedinformativeSNPs <- sort(informativeSNPs, decreasing = FALSE)
sortedinformativeSNPnames <- names(sortedinformativeSNPs[1:backSNPsb])
names(backSNPlist) <- sortedinformativeSNPnames
write.csv(backSNPlist, file = paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type,"_Background_SNPs_","v",MOATversion,".csv"))

#count informative SNPs
try (if (hethomSNPlist == 0) {
  SNPcount2 <- 0 #need to assign as 0 otherwise errors
} else {
  SNPcount2 <- length(hethomSNPlist)
}
)

#reiterate normalmixEM with hethomSNPlist
#2nd pass with only 2 component mix adding mean constraint for 1:2 average MAF i.e. het:hom counts #################

#set matrix for results
lambdarep2 <- matrix(0, nrow=iterations, ncol=2, byrow=TRUE)
dimnames(lambdarep2) = list (c(), c("Het", "Hom"))
murep2 <- matrix(0, nrow=iterations, ncol=2, byrow=TRUE)
dimnames(murep2) = list (c(), c("Het", "Hom"))
sigmarep2 <- matrix(0, nrow=iterations, ncol=2, byrow=TRUE)
dimnames(sigmarep2) = list (c(), c("Het", "Hom"))
restartsrep2 <- matrix(0, nrow=iterations, ncol=1, byrow=TRUE)

#running modelling pass2
try(for(i in 1:iterations){
  results2 <- normalmixEM(hethomSNPlist, k=2, mu = c(maxSNP/2, maxSNP), maxrestart = 1500, maxit = 1000, mean.constr = c("a","2a"), sd.constr = c("b","2b"))
  lambdarep2[i,] <- (results2[["lambda"]])
  murep2[i,] <- (results2[["mu"]])
  sigmarep2[i,] <- (results2[["sigma"]])
  restartsrep2[i,] <- (results2[["restarts"]])
}
)

for (i in 1:iterations) {
  if  (lambdarep2[i,2] == "0") {
    lambdarep2 <- lambdarep2[-(i:iterations),]
    sigmarep2 <- sigmarep2[-(i:iterations),]
    murep2 <- murep2[-(i:iterations),]
    break
  }
  
}

#number of competed iterations
Nopass2 <- length(lambdarep2)/2

#error need to catch 1 SNP pass
if (Nopass2 == 1) {
  lambdalist2 <- lambdarep2
  mulist2 <- murep2
  sigmalist2 <- sigmarep2
} else {
  #averaging data from 2nd pass iterations and add to new lists
  lambdalist2 <- c(mean(lambdarep2[,1]), mean(lambdarep2[,2]))
  mulist2 <- c(mean(murep2[,1]), mean(murep2[,2]))
  sigmalist2 <- c(mean(sigmarep2[,1]), mean(sigmarep2[,2]))
}

hetSNPlist2 <- 0
homSNPlist2 <- 0

hetSNPs2 <- round(lambdalist2[1]*SNPcount2, 0) #number of hetSNPs

#error handling for failed pass2 modelling
if (is.nan(hetSNPs2)) {
  hetSNPs2 <- (length(hethomSNPlist))
  }

homSNPs2 <- 0
homSNPs2 <- round(lambdalist2[2]*SNPcount2, 0) #number of homSNPs

#error handling for failed pass2 modelling
if (is.nan(homSNPs2)) {
  homSNPs2 <- 0
  error <- 0
} else 
error <- sum(length(hetSNPs2+homSNPs2))

if (hetSNPs2==0) {
  hetSNPlist2 <- 0
} else
  hetSNPlist2 <- hethomSNPlist[1:hetSNPs2]

#error handling for failed pass2 modelling
if (error == 0) {
  homSNPlist2 <- 0
} else 
homSNPlist2 <- hethomSNPlist[(hetSNPs2+1):SNPcount2]

ddcfDNA <- 0

error1 <- ""
error2 <- ""
error3 <- ""

#sexing of sample required for ddcfDNA calcualtion
if (sex == "XX/XY") {
  Male <- 1
} else {
  Male <- 0
}
if (error == 0) {
  error1 <- "not enough SNPs for Pass2 genotyping"
  error2 <- "cannot determine if Het/Hom"
  error3 <- "not enough SNPs for Pass2 genotyping; cannot determine if Het/Hom; cannot accurately determine ddcfDNA fraction"
  hetSNPs2 <- SNPcount2
  ddcfDNA <- (((mean(hethomSNPlist)*2*SNPcount2) + (Yfraction*Male)) / (SNPcount2+Male)) #if pass2 modelling fails assume informative SNPs are HOM/HET as most likely
  medcfDNA<- median(hethomSNPlist*2) #if pass2 modelling fails assume informative SNPs are HOM/HET as most likely
} else if (sex == "XX/XY") {
  ddcfDNA <- ((Yfraction)+(sum(hetSNPlist2)*2)+(sum(homSNPlist2)))/((SNPcount2)+1)
  medcfDNA <- median(c(Yfraction,(hetSNPlist2*2),(homSNPlist2)))
} else {
  ddcfDNA <- ((sum(hetSNPlist2)*2)+(sum(homSNPlist2)))/SNPcount2
  medcfDNA <- median(c((hetSNPlist2*2),(homSNPlist2)))
}
#error handling####################################
if (is.nan (ddcfDNA)){
  ddcfDNA <- 0
  medcfDNA <- 0
}

if (ddcfDNA < 0.1168) {
  error4 <- "WARNING - ddcfDNA is within background level"
} else 
  error4 <- ""

error5 <- ""
if (sum(as.numeric(SNPcsvmat[24:114,2]))/(SNPs) < 10000) {
  error5 <- "low SNP count"
}
###################################################

#Plots########################################################################################################################################################
pdf(paste0(Worksheet,"_",ID,"_", Episode, "_", PtName, "_", Type, "_plots_","v",MOATversion,".pdf"))
par(mfrow=c(2,2), mar=c(0.8,3.5,3,1.5), mgp=c(2.3,1,0), oma=c(0,0,2,0))

setx1 <- rep(1,SNPcount)

#setting Y-axis height limit based on highest SNP
if (sex == "XX/XY") {
  ymax <- (max(c(ddcfDNA, 0.1168, 0.052, (hetSNPlist2*2), homSNPlist2, Yfraction, infSNPs)))*1.05
} else 
  ymax <- (max(c(ddcfDNA, 0.1168, 0.052, (hetSNPlist2*2), homSNPlist2, infSNPs)))*1.05

#plot 1 - top left####################
AfreqallSNPs <- SNPoutput
AfreqallSNPs <- cbind(SNPoutput,Afreq = 0)
for (i in 1:(SNPs)) {
  if (AfreqallSNPs[i,3] == "AB") {
    AfreqallSNPs[i,5] <- AfreqallSNPs[i,4]
  } else 
    AfreqallSNPs[i,5] <- (as.numeric(AfreqallSNPs[i,1]) / (as.numeric(AfreqallSNPs[i,1]) + as.numeric(AfreqallSNPs[i,2]))*100)
}

#All SNPs
plot(jitter(rep(1,(SNPs)),4),as.numeric(AfreqallSNPs[,5]), col = rgb(0,0,0,0.6), xlim = c(0.1,2.1), ylim = c(0,100), cex.lab=0.8, xaxt = 'n', xlab="", ylab= "A frequency", pch=20, title(main="All Genotyped SNPs", cex.main=0.8, line=0.5))

#Plot2 - top right###################
#RAW SNP DATA
plot(jitter(setx1,1),infSNPs, col = rgb(0,0,0,0.7), pch=20, xaxt = 'n', xlim = c(0.1,2.2), ylim = c(0,ymax), ylab= "Minor allele %", cex.lab=0.8, title(main="Recipient Homozygous SNPs", cex.main=0.8, line=0.5))
points(2,Yfraction, pch=20, col="blue", cex=0.8)
text(1, mean(infSNPs), cex = 0.7, pos = 4, labels = paste0("n=", length(infSNPs)))
abline(0.1168,0, col="red") #based on actual background SNPs from pooled spiked data (average+2SD)=99% - ***REVISIT WITH MORE DATA***
abline(0.052,0, col="green") #based on modelled background SNPs from pooled spiked data (average+2SD)=99%
legend("topleft", legend =c("Unallocated SNP","Background SNP", "Het SNP", "Hom-diff SNP", "ZFXY","LOQ", "LOQ - modelled") ,pch=c(20,20,20,20,20,3,3), col=c("black","red","orange","purple", "blue", "red", "green"), cex = 0.7)

#plot 3 - bottom left
#sorted SNP data
#plot and colour code RecipientHom SNPs (Hom-same=backgound, het, hom-diff) #backgound from pass1 #het/hom from pass2
plot(jitter(rep(1,backSNPsb),3),backSNPlist, col = rgb(1,0,0,0.5), xlim = c(0.1,4.2),cex.lab=0.8, ylim = c(0,ymax), xaxt = 'n', ylab= "Minor allele %", pch=20, title(main=paste0("Assigned SNPs ","(Pass1 n=",Nopass1, ", Pass2 n=",Nopass2,")"), cex.main=0.8,line=0.5), xlab="")
try(points(jitter(rep(2,hetSNPs2),1),hetSNPlist2, col = rgb(1,0.5,0,0.65), pch=20))
try(points(jitter(rep(3,homSNPs2),1),homSNPlist2, col = rgb(0.55,0,0.65,0.65), pch=17))
points(4,Yfraction, pch=20, col="blue", cex=0.8)
points(1.2, mulist[1], pch="---", cex=1, col="black")
points(2.2, mulist2[1], pch="---", cex=1, col="black")
points(3.2, mulist2[2], pch="---", cex=1, col="black")
points(1.2, mulist[1]-sigmalist[1], pch="---", cex=1, col="green")
points(2.2, mulist2[1]-sigmalist2[1], pch="---", cex=1, col="green")
points(3.2, mulist2[2]-sigmalist2[2], pch="---", cex=1, col="green")
points(1.2, mulist[1]+sigmalist[1], pch="---", cex=1, col="green")
points(2.2, mulist2[1]+sigmalist2[1], pch="---", cex=1, col="green")
points(3.2, mulist2[2]+sigmalist2[2], pch="---", cex=1, col="green")
points(1.2, mulist[1]-2*sigmalist[1], pch="---", cex=1, col="red")
points(2.2, mulist2[1]-2*sigmalist2[1], pch="---", cex=1, col="red")
points(3.2, mulist2[2]-2*sigmalist2[2], pch="---", cex=1, col="red")
points(1.2, mulist[1]+2*sigmalist[1], pch="---", cex=1, col="red")
points(2.2, mulist2[1]+2*sigmalist2[1], pch="---", cex=1, col="red")
points(3.1, mulist2[2]+2*sigmalist2[2], pch="---", cex=1, col="red")
points(0.8, mean(backSNPlist), pch="---", cex=1, col="black")
points(1.8, mean(hetSNPlist2), pch="---", cex=1, col="black")
points(2.8, mean(homSNPlist2), pch="---", cex=1, col="black")
points(0.8, mean(backSNPlist)-sd(backSNPlist), pch="---", cex=1, col="green")
points(1.8, mean(hetSNPlist2)-sd(hetSNPlist2), pch="---", cex=1, col="green")
points(2.8, mean(homSNPlist2)-sd(homSNPlist2), pch="---", cex=1, col="green")
points(0.8, mean(backSNPlist)+sd(backSNPlist), pch="---", cex=1, col="green")
points(1.8, mean(hetSNPlist2)+sd(hetSNPlist2), pch="---", cex=1, col="green")
points(2.8, mean(homSNPlist2)+sd(homSNPlist2), pch="---", cex=1, col="green")
points(0.8, mean(backSNPlist)-2*sd(backSNPlist), pch="---", cex=1, col="red")
points(1.8, mean(hetSNPlist2)-2*sd(hetSNPlist2), pch="---", cex=1, col="red")
points(2.8, mean(homSNPlist2)-2*sd(homSNPlist2), pch="---", cex=1, col="red")
points(0.8, mean(backSNPlist)+2*sd(backSNPlist), pch="---", cex=1, col="red")
points(1.8, mean(hetSNPlist2)+2*sd(hetSNPlist2), pch="---", cex=1, col="red")
points(2.8, mean(homSNPlist2)+2*sd(homSNPlist2), pch="---", cex=1, col="red")
text(1,max(backSNPlist), cex=0.7, pos=3, labels = paste0("n=", backSNPsb))
text(2,max(hetSNPlist2), cex=0.7, pos=3, labels = paste0("n=", hetSNPs2))
text(3,max(homSNPlist2), cex=0.7, pos=3, labels = paste0("n=", homSNPs2))
abline(0.1168,0, col="red") #based on actual background SNPs from pooled spiked data (average+2SD)=99% - ***REVISIT WITH MORE DATA***
abline(0.052,0, col="green") #based on modelled background SNPs from pooled spiked data (average+2SD)=99% probably not useful
legend("topleft", legend =c("Unallocated SNP","Background SNP", "Het SNP", "Hom-diff SNP", "ZFXY","LOQ", "LOQ - modelled") ,pch=c(20,20,20,20,20,3,3), col=c("black","red","orange","purple", "blue", "red", "green"), cex = 0.7)

#plot4 - bottom right
#combined and averaged data
#plot and colour code RecipientHom SNPs (Hom-same=backgound, het, hom-diff)
plot(jitter(rep(2,backSNPsb),1),backSNPlist, col = rgb(1,0,0,0.3), xlim = c(0,4.2),cex.lab=0.8, ylim = c(0,ymax), xaxt = 'n', ylab= "ddcfDNA %", pch=20, title(main="ddcfDNA", cex.main=0.8,line=0.5),xlab="")
try(points(jitter(rep(2,hetSNPs2),1),hetSNPlist2*2, col = rgb(1,0.5,0,0.65), pch=20))
try(points(jitter(rep(2,homSNPs2),1),homSNPlist2, col = rgb(0.55,0,0.65,0.65), pch=20))
points(3,ddcfDNA, col= "green", pch=20)
points(3,medcfDNA, col= "blue", pch=20)
points(2,Yfraction, pch=20, col="blue", cex=0.8)
text(3.2,ddcfDNA, cex=0.7, pos=3, labels = paste0((format(round(ddcfDNA,2),nsmall = 2)),"% ddcfDNA (mean)"))
text(3.2,medcfDNA, cex=0.7, pos=1, labels = paste0((format(round(medcfDNA,2),nsmall = 2)),"% ddcfDNA (median)"))
abline(0.1168,0, col="red") #based on actual background SNPs from pooled spiked data (average+2SD)=99% - ***REVISIT WITH MORE DATA***
abline(0.052,0, col="green") #based on modelled background SNPs from pooled spiked data (average+2SD)=99% probably not useful
legend("topleft", legend =c("Unallocated SNP","Background SNP", "Het SNP", "Hom-diff SNP", "ZFXY","LOQ", "LOQ - modelled") ,pch=c(20,20,20,20,20,3,3), col=c("black","red","orange","purple", "blue", "red", "green"), cex = 0.7)

mtext((paste0(Episode, ": ", PtName, " ", "(",Type,")", "                " ,Worksheet,"       ")), outer=TRUE,  cex=1, adj=1)

dev.off()

######################################################################################################################################################################
#create analysis sheet with metrics


QClistnames <- c("",
                 "Recipient Hom SNPs","Pass1 (3 component)", "SNP Type", "SD", "Group Average", "Mixing Proportions", 
                 "",
                 "Pass2 (2 component)", "SNP Type", "SD", "Group Average", "Mixing Proportions",
                 "",
                 "Summary","SNP Counts","SRY Count",
                 "",
                 "",
                 "Background", "Het", "Hom","ZFXY","Average",
                 "",
                 "Conclusion:",
                 "",
                 "First Reader:","Second Reader:",
                 "",
                 "Version","Updated")

QCmat <- matrix("",0,3)

QCmat <- rbind(QCmat,c("","","")) #blank row

#Pass1 modelled data
QCmat <- rbind(QCmat,c(SNPcount,"Rejected SNPs",lenpoorSNPs))
QCmat <- rbind(QCmat,c("Iterations:", Nopass1, paste0("Restarts: ", sum(restartsrep))))
QCmat <- rbind(QCmat,c("Hom-same","Het","Hom-diff"))
QCmat <- rbind(QCmat,c(sigmalist))
QCmat <- rbind(QCmat,c(mulist))
QCmat <- rbind(QCmat,c(lambdalist))

QCmat <- rbind(QCmat,c("","","")) #blank row

#Pass2 Modelled data
QCmat <- rbind(QCmat,c("Iterations:",Nopass2, paste0("Restarts: ", sum(restartsrep2))))
QCmat <- rbind(QCmat,c("Hom-same","Het","Hom-diff"))
QCmat <- rbind(QCmat,c("",sigmalist2))
QCmat <- rbind(QCmat,c("",mulist2))
QCmat <- rbind(QCmat,c("",lambdalist2))

QCmat <- rbind(QCmat,c("","","")) #blank row

QCmat <- rbind(QCmat,c("Per SNP","Total",""))
QCmat <- rbind(QCmat,c((sum(as.numeric(SNPcsvmat[24:114,2])))/(SNPs),sum(as.numeric(SNPcsvmat[24:114,2])),error5)) #average per SNP #total SNPcounts
QCmat <- rbind(QCmat,c(SRYcounts,error5,""))

QCmat <- rbind(QCmat,c("","","")) #blank row

QCmat <- rbind(QCmat,c("% ddcfDNA (mean)","% ddcfDNA (median)","Informative SNPs")) 
QCmat <- rbind(QCmat,c(mean(backSNPlist),median(backSNPlist) , backSNPsb)) #average backgound #nunber of SNPs
QCmat <- rbind(QCmat,c(mean(hetSNPlist2)*2,median(hetSNPlist2)*2 , hetSNPs2))#hetSNP ddcfDNA #hetSNPs
QCmat <- rbind(QCmat,c(mean(homSNPlist2), median(homSNPlist2), homSNPs2)) #homSNP ddcfDNA #homSNPs
QCmat <- rbind(QCmat,c(Yfraction,Yfraction,error3))
QCmat <- rbind(QCmat,c(ddcfDNA,medcfDNA,sum(hetSNPs2,homSNPs2)))

QCmat <- rbind(QCmat,c("",error4,"")) #blank row

QCmat <- rbind(QCmat,c("","","")) #conclusion

QCmat <- rbind(QCmat,c("","","")) #blank row

QCmat <- rbind(QCmat,c("","","Date:")) #1st reader
QCmat <- rbind(QCmat,c("","","Date:")) #2nd reader

QCmat <- rbind(QCmat,c("","","")) #blank row

QCmat <- rbind(QCmat,c(MOATversion,"",""))
QCmat <- rbind(QCmat,c(Updated,"",""))

#rownames(QCmat) <- QClistnames
colnames(QCmat) <- c(resultsfolder,"","")
rownames(QCmat) <- QClistnames

write.csv(QCmat, file = paste0(resultsfolder,"_v",MOATversion,".csv"))

#?return to which folder for next repeat
setwd(wd)
}

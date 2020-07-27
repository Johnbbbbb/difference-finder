FirstSubstrateSet<- read.csv("input1.csv", stringsAsFactors=FALSE, header = FALSE)
Firstsubbackfreq<- read.csv("input2.csv", header=FALSE, stringsAsFactors=FALSE)
SubstrateHeader<-FirstSubstrateSet[1,]
FirstSubstrateSet<- FirstSubstrateSet[2:nrow(FirstSubstrateSet),]
if(nrow(Firstsubbackfreq[1,]>35)){
  if(grepl(pattern = "Properties", x=Firstsubbackfreq[1,22])){
    Firstsubbackfreq<-t(Firstsubbackfreq)
  }
}


SecondSubstrateSet<- read.csv("input3.csv", stringsAsFactors=FALSE, header = FALSE)
Secondsubbackfreq<- read.csv("input4.csv", header=FALSE, stringsAsFactors=FALSE)
SecondSubstrateSet<- SecondSubstrateSet[2:nrow(SecondSubstrateSet),]
if(nrow(Secondsubbackfreq[1,]>35)){
  if(grepl(pattern = "Properties", x=Secondsubbackfreq[1,22])){
    Secondsubbackfreq<-t(Secondsubbackfreq)
  }
}


ThirdSubstrateSet<- read.csv("input5.csv", stringsAsFactors=FALSE, header = FALSE)
Thirdsubbackfreq<- read.csv("input6.csv", header=FALSE, stringsAsFactors=FALSE)
ThirdSubstrateSet<- ThirdSubstrateSet[2:nrow(ThirdSubstrateSet),]
if(nrow(Thirdsubbackfreq[1,]>35)){
  if(grepl(pattern = "Properties", x=Thirdsubbackfreq[1,22])){
    Thirdsubbackfreq<-t(Thirdsubbackfreq)
  }
}

#the above three sections bring in the input files and ensure they are properly aligned.  The if statements align them if they are misaligned

First_unshared_motifs_table<-"R1 substrates.csv"
First_unshared_subbackfreq<-"R1 SBF.csv"

Second_unshared_motifs_table<-"R2 subs.csv"
Second_unshared_subbackfreq<-"R2 SBf.csv"

Third_unshared_motifs_table<-"R3 subs.csv"
Third_unshared_subbackfreq<-"R3 SBF.csv"

#the above 4 sections create the names of the output files that this tool can create




FirstxY<-rep("xY",times=nrow(FirstSubstrateSet))
FirstSubstrateSet[,11]<-FirstxY

SecondxY<-rep("xY",times=nrow(SecondSubstrateSet))
SecondSubstrateSet[,11]<-SecondxY

ThirdxY<-rep("xY",times=nrow(ThirdSubstrateSet))
ThirdSubstrateSet[,11]<-ThirdxY

#for each input file, mark the phospho-Tyrosine (which is always housed in column 11) 
#with an xY to denote that is it a phospho and not regular tyrosine


#currently the substrates are a dataframe with many values, I want to collapse each substrate into a single variable, 
#but if I simply use the paste() function then I lose any information about whether the substrate was truncated on the C or N terminal.  
#so the for loops below then are constructed so as to retain that information on how any substrate is truncated.

FTLwtmotifs=matrix(,nrow = nrow(FirstSubstrateSet),ncol=1)
FTLwtAccessionNumbers=matrix(data = Firstsubbackfreq[1,],ncol=1)
#create the vectors which will house the first set of substrates and the first set of accession numbers


for (i in 1:nrow(FirstSubstrateSet)){
  FTLwtletters<-FirstSubstrateSet[i,4:18]
  FTLwtletters<-FTLwtletters[FTLwtletters !="XXXXX"]
  FTLwtletters<-paste(FTLwtletters, sep="", collapse="")
  leftspaces<-c()
  rightspaces<-c()
  
  YYYmotif <- unlist(strsplit(FTLwtletters, split = ""))
  YYYposition <- match(x = "x", table = YYYmotif)
  #position itself tells me how much is to the left of that X by what it's number is.  x at position 4 tells me that there are
  #just 3 letters to the left of x
  
  YYYLettersToTheLeft <- YYYposition - 1
  YYYLettersToTheRight <- length(YYYmotif) - YYYposition - 1
  #how many letters to the right SHOULD just be length(motif)-position-1 if it's 5 long and x is at 3 then Y is at 4 and there is
  #just 1 spot to the right of Y so LettersToTheRight<-1 because 5-3-1=1
  
  
  
  if (YYYLettersToTheLeft < 7 | YYYLettersToTheRight < 7) {
    leftspaces<-rep(" ",times=(7-YYYLettersToTheLeft))
    rightspaces<-rep(" ",times=7-(YYYLettersToTheRight))
    #add blank spaces if the motif has less than 4 letters to the left/right
    motif<-c(leftspaces,YYYmotif,rightspaces)
    #save that motif, which is the Y and +/- 4 amino acids, including truncation
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    FTLwtletters<-motif
    FTLwtmotifs[i,1]<-FTLwtletters
  }
  
  if(YYYLettersToTheLeft>6 && YYYLettersToTheRight>6){
    motif<-YYYmotif
    #add blank spaces if the motif has less than 4 letters to the left/right
    motif<-c(leftspaces,YYYmotif,rightspaces)
    #save that motif, which is the Y and +/- 4 amino acids, including truncation
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    FTLwtletters<-motif
    FTLwtmotifs[i,1]<-FTLwtletters
  }
}

D835Ymotifs=matrix(,nrow = nrow(SecondSubstrateSet),ncol=1)
D835YAccessionNumbers<-matrix(data = Secondsubbackfreq[1,],ncol = 1)
#vectors to house the second set of substrates and accession numbers

for (i in 1:nrow(SecondSubstrateSet)){
  D835letters<-SecondSubstrateSet[i,4:18]
  D835letters<-D835letters[D835letters !="XXXXX"]
  D835letters<-paste(D835letters, sep="", collapse="")
  leftspaces<-c()
  rightspaces<-c()
  
  YYYmotif <- unlist(strsplit(D835letters, split = ""))
  YYYposition <- match(x = "x", table = YYYmotif)
  #position itself tells me how much is to the left of that X by what it's number is.  x at position 4 tells me that there are
  #just 3 letters to the left of x
  
  YYYLettersToTheLeft <- YYYposition - 1
  #how many letters to the right SHOULD just be length(motif)-position-1 if it's 5 long and x is at 3 then Y is at 4 and there is
  #just 1 spot to the right of Y so LettersToTheRight<-1 because 5-3-1=1
  YYYLettersToTheRight <- length(YYYmotif) - YYYposition - 1
  #then sanity check, we're currently looking only at +/-4, but this spot allows for up to +/- 7 as well, just depends on what the
  #variable the user puts in is
  if (YYYLettersToTheLeft < 7 | YYYLettersToTheRight < 7) {
    leftspaces<-rep(" ",times=(7-YYYLettersToTheLeft))
    rightspaces<-rep(" ",times=7-(YYYLettersToTheRight))
    #add blank spaces if the motif has less than 4 letters to the left/right
    motif<-c(leftspaces,YYYmotif,rightspaces)
    #save that motif, which is the Y and +/- 4 amino acids, including truncation
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    D835letters<-motif
    D835Ymotifs[i,1]<-D835letters
  }
  
  if(YYYLettersToTheLeft>6 && YYYLettersToTheRight>6){
    motif<-YYYmotif
    #add blank spaces if the motif has less than 4 letters to the left/right
    motif<-c(leftspaces,YYYmotif,rightspaces)
    #save that motif, which is the Y and +/- 4 amino acids, including truncation
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    D835letters<-motif
    D835Ymotifs[i,1]<-D835letters
  }
}


ITDmotifs=matrix(,nrow = nrow(ThirdSubstrateSet),ncol=1)
ITDAccessionNumbers<-matrix(data = Thirdsubbackfreq[1,],ncol = 1)
#third set of substrates and accession numbers


for (i in 1:nrow(ThirdSubstrateSet)){
  ITDletters<-ThirdSubstrateSet[i,4:18]
  ITDletters<-ITDletters[ITDletters !="XXXXX"]
  ITDletters<-paste(ITDletters, sep="", collapse="")
  YYYmotif <- unlist(strsplit(ITDletters, split = ""))
  leftspaces<-c()
  rightspaces<-c()
  YYYposition <- match(x = "x", table = YYYmotif)
  #position itself tells me how much is to the left of that X by what it's number is.  x at position 4 tells me that there are
  #just 3 letters to the left of x
  
  YYYLettersToTheLeft <- YYYposition - 1
  #how many letters to the right SHOULD just be length(motif)-position-1 if it's 5 long and x is at 3 then Y is at 4 and there is
  #just 1 spot to the right of Y so LettersToTheRight<-1 because 5-3-1=1
  YYYLettersToTheRight <- length(YYYmotif) - YYYposition - 1
  #then sanity check, we're currently looking only at +/-4, but this spot allows for up to +/- 7 as well, just depends on what the
  #variable the user puts in is
  if (YYYLettersToTheLeft < 7 | YYYLettersToTheRight < 7) {
    leftspaces<-rep(" ",times=(7-YYYLettersToTheLeft))
    rightspaces<-rep(" ",times=7-(YYYLettersToTheRight))
    #add blank spaces if the motif has less than 4 letters to the left/right
    motif<-c(leftspaces,YYYmotif,rightspaces)
    #save that motif, which is the Y and +/- 4 amino acids, including truncation
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    ITDletters<-motif
    ITDmotifs[i,1]<-ITDletters
    # ITDAccessionNumbers[i,1]<-FirstSubstrateSet[i,3]
  }
  
  if(YYYLettersToTheLeft>6 && YYYLettersToTheRight>6){
    motif<-YYYmotif
    #add blank spaces if the motif has less than 4 letters to the left/right
    motif<-c(leftspaces,YYYmotif,rightspaces)
    #save that motif, which is the Y and +/- 4 amino acids, including truncation
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    ITDletters<-motif
    ITDmotifs[i,1]<-ITDletters
    # ITDAccessionNumbers[i,1]<-FirstSubstrateSet[i,3]
  }
}


SubstrateOverlap1<-intersect(D835Ymotifs,ITDmotifs)
SubstrateOverlap1<-as.matrix(SubstrateOverlap1)
SubstrateOverlapFINAL<-intersect(FTLwtmotifs,SubstrateOverlap1)
#this tool is essentially looking for the intersection of three sets, and to do so it performs 1 intersection at a time
#this find the intersection of the substrate sets

AccessionOverlap1<-intersect(D835YAccessionNumbers,ITDAccessionNumbers)
AccessionOverlapFinal<-intersect(AccessionOverlap1,FTLwtAccessionNumbers)
AccessionOverlapFinal<-unlist(AccessionOverlapFinal)
#this tool is essentially looking for the intersection of three sets, and to do so it performs 1 intersection at a time
#this find the intersection of the accession number sets


for (x in 1:length(AccessionOverlapFinal)) {
  for (y in 1:ncol(Firstsubbackfreq)) {
    Acc<-AccessionOverlapFinal[x]
    SBF<-Firstsubbackfreq[1,y]
    if(Acc==SBF){
      FinalMatrix<-cbind(FinalMatrix,Firstsubbackfreq[,y])
    }
  }
}
#for every accession number, go back to the original substrate background frequency file and find the protein statistics associated with that number
#columnbind all those protein statistics together
FinalMatrix<-FinalMatrix[,2:ncol(FinalMatrix)]


#write all files, in the proper format so that the files coming out look like the files coming in

if(grepl(pattern = "Properties", x=FinalMatrix[22,1])==FALSE){
  Outputmatrix<-cbind(Firstsubbackfreq[,1],FinalMatrix)
  write.table(x=Outputmatrix,file = Shared_subbackfreq_table,quote = FALSE,sep = ",",row.names = FALSE,col.names = FALSE,na="")
} else {
  write.table(x=FinalMatrix,file = Shared_subbackfreq_table,quote = FALSE,sep = ",",row.names = FALSE,col.names = FALSE,na="")
}

SubstrateMatrix<-SubstrateHeader
if(ncol(SubstrateMatrix)>18){
  SubstrateMatrix<-SubstrateMatrix[,1:18]
}

for (z in 1:length(SubstrateOverlapFINAL)) {
  motif<-SubstrateOverlapFINAL[z]
  newmotif<-unlist(strsplit(motif,split = ""))
  
  Addition<-""
  outputmotif<-c(Addition,Addition,Addition,newmotif)
  SubstrateMatrix<-rbind(SubstrateMatrix,outputmotif)
}
write.table(x=SubstrateMatrix,file = Shared_motifs_table,quote = FALSE,sep = ",",row.names = FALSE,col.names = FALSE,na="")

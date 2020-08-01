#this brings in the input files
FirstSubstrateSet<- read.csv("S1.csv", stringsAsFactors=FALSE,colClasses = "character")
Firstsubbackfreq<- read.csv("SBF1.csv", header=FALSE, stringsAsFactors=FALSE)
SecondSubstrateSet<- read.csv("S2.csv", stringsAsFactors=FALSE,colClasses = "character")
Secondsubbackfreq<- read.csv("SBF2.csv", header=FALSE, stringsAsFactors=FALSE)


#this names the output files
First_unshared_motifs_table<-"1RS.csv"
First_unshared_subbackfreq<-"1RSBF.csv"
Second_unshared_motifs_table<-"2RS.csv"
Second_unshared_subbackfreq<-"2RSBF.csv"

#this creates the headers which comes from the input files, so that the output files can be given this header so that they will look identical to the input files
EmptySubHeader<-colnames(FirstSubstrateSet)
EmptySubHeader<-matrix(EmptySubHeader, nrow=1)
EmptySBFHeader<-Firstsubbackfreq[,1]


#ensure that all phospho-amino acids get marked with an "x" to denote their phosphoness
FirstCentralLetters<-FirstSubstrateSet[,11]
SecondCentralLetters<-SecondSubstrateSet[,11]

#use an 3 apply functions to create vectors, these vetors have true values where they find an S, T or Y.  
#so FirstEsses has a True anywhere it sees an S, and FirstTees has a True anywhere it sees a Y
FirstEsses<-sapply(FirstCentralLetters, grepl, pattern="S", ignore.case=TRUE)
FirstTees<-sapply(FirstCentralLetters, grepl, pattern="T", ignore.case=TRUE)
FirstWys<-sapply(FirstCentralLetters, grepl, pattern="Y", ignore.case=TRUE)

#do the same for the second substrate set's central letters
SecondEsses<-sapply(SecondCentralLetters, grepl, pattern="S", ignore.case=TRUE)
SecondTees<-sapply(SecondCentralLetters, grepl, pattern="T", ignore.case=TRUE)
SecondWys<-sapply(SecondCentralLetters, grepl, pattern="Y", ignore.case=TRUE)

#where there is a True value in FirstEsses, replace that value in the original substrate set with an xS.  This is because there was originally an S in that
#position, and the programs wants that S to be marked with an x, denoting phospho
FirstCentralLetters<-replace(FirstCentralLetters,FirstEsses,"xS")
FirstCentralLetters<-replace(FirstCentralLetters,FirstTees,"xT")
FirstCentralLetters<-replace(FirstCentralLetters,FirstWys,"xY")

SecondCentralLetters<-replace(SecondCentralLetters,SecondEsses,"xS")
SecondCentralLetters<-replace(SecondCentralLetters,SecondTees,"xT")
SecondCentralLetters<-replace(SecondCentralLetters,SecondWys,"xY")

#then put these x-marked letters back where they were found, in position 11 of the substrate sets
FirstCentralLetters->FirstSubstrateSet[,11]
SecondCentralLetters->SecondSubstrateSet[,11]

FTLwtmotifs=matrix(,nrow = nrow(FirstSubstrateSet),ncol=1)
FTLwtAccessionNumbers=matrix(,nrow = nrow(FirstSubstrateSet),ncol=1)

for (i in 1:nrow(FirstSubstrateSet)){
  FTLwtletters<-FirstSubstrateSet[i,4:18]
  FTLwtletters<-FTLwtletters[FTLwtletters !="XXXXX"]
  FTLwtletters<-paste(FTLwtletters, sep="", collapse="")
  leftspaces<-c()
  rightspaces<-c()
  
  #position itself tells one how much is to the left of that X by what it's number is.  x at position 4 tells me that there are
  #just 3 letters to the left of x
  YYYmotif <- unlist(strsplit(FTLwtletters, split = ""))
  YYYposition <- match(x = "x", table = YYYmotif)
  
  #how many letters to the right SHOULD just be length(motif)-position-1 if it's 5 long and x is at 3 then Y is at 4 and there is
  #just 1 spot to the right of Y so LettersToTheRight<-1 because 5-3-1=1
  YYYLettersToTheLeft <- YYYposition - 1
  YYYLettersToTheRight <- length(YYYmotif) - YYYposition - 1
  
  
  if (YYYLettersToTheLeft < 7 | YYYLettersToTheRight < 7) {
    #add blank spaces if the motif has less than 4 letters to the left/right
    leftspaces<-rep(" ",times=(7-YYYLettersToTheLeft))
    rightspaces<-rep(" ",times=7-(YYYLettersToTheRight))
    motif<-c(leftspaces,YYYmotif,rightspaces)
    #save that motif, which is the Y and +/- 4 amino acids, including truncation
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    FTLwtletters<-motif
    FTLwtmotifs[i,1]<-FTLwtletters
    FTLwtAccessionNumbers[i,1]<-FirstSubstrateSet[i,3]
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
    FTLwtAccessionNumbers[i,1]<-FirstSubstrateSet[i,3]
    
    
  }
  
}

D835Ymotifs=matrix(,nrow = nrow(SecondSubstrateSet),ncol=1)
D835YAccessionNumbers<-matrix(,nrow = nrow(SecondSubstrateSet),ncol = 1)

for (i in 1:nrow(SecondSubstrateSet)){
  D835letters<-SecondSubstrateSet[i,4:18]
  D835letters<-D835letters[D835letters !="XXXXX"]
  D835letters<-paste(D835letters, sep="", collapse="")
  leftspaces<-c()
  rightspaces<-c()
  
  YYYmotif <- unlist(strsplit(D835letters, split = ""))
  YYYposition <- match(x = "x", table = YYYmotif)
  
  YYYLettersToTheLeft <- YYYposition - 1
  YYYLettersToTheRight <- length(YYYmotif) - YYYposition - 1
  if (YYYLettersToTheLeft < 7 | YYYLettersToTheRight < 7) {
    leftspaces<-rep(" ",times=(7-YYYLettersToTheLeft))
    rightspaces<-rep(" ",times=7-(YYYLettersToTheRight))
    motif<-c(leftspaces,YYYmotif,rightspaces)
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    D835letters<-motif
    D835Ymotifs[i,1]<-D835letters
    D835YAccessionNumbers[i,1]<-SecondSubstrateSet[i,3]
  }
  
  if(YYYLettersToTheLeft>6 && YYYLettersToTheRight>6){
    motif<-YYYmotif
    motif<-c(leftspaces,YYYmotif,rightspaces)
    motif<-motif[!motif %in% "x"]
    motif<-paste(motif, sep="", collapse="")
    D835letters<-motif
    D835Ymotifs[i,1]<-D835letters
    D835YAccessionNumbers[i,1]<-SecondSubstrateSet[i,3]
  }
}

names(FTLwtmotifs)<-FTLwtAccessionNumbers
names(D835Ymotifs)<-D835YAccessionNumbers


FTLwtmotifsFINAL<-FTLwtmotifs[!FTLwtmotifs %in% D835Ymotifs]
FTLwtmotifsFINAL<-FTLwtmotifsFINAL[!duplicated(FTLwtmotifsFINAL)]

D835YmotifsFINAL<-D835Ymotifs[!D835Ymotifs %in% FTLwtmotifs]
D835YmotifsFINAL<-D835YmotifsFINAL[!duplicated(D835YmotifsFINAL)]


columnalheader<-c(rep(NA,36))
FTLFinalMatrix<-matrix(data =columnalheader,nrow = 1)

FLTheader<-c("Substrate","Species","Reference","-7","-6","-5","-4","-3","-2","-1","0","1","2","3","4","5","6","7","Phosphite")

if (length(FTLwtmotifsFINAL)>0){
  for (k in 1:length(FTLwtmotifsFINAL)) {
    AN<-00000
    #it is necessary to destroy the accession number multiple times to ensure it is
    #destroyed immediately after use
    for (m in 1:ncol(Firstsubbackfreq)) {
      AN <- as.character(Firstsubbackfreq[1, m])
      if (grepl(pattern = AN,
                x = names(FTLwtmotifsFINAL[k]),
                fixed = TRUE) == TRUE) {
        outputmatrix <- as.character(Firstsubbackfreq[, m])
        outputmatrix <- matrix(outputmatrix, nrow = 1)
        #with that accession number, find a match in the subbackfreq file and save it here
        FTLFinalMatrix<-rbind(FTLFinalMatrix,outputmatrix)
      }
    }
  }
  FTLFinalMatrix<-FTLFinalMatrix[!duplicated(FTLFinalMatrix),]
  FTLFinalMatrix<-FTLFinalMatrix[2:nrow(FTLFinalMatrix),]
  
  
  FTLoutputmatrix<-matrix(data=c(FTLwtmotifsFINAL,names(FTLwtmotifsFINAL)),ncol = 2)
  # FLTheader<-unlist(FLTheader)
  lefthandFLT<-matrix(data = rep(NA,times=2*nrow(FTLoutputmatrix)),nrow=nrow(FTLoutputmatrix))
  righthandFLT<-matrix(data = rep(NA,times=1*nrow(FTLoutputmatrix)),nrow=nrow(FTLoutputmatrix))
  FLTaccessionset<-FTLoutputmatrix[,2]
  FTLmeat<-sapply(FTLoutputmatrix[,1], strsplit, "")
  FTLmeat<-sapply(FTLmeat, unlist)
  colnames(FTLmeat)<-NULL
  FTLmeat<-t(FTLmeat)
  
  FTLoutputmatrix2<-cbind(lefthandFLT,FLTaccessionset,FTLmeat,righthandFLT)
  colnames(FTLoutputmatrix2)<-NULL
  rownames(FTLoutputmatrix2)<-NULL
  colnames(FLTheader)<-NULL
  rownames(FLTheader)<-NULL
  
  
  FirstCentralLettersAGAIN<-FTLoutputmatrix2[,11]
  
  FirstEsses<-sapply(FirstCentralLettersAGAIN, grepl, pattern="S", ignore.case=TRUE)
  FirstTees<-sapply(FirstCentralLettersAGAIN, grepl, pattern="T", ignore.case=TRUE)
  FirstWys<-sapply(FirstCentralLettersAGAIN, grepl, pattern="Y", ignore.case=TRUE)
  
  FirstCentralLettersAGAIN<-replace(FirstCentralLettersAGAIN,FirstEsses,"xS")
  FirstCentralLettersAGAIN<-replace(FirstCentralLettersAGAIN,FirstTees,"xT")
  FirstCentralLettersAGAIN<-replace(FirstCentralLettersAGAIN,FirstWys,"xY")
  
  FirstCentralLettersAGAIN->FTLoutputmatrix2[,11]
  
  FTLoutputmatrix2<-rbind(FLTheader,FTLoutputmatrix2)
  
  write.table(x=FTLoutputmatrix2,
              file=First_unshared_motifs_table,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
  
  columnalheader<-c(as.character(Firstsubbackfreq[1:36,1]))
  columnalheader<-matrix(columnalheader,nrow = 1)
  write.table(x=columnalheader,
              file=First_unshared_subbackfreq,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
  
  write.table(x=FTLFinalMatrix,
              file=First_unshared_subbackfreq,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
} else{
  FTLFinalMatrix<-columnalheader
  write.table(x=EmptySubHeader,
              file=First_unshared_motifs_table,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
  
  columnalheader<-c(as.character(Firstsubbackfreq[1:36,1]))
  columnalheader<-matrix(columnalheader,nrow = 1)
  write.table(x=columnalheader,
              file=First_unshared_subbackfreq,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
}


columnalheader<-c(rep(NA,36))
D835YFinalMatrix<-matrix(data =columnalheader,nrow = 1)

if (length(D835YmotifsFINAL)>0){
  for (k in 1:length(D835YmotifsFINAL)) {
    #it is necessary to destroy the accession number multiple times to ensure it is
    #destroyed immediately after use
    for (m in 1:ncol(Secondsubbackfreq)) {
      AN <- as.character(Secondsubbackfreq[1, m])
      if (grepl(pattern = AN,
                x = names(D835YmotifsFINAL[k]),
                fixed = TRUE) == TRUE) {
        outputmatrix <- as.character(Secondsubbackfreq[, m])
        outputmatrix <- matrix(outputmatrix, nrow = 1)
        #with that accession number, find a match in the subbackfreq file and save it here
        D835YFinalMatrix<-rbind(D835YFinalMatrix,outputmatrix)
      }
    }
  }
  D835YFinalMatrix<-D835YFinalMatrix[!duplicated(D835YFinalMatrix),]
  D835YFinalMatrix<-D835YFinalMatrix[2:nrow(D835YFinalMatrix),]
  
  D835Youtputmatrix<-matrix(data=c(D835YmotifsFINAL,names(D835YmotifsFINAL)),ncol = 2)
  
  D835Yheader<-c("Substrate","Species","Reference","-7","-6","-5","-4","-3","-2","-1","0","1","2","3","4","5","6","7","Phosphite")
  # D835Yheader<-unlist(D835Yheader)
  lefthandD835<-matrix(data = rep(NA,times=2*nrow(D835Youtputmatrix)),nrow=nrow(D835Youtputmatrix))
  righthandD835<-matrix(data = rep(NA,times=1*nrow(D835Youtputmatrix)),nrow=nrow(D835Youtputmatrix))
  D835Yaset<-D835Youtputmatrix[,2]
  D835meat<-sapply(D835Youtputmatrix[,1], strsplit, "")
  D835meat<-sapply(D835meat, unlist)
  colnames(D835meat)<-NULL
  D835meat<-t(D835meat)
  
  D835Youtputmatrix2<-cbind(lefthandD835,D835Yaset,D835meat,righthandD835)
  colnames(D835Youtputmatrix2)<-NULL
  rownames(D835Youtputmatrix2)<-NULL
  colnames(D835Yheader)<-NULL
  rownames(D835Yheader)<-NULL
  
  
  SecondCentralLettersAGAIN<-D835Youtputmatrix2[,11]
  
  SecondEsses<-sapply(SecondCentralLettersAGAIN, grepl, pattern="S", ignore.case=TRUE)
  SecondTees<-sapply(SecondCentralLettersAGAIN, grepl, pattern="T", ignore.case=TRUE)
  SecondWys<-sapply(SecondCentralLettersAGAIN, grepl, pattern="Y", ignore.case=TRUE)
  
  SecondCentralLettersAGAIN<-replace(SecondCentralLettersAGAIN,SecondEsses,"xS")
  SecondCentralLettersAGAIN<-replace(SecondCentralLettersAGAIN,SecondTees,"xT")
  SecondCentralLettersAGAIN<-replace(SecondCentralLettersAGAIN,SecondWys,"xY")
  
  SecondCentralLettersAGAIN->D835Youtputmatrix2[,11]
  
  D835Youtputmatrix2<-rbind(D835Yheader,D835Youtputmatrix2)
  
  write.table(x=D835Youtputmatrix2,
              file=Second_unshared_motifs_table,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
  
  columnalheader<-c(as.character(Firstsubbackfreq[1:36,1]))
  columnalheader<-matrix(columnalheader,nrow = 1)
  write.table(x=columnalheader,
              file=Second_unshared_subbackfreq,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
  
  write.table(x=D835YFinalMatrix,
              file=Second_unshared_subbackfreq,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
} else {
  D835YFinalMatrix<- columnalheader
  write.table(x=EmptySubHeader,
              file=Second_unshared_motifs_table,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
  
  columnalheader<-c(as.character(Firstsubbackfreq[1:36,1]))
  columnalheader<-matrix(columnalheader,nrow = 1)
  write.table(x=columnalheader,
              file=Second_unshared_subbackfreq,
              quote=FALSE, sep=",",
              row.names=FALSE,col.names = FALSE, na="", append=TRUE)
}

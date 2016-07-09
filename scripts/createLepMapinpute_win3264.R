createLepMapinpute <- function(pedfaminfo,genoplinkped,outname){
  cat('\n')
  pedi <- read.table(paste(pedfaminfo,sep=''),header=F,stringsAsFactors=F)
  cat('.... FAM file with pedigree information imported .... \n')
  
  colnames(pedi) <- c('FIID','IID','Sire','Dam','Sex','Pheno')
  pedi <-  pedi[order(pedi$Sire,pedi$Dam),]
  
  hmsires <- nrow(pedi[which(pedi$Sire==0 & pedi$Dam!=0),])
  hmdams <- nrow(pedi[which(pedi$Sire!=0 & pedi$Dam==0),])
  if(hmsires>0){DumSires <- paste('DumS',1:hmsires,sep='')}
  if(hmdams>0){DumDams <- paste('DumD',1:hmdams,sep='')}
  
  if(hmsires>0 & hmdams>0){
    dumped <- data.frame(IID=c(DumSires,DumDams))
    dumped <- cbind.data.frame(dumped,dumped,0,0,0,0)
    colnames(dumped) <- c('FIID','IID','Sire','Dam','Sex','Pheno')
    pedi <- rbind.data.frame(dumped,pedi)
    pedi[which(pedi$Sire==0 & pedi$Dam!=0),'Sire'] <- DumSires
    pedi[which(pedi$Sire!=0 & pedi$Dam==0),'Dam'] <- DumDams
  } else if(hmsires>0 & hmdams<1){
    dumped <- data.frame(IID=DumSires)
    dumped <- cbind.data.frame(dumped,dumped,0,0,0,0)
    colnames(dumped) <- c('FIID','IID','Sire','Dam','Sex','Pheno')
    pedi <- rbind.data.frame(dumped,pedi)
    pedi[which(pedi$Sire==0 & pedi$Dam!=0),'Sire'] <- DumSires
  } else if(hmsires<1 & hmdams>0){
    dumped <- data.frame(IID=DumDams)
    dumped <- cbind.data.frame(dumped,dumped,0,0,0,0)
    colnames(dumped) <- c('FIID','IID','Sire','Dam','Sex','Pheno')
    pedi <- rbind.data.frame(dumped,pedi)
    pedi[which(pedi$Sire!=0 & pedi$Dam==0),'Dam'] <- DumDams
  }
  
  offsprings <- pedi[which(pedi$Sire!=0 & pedi$Dam!=0),]
  sires <- sort(unique(offsprings$Sire)); dams <- sort(unique(offsprings$Dam))
  
  offsprings$famID <-  paste(offsprings$Sire,'><',offsprings$Dam,sep='')
  famILY <- data.frame(FSFAM=sort(table(offsprings$famID)))
  famILY$Sire <- paste('S000',1:nrow(famILY),sep='')
  famILY$Dam <- paste('D000',1:nrow(famILY),sep='')
  famILY$famID <- rownames(famILY)
  famILY$famIDnr <- 1:nrow(famILY)
  
  offsprings.ped <- matrix(0,nrow=nrow(offsprings),ncol=8)
  colnames(offsprings.ped) <- c('FIID','IID','Sire','Dam','Sex','oldSire','oldDAM','oldIID')
  for(i in 1:nrow(offsprings)){
    offsprings.ped[i,'IID'] <- as.vector(offsprings[i,'IID'])
    offsprings.ped[i,'Sex'] <- offsprings[i,'Sex']
    offsprings.ped[i,'oldSire'] <- offsprings[i,'Sire']
    offsprings.ped[i,'oldDAM'] <- offsprings[i,'Dam']
    offsprings.ped[i,'oldIID'] <- as.vector(offsprings[i,'IID'])
    getpar <- merge(offsprings[i,],famILY,by='famID')
    offsprings.ped[i,c('FIID')] <- getpar[,c('famIDnr')]
    offsprings.ped[i,c('Sire')] <- getpar[,c('Sire.y')]
    offsprings.ped[i,c('Dam')] <- getpar[,c('Dam.y')]
  }
  
  parent.ped <- matrix(0,nrow=nrow(famILY)*2,ncol=6)
  colnames(parent.ped) <- c('FIID','IID','Sire','Dam','Sex','oldIID')
  parent <- data.frame(IID=c(famILY$Sire,famILY$Dam))
  parent$Sex <- c(rep('1',length(famILY$Sire)),rep('2',length(famILY$Dam)))
  parent$famIDnr <- c(famILY$famIDnr,famILY$famIDnr)
  parent$famID <- c(famILY$famID,famILY$famID)
  
  for(j in 1:nrow(parent)){
    parent.ped[j,'IID'] <- as.vector(parent[j,'IID'])
    parent.ped[j,'Sex'] <- parent[j,'Sex']
    parent.ped[j,'FIID'] <- parent[j,'famIDnr']
    
    if(parent[j,'Sex']=='1'){
      getoldid <- unlist(strsplit(parent[j,'famID'],split='><'))[1]
    } else if(parent[j,'Sex']=='2'){
      getoldid <- unlist(strsplit(parent[j,'famID'],split='><'))[2]
    } 
    parent.ped[j,c('oldIID')] <- as.vector(getoldid)
  }
  
  parentoffsring <- rbind.data.frame(parent.ped,offsprings.ped[,c(1:5,8)])
  parentoffsring <- parentoffsring[order(as.numeric(as.vector(parentoffsring$FIID))),]
  
  cat('.... Data editing on FAM (pedigree information) done .... \n')
  
  cat('.... importing genotype data (PLINK ped file format) .... \n')
  genodata <- read.table(paste(genoplinkped),stringsAsFactors = F)
  cat('.... Genotypes (PLINK ped file format) imported .... \n')
  genoiid <- data.frame(IID=genodata[,1]); genoiid$nr <- 1:nrow(genoiid)
  
  cat('\n')
  cat('.... Preparing LepMAP file format .... \n')
  cat('.... Note : This might take a while .... \n')
  cat('\n')
  
  genofile <- data.frame(matrix(0,nrow(parentoffsring),ncol(genodata)))
  iterchecks.anim <- round(nrow(genofile)/10,digits=0)
  for(m in 1:nrow(parentoffsring)){
    anim <- data.frame(IID=as.vector(parentoffsring[m,c('oldIID')]))
    anim.dat <- merge(anim,genoiid,by=1)
    if(nrow(anim.dat)==0){
      genofile[m,1] <- as.vector(parentoffsring[m,1])
      genofile[m,2] <- as.vector(parentoffsring[m,2])
      genofile[m,3] <- as.vector(parentoffsring[m,3])
      genofile[m,4] <- as.vector(parentoffsring[m,4])
      genofile[m,5] <- as.vector(parentoffsring[m,5])
      
    } else {
      genofile[m,1] <- as.vector(parentoffsring[m,1])
      genofile[m,2] <- as.vector(parentoffsring[m,2])
      genofile[m,3] <- as.vector(parentoffsring[m,3])
      genofile[m,4] <- as.vector(parentoffsring[m,4])
      genofile[m,5] <- as.vector(parentoffsring[m,5])
      genofile[m,7:ncol(genofile)] <- genodata[anim.dat$nr,7:ncol(genofile)]
    }
    if(m %% iterchecks.anim==0){
      cat(paste('-- now at animal ... ',m,' ... out of ... ',nrow(parentoffsring),sep=''),' \n')
    }
  }
  cat('\n')
  cat('.... LepMAP file format done and writing data to disk .... \n')
  
  write.table(genofile,paste(outname,'.linkage',sep=''),quote=F,col.names=F,row.names=F)
  write.table(parentoffsring,paste(outname,'.oldids',sep=''),quote=F,col.names=T,row.names=F)
  
  cat('.... Completed enjoy LepMAPing .... \n')
  cat('\n\n')
  
  cat('@--------------------------------------------------------@\n')
  cat('@                                                        @\n')
  cat('@                Created by S.A. Boison ::: May 10, 2016 @\n')
  cat('@               please report bugs to soloboan@yahoo.com @\n')
  cat('@--------------------------------------------------------@\n')
}


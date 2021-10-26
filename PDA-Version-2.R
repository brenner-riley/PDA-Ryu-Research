####Read data
#This data is from Supplementary Table 3 (HCD data) from Marx et al. 2013
data <- read.csv(file='/Users/rileybrenner/Desktop/RStudio Files/41587_2013_BFnbt2585_MOESM14_ESM.csv', sep=',')
head(data)
names(data)

setwd("/Users/rileybrenner/Desktop/HCD Data")
files <- lapply(list.files(path = "/Users/rileybrenner/Desktop/HCD Data",pattern = ".csv"), read.csv)

oxidation_data <- data.frame()
for(i in 1:length(files))
{
  oxidation_data <- rbind(data.frame("Sequence" = files[[i]]$Sequence, "First.Scan" = files[[i]]$First.Scan, "Spectrum.File" = files[[i]]$Spectrum.File),oxidation_data)
}



####Note
#Problems with data: Multiple different peptides were found in the same scan, and we cannot differentiate them due to technical difficulties.
#Thus, fitting a linear regression model or doing hypothesis tests, observations are not independent. 
#Visualization will be OK with this problem. 
#Later, for a linear regression model and hypothesis tests, we will use a different dataset for a training set. (Ok to use this for the test set.) 
#For now, let students practice coding. 

####Keep peptides with E.value < 0.01, Classifications..FDR.=TRUE, and (Classifications..FLR.=TRUE or unmodified peptide)
f.index=which((data$E.Value<0.01) & (data$Classifications..FDR.==TRUE) & ((data$Classifications..FLR.==TRUE) | (data$Phosphosite.Position=="")))
f.data <- data[f.index,]

boxplot(f.data$Retention.Time..min. ~ nchar(f.data$Peptide.Sequence),  horizontal = TRUE, xlab = "Retention Time", ylab = "Peptide Length", main = " ", col = "#00A5E0")
hist(f.data$Retention.Time..min. , xlab = "Retention Time", main = " ", col = "#F2AF29")

f.data=select(f.data , c("Spectrum","Retention.Time..min.","Intensity..arb..Units.","Peptide.Sequence","Phosphosite.Position"))
f.data$Oxidation.Location<-c(NA)

####lowercase for phosphosite


table(f.data$Phosphosite.Position) #note that all has only single-phosphorylation in this dataset
phospho <- function(x){ 
  seq <- x[1]; phospho.loc <- x[2]
  phospho.loc <- as.numeric(phospho.loc) #now, empty phospho.loc is NA
  if (!is.na(phospho.loc)){
    temp = strsplit(seq, "")[[1]] #another way instead of unlist()
    temp[phospho.loc] = tolower(temp[phospho.loc]) #change to lowercase for modified amino acid
    seq=paste(temp, collapse="")
  } 
  return(seq)
}

####lowercase for oxidation
unoxidized <- function(x){
  seq = x[2]
  temp <- c()
  if(str_detect(seq, "s")){
    temp = strsplit(seq, "s") %>%
      unlist()
    seq = paste(toupper(temp[1]),"s", toupper(temp[2]), sep="")
  }else if(str_detect(seq, "t")){
    temp = strsplit(seq, "t") %>%
      unlist()
    seq = paste(toupper(temp[1]),"t", toupper(temp[2]), sep="")
  }else if(str_detect(seq, "y")){
    temp = strsplit(seq, "y") %>%
      unlist()
    seq = paste(toupper(temp[1]),"y", toupper(temp[2]), sep="")
  }else{
    seq=toupper(seq)
  }
  return(seq)
}
oxidation <- function(x){ 
  seq <- x[1]; oxi.loc <- x[2]
  oxi.loc <- as.numeric(oxi.loc) #now, empty phospho.loc is NA
  if (!is.na(oxi.loc)){
    temp = strsplit(seq, "")[[1]] #another way instead of unlist()
    temp[oxi.loc] = tolower(temp[oxi.loc]) #change to lowercase for modified amino acid
    seq=paste(temp, collapse="")
  } 
  return(seq)
}

{
  spec_id <- f.data$Spectrum %>%
    strsplit(" ") %>%
    unlist()
  f.data$Spectrum = spec_id[seq(3, length(spec_id), 3)]
  oxidation_subset = oxidation_data[which(str_detect(oxidation_data$Sequence, "m")),c(2,1,3)]
  oxidation_subset <- oxidation_subset[which(oxidation_subset$First.Scan %in% f.data$Spectrum),]
  
  seq <- data.frame("Spectrum" = oxidation_subset$First.Scan,"Retention.Time..min."= NA, 
                    "Intensity..arb..Units." = NA,"Peptide.Sequence" = toupper(oxidation_subset$Sequence),
                    "Phosphite.Position" = NA,"Oxidation.Location" = sapply(oxidation_subset$Sequence , regexpr , pattern = "m"),
                    "Peptide.Sequence2" = apply(oxidation_subset, 1, unoxidized ))
  f.data$Peptide.Sequence2 <- apply(cbind(f.data$Peptide.Sequence,f.data$Phosphosite.Position), 1, phospho)
  names(seq) <- names(f.data)
  f.data <- rbind(f.data, seq)
  
  f.data.oxi <- f.data %>%
    group_by(Spectrum, Peptide.Sequence2) %>%
    summarise(Peptide.Sequence,Retention.Time..min., max(Oxidation.Location, na.rm = T), Intensity..arb..Units.)
  
  f.data <- f.data.oxi[-which(is.na(f.data.oxi$Retention.Time..min.)),]
  f.data[which(f.data$`max(Oxidation.Location, na.rm = T)` == -Inf),5]<-NA 
  
  
  f.data$Peptide.Sequence2 <- apply(cbind(f.data$Peptide.Sequence2,f.data$`max(Oxidation.Location, na.rm = T)`), 1, oxidation)
  f.data <- f.data[-5]
  
}

library_data <- oxidation_data[which(oxidation_data$Sequence %in% f.data$Peptide.Sequence2),]

raw_files <- library_data$Spectrum.File%>%
  strsplit("~") %>%
  unlist() 
raw_files <- raw_files[seq(7, length(raw_files), 7)] %>%
  strsplit(".", fixed = TRUE) %>%
  unlist()
library_data$Spectrum.File <- raw_files[seq(1, length(raw_files), 4)]

f.data$Library <- NA
library_data <- data.frame("Spectrum" = as.character(library_data$First.Scan),"Peptide.Sequence2" = library_data$Sequence, 
                           "Peptide.Sequence" =NA, "Retention.Time..min." = NA,"Intensity..arb.Units" = NA, 
                           "Library" = library_data$Spectrum.File)
f.data <- rbind(f.data, library_data)

f.data.library <- f.data %>%
  group_by(Spectrum, Peptide.Sequence2) %>%
  summarise(Peptide.Sequence,Retention.Time..min., Intensity..arb..Units., max(as.numeric(Library), na.rm = T))

f.data <- f.data.library[-which(is.na(f.data.library$Retention.Time..min.)),]
f.data <- f.data[-which(f.data$`max(as.numeric(Library), na.rm = T)`== -Inf),]

colnames(f.data) <- colnames(library_data)

seq <- c()
rt <- c()
set.seed(3)
for(x in 1:96)
{
  t_1 <- NA
  t_2 <- NA
  group <- f.data[which(f.data$Peptide.Sequence2 == f.data$Peptide.Sequence),]
  group <- group[which(group$Library == x),]
  if(nrow(group)>0)
  {
    loc <- sample(1:nrow(group),1)
    t_1 <- group$Peptide.Sequence[loc]
    t_2 <- group$Retention.Time..min.[loc]
  }
  seq <- append(seq, t_1, length(seq))
  rt <- append(rt, t_2, length(rt))
}

ref_peps <- data.frame("Sequence" = seq, "RT..min." = rt)

##Finds abundance of each Amino Acid in a peptide sequence.
count_each <- function(x)
{
  A <- str_count(x,"A")
  R <- str_count(x,"R")
  N <- str_count(x,"N")
  D <- str_count(x,"D")
  C <- str_count(x,"C")
  Q <- str_count(x,"Q")
  E <- str_count(x,"E")
  G <- str_count(x,"G")
  H <- str_count(x,"H")
  I <- str_count(x,"I")
  L <- str_count(x,"L")
  K <- str_count(x,"K")
  M <- str_count(x,"M")
  F_ <- str_count(x,"F")
  P <- str_count(x,"P")
  S <- str_count(x,"S")
  T_ <- str_count(x,"T")
  W <- str_count(x,"W")
  Y <- str_count(x,"Y")
  V <- str_count(x,"V")
  s <- str_count(x,"s")
  t <- str_count(x,"t")
  y <- str_count(x,"y")
  m <- str_count(x,"m")
  return(c(A,R,N,D,C,Q,E,G,H,I,L,K,M,F_,P,S,T_,W,Y,V,s,t,y,m))
}
##finds the difference in Amino Acid abundance between two peptide sequence,
## and difference in retention time
dif_each <- function(x,y)
{
  seq_count <- x[c(1:24)]; loc<-x[25]
  seed_count<-y[loc,c(1:24)]
  dif_count <- as.numeric(seq_count) - seed_count
  if(!is.na(y[loc,25])){
    dif_count <- append(dif_count, as.numeric(x[26])-y[loc,25], length(dif_count))
    return(dif_count)
  } else {
    temp <- rep(NA, 25)
    return(temp)
  }
}


##Prepares count and ms_count for processing in dif_each,
##done by using count_each and binding some aditional columns from original files
{
  count <- t(sapply(ref_peps$Sequence, count_each))
  rownames(count)<-c(1:nrow(count))
  count<- cbind(count, ref_peps$RT..min.)
  
  ms_count <- t(sapply(f.data$Peptide.Sequence2, count_each))
  ms_count <- cbind(ms_count, f.data$Library)
  ms_count <-cbind(ms_count, f.data$Retention.Time..min.)
}

temp<- t(apply(ms_count, 1, dif_each, count))

colnames(temp)<-c("A","R","N","D","C","Q","E","G","H",
                  "I","L","K","M","F","P","S","T","W",
                  "Y","V","s","t","y","m","RT.change")
rownames(temp)<-c(1:nrow(temp))

##Linear Regression and values associated with them
lin_Reg <- lm(temp[,25]~temp[,c(1:24)])
summary(lin_Reg)
lin_Reg
confint(lin_Reg)




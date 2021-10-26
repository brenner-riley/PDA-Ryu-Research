library(stringr)
library(dplyr)

##Reading the original seed peptides from the study, if we use another set of data for verification then
##this will probably be completely unnecessary and we will just look for the first or a random unmodified 
##peptide to calculate the differential values.
seed_peps <- read.csv("/Users/rileybrenner/Desktop/RStudio Files/STable2_41587_2013_BFnbt2585_MOESM13_ESM.csv")

##How I initially got around not having all the seed peptides occur in the data files, again if we use another
##set of data this will be unnecessary.
createAltSeeds <- function(x)
{
  temp <- x
  loc <- str_locate(x, "p")+1
  substr(x, loc, loc) <- "S"
  alt_1 <- x
  substr(x, loc, loc) <- "T"
  alt_2 <- x
  substr(x, loc, loc) <- "Y"
  alt_3 <-x
  ret <- c(alt_1, alt_2, alt_3)
  ret <- ret[-which(ret == temp)]
  return(ret)
}


#Reading in all the files from a directory.
setwd("/Users/rileybrenner/Desktop/HCD Data")
files <- lapply(list.files(path = "/Users/rileybrenner/Desktop/HCD Data",pattern = ".csv"), read.csv)

ms_collected_peps <- data.frame()
for(i in 1:length(files))
{
  ms_collected_peps <- rbind(ms_collected_peps,data.frame(files[[i]]$Sequence,files[[i]]$RT..min.,files[[i]]$Spectrum.File))
}
colnames(ms_collected_peps)<-c("Sequence","RT..min.","Spectrum.File")


##Part of the process by which I found alternative reference peptides for the analysis.
alt_peps <- data.frame(t(sapply(seed_peps$Seed.Peptide.Sequence, createAltSeeds)))
seed_peps$Seed.Peptide.Sequence = sapply(seed_peps$Seed.Peptide.Sequence, str_remove, pattern = "p")
row.names(alt_peps)<- c(1:nrow(alt_peps))
colnames(alt_peps)<-c("FirstAlt", "SecondAlt")
alt_peps$FirstAlt = sapply(alt_peps$FirstAlt, str_remove, pattern = "p")
alt_peps$SecondAlt = sapply(alt_peps$SecondAlt, str_remove, pattern = "p")
##Looks for the alternate peptides in ms_collected to have a reference retention time.
for(i in 1:nrow(seed_peps))
{
  locs = which(ms_collected_peps$Sequence  == alt_peps$FirstAlt[i] )
  alt_peps$First.RT..min.[i] <- median(ms_collected_peps$RT..min.[locs])
  locs = which(ms_collected_peps$Sequence == alt_peps$SecondAlt[i] )
  alt_peps$Second.RT..min.[i] <- median(ms_collected_peps$RT..min.[locs])
}
alt_peps$Seed.Num <-c(1:96)
##Looks for seed peptides in ms_collected to grab their retention time.
for(i in 1:nrow(seed_peps))
{
  locs = which(ms_collected_peps$Sequence == seed_peps$Seed.Peptide.Sequence[i])
  seed_peps$RT..min.[i] <- median(ms_collected_peps$RT..min.[locs])
}


##Turns Spectrum.File into the associated reference peptide value
raw_files <- ms_collected_peps$Spectrum.File%>%
  strsplit("~") %>%
  unlist() 
raw_files <- raw_files[seq(7, length(raw_files), 7)] %>%
  strsplit(".", fixed = TRUE) %>%
  unlist()
ms_collected_peps$Spectrum.File <- raw_files[seq(1, length(raw_files), 4)]

##Following 60ish lines is the innefective way I collected the reference peptides in order to make the 
##call of dif_each easier.
{
  alt_peps <- alt_peps[which(is.na(seed_peps$RT..min.)),]
  temp <- alt_peps[which(is.na(alt_peps$First.RT..min.)&is.na(alt_peps$Second.RT..min.)),]
  alt_peps <- alt_peps[-which(is.na(alt_peps$First.RT..min.)&is.na(alt_peps$Second.RT..min.)),]
  seed_peps <- seed_peps[-which(is.na(seed_peps$RT..min.)),]
  alt_1 <-alt_peps[-which(is.na(alt_peps$First.RT..min.)),]
  alt_2<-alt_peps[which(is.na(alt_peps$First.RT..min.)),]

  seq <- c()
  rt <- c()
  seed <- c()


  seq <- append(seq, alt_1$FirstAlt, length(seq))
  rt <- append(rt, alt_1$First.RT..min., length(rt))
  seed <- append(seed, alt_1$Seed.Num, length(seed))

  seq <- append(seq, alt_2$SecondAlt, length(seq))
  rt <- append(rt, alt_2$Second.RT..min., length(rt))
  seed <- append(seed, alt_2$Seed.Num, length(seed))

  seq <- append(seq, seed_peps$Seed.Peptide.Sequence, length(seq))
  rt <- append(rt, seed_peps$RT..min., length(rt))
  seed <- append(seed, seed_peps$Library, length(seed))

  ##Looks for an unmodified reference peptide for the files without the seed, or alternative seed.
  for(x in temp$Seed.Num)
  {
    group <- which(ms_collected_peps$Spectrum.File == x)
    for(y in group)
    {
      if(toupper(ms_collected_peps$Sequence[y])==ms_collected_peps$Sequence[y])
      {
        seq <- append(seq, ms_collected_peps$Sequence[y], length(seq))
        rt <- append(rt, ms_collected_peps$RT..min.[y], length(rt))
        seed <- append(seed, ms_collected_peps$Spectrum.File[y], length(seed))
        break
      }
    }
  }



  ref_peps <-data.frame("Seed.Peptide.Sequence" = seq, "RT..min." = rt, "Seed" = seed)
  rownames(ref_peps)<- as.numeric(ref_peps$Seed)

  ##Remakes the seed_peps files, to make creation of the ordered reference peps data frame a bit easier.
  ##Somewhat redundant.
  seed_peps <- read.csv("/Users/rileybrenner/Desktop/RStudio Files/STable2_41587_2013_BFnbt2585_MOESM13_ESM.csv")
  seed_peps$RT..min. <- NA

  seed_peps <- seed_peps[c(2,7)]
  ##Collects final reference peptides and their associated retention times.
  for(i in 1:nrow(seed_peps))
  {
      loc = which(ref_peps$Seed == i)
      seed_peps$Seed.Peptide.Sequence[i] <- ref_peps$Seed.Peptide.Sequence[i]
      seed_peps$RT..min.[i] <- ref_peps$RT..min.[i]
  }
}

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
  count <- t(sapply(seed_peps$Seed.Peptide.Sequence, count_each))
  rownames(count)<-c(1:nrow(count))
  count<- cbind(count, seed_peps$RT..min.)

  ms_count <- t(sapply(ms_collected_peps$Sequence, count_each))
  ms_count <- cbind(ms_count, ms_collected_peps$Spectrum.File)
  ms_count <-cbind(ms_count, ms_collected_peps$RT..min.)
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

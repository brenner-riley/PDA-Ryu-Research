seed_peps <- read.csv("/Users/rileybrenner/Desktop/RStudio Files/STable2_41587_2013_BFnbt2585_MOESM13_ESM.csv")

seed_peps$Seed.Peptide.Sequence = sapply(seed_peps$Seed.Peptide.Sequence, str_remove, pattern = "p")




setwd("/Users/rileybrenner/Desktop/HCD Data")
files <- lapply(list.files(path = "/Users/rileybrenner/Desktop/HCD Data",pattern = ".csv"), read.csv)

ms_collected_peps <- data.frame()
for(i in 1:length(files))
{
  ms_collected_peps <- rbind(files[[i]],ms_collected_peps)
}
seed_peps$RT..min. <- c()
for(i in 1:nrow(seed_peps))
{
  locs = which(ms_collected_peps$Sequence == seed_peps$Seed.Peptide.Sequence[i])
  seed_peps$RT..min.[i] <- median(ms_collected_peps$RT..min.[locs])
}




ms_collected_peps <- ms_collected_peps[-which(is.na(ms_collected_peps$q.Value)),]
ms_collected_peps <- ms_collected_peps[which(ms_collected_peps$q.Value<=0.05),]



raw_files <- ms_collected_peps$Spectrum.File%>%
  strsplit("~") %>%
  unlist() 
raw_files <- raw_files[seq(7, length(raw_files), 7)] %>%
  strsplit(".", fixed = TRUE) %>%
  unlist()
ms_collected_peps$Spectrum.File <- raw_files[seq(1, length(raw_files), 4)]






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



count <- t(sapply(seed_peps$Seed.Peptide.Sequence, count_each))
rownames(count)<-c(1:nrow(count))
count<- cbind(count, seed_peps$RT..min.)

ms_count <- t(sapply(ms_collected_peps$Sequence, count_each))
ms_count <- cbind(ms_count, ms_collected_peps$Spectrum.File)
ms_count <-cbind(ms_count, ms_collected_peps$RT..min.)


temp<- t(apply(ms_count, 1, dif_each, count))

temp_2 <- na.omit(temp)
colnames(temp_2)<-c("A","R","N","D","C","Q","E","G","H",
                    "I","L","K","M","F","P","S","T","W",
                    "Y","V","s","t","y","m","RT.change")
rownames(temp_2)<-c(1:nrow(temp_2))

lin_Reg <- lm(temp_2[,25]~temp_2[,c(1:24)])
summary(lin_Reg)
lin_Reg


confint(lin_Reg)
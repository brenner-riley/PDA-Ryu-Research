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
  seq_count <- x[c(1:24)]; seed_count<-y[1,]
  dif_count <- as.numeric(seq_count) - seed_count
  return(dif_count)
}


count <- t(sapply(seed_peps$Seed.Peptide.Sequence, count_each))
colnames(count)<-c("A","R","N","D","C","Q","E","G","H",
                   "I","L","K","M","F","P","S","T","W",
                   "Y","V","s","t","y","m")

rownames(count)<-c(1:nrow(count))


ms_count <- t(sapply(ms_collected_peps$Sequence, count_each))
ms_count <- cbind(ms_count, ms_collected_peps$Spectrum.File)

temp<- t(apply(ms_count, 1, dif_each, count))


rownames(temp)<-c(1:nrow(temp))

lin_Reg <- lm(ms_collected_peps$RT..min.~temp)
summary(lin_Reg)
lin_Reg

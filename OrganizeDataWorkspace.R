
library(ggplot2)
library(dplyr)
library(stringr)
install.packages("quantreg")     # may be needed if not already installed
library("quantreg")
library(MASS)
install.packages("ggpubr")
library(ggpubr)
library("cowplot")

sessionInfo()

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


####Retention Time Calculation
#Extra retention time of scan with the highest precursor intensity per Peptide.Sequence2
#Included Peptide.Sequence in group for the next step calculating retention time change. 
RetentionTimeByPep <- new.data %>%
  group_by(Peptide.Sequence.mod, Peptide.Sequence) %>%
  summarize(RetentionTime = Retention.Time..min.[which.max(Intensity..arb..Units.)])



write.table(RetentionTimeByPep[,c(1,3)], file='/Users/rileybrenner/Desktop/RStudio Files/RetentionTime_HCD_Marx2013_SuppT3.csv',
            sep=",",row.names = FALSE) 

#Kurtis can stop here. The below are for Riley and Alex to calculate Retention time changes. 
####Retention Time Change
#First, find retention time for unmodified peptides. 
RetentionTimeUnmod <- RetentionTimeByPep %>%
  group_by(Peptide.Sequence) %>%
  summarize(rep.rt = RetentionTime[Peptide.Sequence.mod==Peptide.Sequence])

#Second, attach retention time of unmodified peptides to corresponding modified peptides. 
Template = merge(RetentionTimeByPep, RetentionTimeUnmod, by.x="Peptide.Sequence", by.y="Peptide.Sequence", all.x = TRUE ) 
Template = Template[Template$Peptide.Sequence != Template$Peptide.Sequence.mod,] #keep only modified peptides
Template = Template[!is.na(Template$rep.rt),] #keep only modified peptides that have corresponding unmodified peptides
Template$rt.change = Template$RetentionTime-Template$rep.rt  #phospho - unmodified retention time

write.table(Template, file='/Users/rileybrenner/Desktop/RStudio Files/RetentionTime_HCD_Marx2013_SuppT3_RTchange.csv',
            sep=",",row.names = FALSE) 



####Summary Statistics


summary(Template$rt.change)
elution.window=25/60/2 #25 seconds; +-25/2 seconds
mean(Template$rt.change>0+elution.window) #0.7108509 #In publication, 70% (close enough since we don't know filtering criteria for the publication)
mean(abs(Template$rt.change)<elution.window) #0.04138363 #In publication, 4%
mean(Template$rt.change < -elution.window) #0.2477655 #In publication, 26%
  
  
  hist(Template$rt.change, breaks = 256, xlim =c(-15,15))  #Alex can make this to have smaller breaks to make it look better.
  rect(xleft=-elution.window,xright = elution.window,ybottom=0,ytop=1200, col = rgb(red = 129,blue=217 ,green=230, alpha =100, maxColorValue = 255))
  
  boxplot(Template$rt.change ~ nchar(Template$Peptide.Sequence),  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States")
  boxplot(Template$rt.change ~ nchar(Template$Peptide.Sequence), outline = F,  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States", ylim= c(-15,15))
  lines(y = c(0,27), x = c(elution.window,elution.window), lty = 2)
  lines(y = c(0,27), x = c(-elution.window,-elution.window), lty = 2)


#Now let's repeat the following but filter only for "s" modifications, 
#attach retention time of unmodified peptides to corresponding modified peptides. 
  Template_s = merge(RetentionTimeByPep, RetentionTimeUnmod, by.x="Peptide.Sequence", by.y="Peptide.Sequence", all.x = TRUE ) 
  Template_s = Template_s[(Template_s$Peptide.Sequence != Template_s$Peptide.Sequence2 & str_detect(Template_s$Peptide.Sequence2, "s")),] #keep only modified peptides
  Template_s = Template_s[!is.na(Template_s$rep.rt),] #keep only modified peptides that have corresponding unmodified peptides
  Template_s$rt.change = Template_s$RetentionTime-Template_s$rep.rt  #phospho - unmodified retention time
  
  hist(Template_s$rt.change, breaks = 256, xlim =c(-15,15))  #Alex can make this to have smaller breaks to make it look better.
  rect(xleft=-elution.window,xright = elution.window,ybottom=0,ytop=250, col = rgb(red = 129,blue=217 ,green=230, alpha =100, maxColorValue = 255))
  
  boxplot(Template_s$rt.change ~ nchar(Template_s$Peptide.Sequence),  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States")
  boxplot(Template_s$rt.change ~ nchar(Template_s$Peptide.Sequence), outline = F,  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States", ylim= c(-15,15))
  lines(y = c(0,27), x = c(elution.window,elution.window), lty = 2)
  lines(y = c(0,27), x = c(-elution.window,-elution.window), lty = 2)


#Now let's repeat the following but filter only for "t" modifications, 
#attach retention time of unmodified peptides to corresponding modified peptides. 
  Template_t = merge(RetentionTimeByPep, RetentionTimeUnmod, by.x="Peptide.Sequence", by.y="Peptide.Sequence", all.x = TRUE ) 
  Template_t = Template_t[(Template_t$Peptide.Sequence != Template_t$Peptide.Sequence2 & str_detect(Template_t$Peptide.Sequence2, "t")),] #keep only modified peptides
  Template_t = Template_t[!is.na(Template_t$rep.rt),] #keep only modified peptides that have corresponding unmodified peptides
  Template_t$rt.change = Template_t$RetentionTime-Template_t$rep.rt  #phospho - unmodified retention time
  
  
  hist(Template_t$rt.change, breaks = 256, xlim =c(-15,15))  #Alex can make this to have smaller breaks to make it look better.
  rect(xleft=-elution.window,xright = elution.window,ybottom=0,ytop=175, col = rgb(red = 129,blue=217 ,green=230, alpha =100, maxColorValue = 255))
  
  boxplot(Template_t$rt.change ~ nchar(Template_t$Peptide.Sequence),  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States")
  boxplot(Template_t$rt.change ~ nchar(Template_t$Peptide.Sequence), outline = F,  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States", ylim= c(-15,15))
  lines(y = c(0,27), x = c(elution.window,elution.window), lty = 2)
  lines(y = c(0,27), x = c(-elution.window,-elution.window), lty = 2)


#Now let's repeat the following but filter only for "y" modifications, 
#attach retention time of unmodified peptides to corresponding modified peptides. 
  Template_y = merge(RetentionTimeByPep, RetentionTimeUnmod, by.x="Peptide.Sequence", by.y="Peptide.Sequence", all.x = TRUE ) 
  Template_y = Template_y[(Template_y$Peptide.Sequence != Template_y$Peptide.Sequence2 & str_detect(Template_y$Peptide.Sequence2, "y")),] #keep only modified peptides
  Template_y = Template_y[!is.na(Template_y$rep.rt),] #keep only modified peptides that have corresponding unmodified peptides
  Template_y$rt.change = Template_y$RetentionTime-Template_y$rep.rt  #phospho - unmodified retention time
  
  hist(Template_y$rt.change, breaks = 256, xlim =c(-15,15))  #Alex can make this to have smaller breaks to make it look better.
  rect(xleft=-elution.window,xright = elution.window,ybottom=0,ytop=1200, col = rgb(red = 129,blue=217 ,green=230, alpha =100, maxColorValue = 255))
  
  boxplot(Template_y$rt.change ~ nchar(Template_y$Peptide.Sequence),  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States")
  boxplot(Template_y$rt.change ~ nchar(Template_y$Peptide.Sequence), outline = F,  horizontal = TRUE, xlab = "Change in Retention Time", ylab = "Peptide Length", main = "Box and Whisker plot for All Phosphorylation States", ylim= c(-15,15))
  lines(y = c(0,27), x = c(elution.window,elution.window), lty = 2)
  lines(y = c(0,27), x = c(-elution.window,-elution.window), lty = 2)
  
  


visualize <- function(data, up = FALSE, low = FALSE, method = "contour", Title = NULL, color = "#F2AF29")
{
  if(up & low)
  {
    data   = data[-which(Template$RetentionTime > up | Template$RetentionTime<low | Template$rep.rt > up | Template$rep.rt<low),]
  }
  
  if(!str_detect(method, "contour"))
  {
    p <- ggplot(data = Template_y, aes(x=rep.rt, y=RetentionTime))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      xlab("Retention Time of Umodified Peptide")  + 
      ylab("Retention Time of Modified Peptide") + 
      ggtitle(Title) + 
      geom_point(color = color, alpha = .1) +
      stat_smooth( method = method, color = "#2708A0", alpha = 0.75, lwd=0.5, level = F) +
      geom_abline(color="#93032E", alpha = 0.75) 
  }else{
    p<-ggplot(data = Template_y, aes(x=rep.rt, y=RetentionTime))+
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
      xlab("Retention Time of Umodified Peptide")  + 
      ylab("Retention Time of Modified Peptide") + 
      ggtitle(Title) + 
      geom_point(color = "#F2AF29", alpha = .1) +
      geom_density2d(aes(colour=..level..)) + 
      scale_colour_gradient(low="#2708A0",high="#93032E") +
      geom_abline(color="black", alpha = 0.9) 
  }
  return(p)
}
visualize_all <- function(data, up = FALSE, low = FALSE, Title = NULL, color = "#F2AF29")
{
  if(up & low)
  {
    data   = data[-which(data$RetentionTime > up | data$RetentionTime<low | data$rep.rt > up | data$rep.rt<low),]
  }
  
  loes <- ggplot(data = data, aes(x=rep.rt, y=RetentionTime))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Retention Time of Umodified Peptide")  + 
    ylab("Retention Time of Modified Peptide") + 
    ggtitle("Loess Fit")+
    geom_point(color = color, alpha = .1) +
    stat_smooth(color = "#2708A0", alpha = 0.75, lwd=0.5, level = F) +
    geom_abline(color="#93032E", alpha = 0.75) 
  
  
  cont <- ggplot(data = data, aes(x=rep.rt, y=RetentionTime))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Retention Time of Umodified Peptide")  + 
    ylab("Retention Time of Modified Peptide") +  
    ggtitle("Contour Plot")+
    geom_point(color = color, alpha = .1) +
    geom_density2d(aes(colour=..level.. ,)) + 
    scale_colour_gradient(low="#2708A0",high="red") +
    geom_abline(color="black", alpha = 0.75) 
  
  rob <- ggplot(data = data, aes(x=rep.rt, y=RetentionTime))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    xlab("Retention Time of Umodified Peptide")  + 
    ylab("Retention Time of Modified Peptide") +  
    ggtitle("Robust Regression")+
    geom_point(color = color, alpha = .1) +
    geom_smooth(method = "rlm", color = "#2708A0", alpha = 0.75, lwd=0.5) +
    geom_abline(color="#93032E", alpha = 0.75) 
  
  
  hist <-ggplot(data, aes(x=rt.change)) + 
    xlab("Change in Retention Time")  + 
    geom_histogram(color="black", fill="white",binwidth=2)
  
  temp <- nchar(data$Peptide.Sequence)
  
  box <- ggplot(data, aes(x=rt.change, y=as.factor(temp))) + 
    xlab("Change in Retention Time")  + 
    ylab("Peptide Length") + 
    geom_boxplot(fill="slateblue", alpha=0.2) 
  
  
  p <- ggdraw() +
    draw_plot(box, x = 0, y = .5, width = 1/3, height = 3/7) +
    draw_plot(hist, x = 1/3, y = .5, width = .6, height = 3/7) +
    draw_plot(loes, x = 0, y = 0, width = 1/3, height = 0.5) +
    draw_plot(rob, x = 1/3, y = 0, width = 1/3, height = 0.5) +
    draw_plot(cont, x = 2/3, y = 0, width = 1/3, height = 0.5) +
    draw_label(Title ,colour = "black", size = 18, y=.95) 
  
  return(p)
}



prelim <- function(data, Title = "Prelim All")
{
  hist <-ggplot(Template, aes(x=rt.change)) + 
    xlab("Change in Retention Time")  + 
    geom_histogram(color="black", fill="white",binwidth=2)
  
  temp <- nchar(Template$Peptide.Sequence)
  
  box <- ggplot(Template, aes(x=rt.change, y=as.factor(temp))) + 
    xlab("Change in Retention Time")  + 
    ylab("Peptide Length") + 
    geom_boxplot(fill="slateblue", alpha=0.2) 
  
  p <- ggdraw() +
    draw_plot(box, x = 0, y = 0, width = .5, height = 6/7) +
    draw_plot(hist, x = 1/2, y = 0, width = .5, height = 6/7) +
    draw_label(Title, colour = "black", size = 18, y=.95)
  return(p)
}



visualize(Template, method = "contour")
visualize_all(Template,  Title = "All Phosphorylation Points")
visualize_all(Template_s, Title = "S Phosphorylation Points", color = "#FF674D")
visualize_all(Template_t, up = 85, low = 20, Title = "T Phosphorylation Points", color = "#89B6A5")


visualize_all(Template_y, up = 85, low = 20, Title = "Y Phosphorylation Points", color = "#C191A1")





prelim(Template, "All Phosphorylation Points")
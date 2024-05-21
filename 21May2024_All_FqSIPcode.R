#Completed March 21 2024

#req packages need to go here
#change back to R Markdown

#From plotting section 

#install.packages("svglite")
#install.packages("hrbrthemes")
#install.packages("viridis")
#install.packages("mcr")
#install.packages("SimplyAgree")
#install.packages("ggridges")
#install.packages("ggbreak")
#install.packages("ggthemes")
#install.packages("egg")
#install.packages("ggtext")
#install.packages("ggpattern")
library(scales)
library(ggtext)
library(rstatix)
library(egg)
library(cowplot)
library(ggthemes)
library(forcats)
library(ggbreak)
library(ggridges)
library(svglite)
library(hrbrthemes)
library(viridis)
library(reshape2)
library(tidyverse)
library(mcr)
library(ggpubr)
library(SimplyAgree)
library(ggridges)
library(ggpubr)
library(vegan)
library(doBy)


#### 1. Internal functions: Calculating WAD, EAF, mean, sd, se, and CI ####
# Note: Make sure you are running the correct EAF eqns. In this example, 
# the isotopes that was used were 15N
# You will also need to edit the columns number according to your df. 
# Once you get there, I have also noted where you would need to make that change! 


# Weighted average density function = WAD.func
WAD.func <- function(y, x){
  WAD <- sum(x*(y/sum(y)))
  WAD
}


# sample.vec <- This function is nested into the bootstrap WAD calculation (below)
sample.vec <- function(x, ...){
  x[sample(length(x), ...)]
}


# boot.WAD.func #
boot.WAD.func <- function(X, vars=c("density.g.ml", "copies", "tube"), CI=0.90, draws=1000, size=NULL){
  
  # Create a dataframe of only x, y, and rep: 
  test.data <- data.frame(x=X[,vars[1]], y=X[,vars[2]], rep=factor(X[,vars[3]]))
  
  # Calculate observed weighted average density (WAD) for each rep:
  obs.wads <- data.frame(matrix(nrow=length(levels(test.data$rep)), ncol=2))
  names(obs.wads) <- c("wad", "rep")
  for (r in 1:length(levels(test.data$rep))){
    obs.wads$rep[r] <- levels(test.data$rep)[r]
    obs.wads$wad[r] <- WAD.func(y=test.data$y[test.data$rep == levels(test.data$rep)[r]], x=test.data$x[test.data$rep == levels(test.data$rep)[r]])
  }
  obs.wads$rep <- factor(obs.wads$rep)
  
  # Bootstrapping: Calculate a bootstrap vector of mean WADs across reps:
  if (is.null(size)){
    size <- dim(obs.wads)[1]
  }
  boot.wads <- numeric()
  for (i in 1:draws){
    boot.wads[i] <- mean(sample.vec(obs.wads$wad, size, replace=TRUE), na.rm=T)
  }
  boot.wads.CI <- quantile(boot.wads, probs=c((1-CI)/2, 1-((1-CI)/2)), na.rm=T)
  reps.NAs <- obs.wads$rep[is.na(obs.wads$wad)]
  if (length(reps.NAs) == 0){
    message <- "none"
  }
  else  message <- paste("Warning: no occurrences in rep ", paste(reps.NAs, collapse=" & "), sep="")
  
  return(list(boot.wads=boot.wads, obs.wads=obs.wads, obs.wad.mean=mean(obs.wads$wad, na.rm=T), boot.wads.mean=mean(boot.wads, na.rm=T), boot.wads.median=median(boot.wads, na.rm=T), boot.wads.CI=boot.wads.CI, message=message))
  
}


# These are used to calculate averages, SD, SE, and CI across your tubes/treatments

summed <- function(x) return(mean(x)*length(x))


trt.sum.fun <- function(x){
  ave <- mean(x, na.rm = TRUE)
  se <- (sd(x, na.rm = TRUE))/(sqrt(length(x)))
  n <- length(x[!is.na(x)])
  func.names <- list(ave = ave, se = se, n = n)
  return(sapply(func.names, unlist))
}


trt.sum.fun.CI <- function(x)  return(c(mean(x, na.rm =TRUE), 
                                        (sd(x, na.rm = TRUE))/(sqrt(length(x[!is.na(x)]))), 
                                        length(x[!is.na(x)]), 
                                        qnorm(0.95)*(sd(x, na.rm = TRUE)/(sqrt(length(x[!is.na(x)]))))
)
)

avg.se.n.fun <- function(x)  return(c(mean(x, na.rm =FALSE), (sd(x, na.rm = FALSE))/(sqrt(length(x))), length(x[!is.na(x)])))
se.sumfun <- function(x)  return(sd(x, na.rm = TRUE)/(sqrt(length(x))))
sd.sumfun <- function(x)  return(sd(x, na.rm = TRUE))


# 90% confidence interval equation
# qnorm(0.95) # = 1.644854
CI.sumfun <- function(x) return(qnorm(0.95)*(sd(x, na.rm = TRUE)/(sqrt(length(x)))))  # <- this calculates 90% C.I.
# CI.95.sumfun <- function(x) return(1.96*(sd(x, na.rm = TRUE)/(sqrt(length(x)))))  # <- this calculates 95% C.I.
len.sumfun <- function(x) return(length(x[!is.na(x)]))  


#### 2. Import data using read.table or read.csv (depending on how you saved it) ####
# Note: Before you start, see example dataframe (df) on how you need to format your data 


# Reading in df
data_all_frac <- read.csv("Lv6-filt-frac.csv")
dim(data_all_frac) #469 710
tail(colnames(data_all_frac), 7) #check col location and names of metadata; here have 6 cols


# Making a metadata based on whole tube (experimental setup)
metadata_tube <- unique(data_all_frac[,705:710])
length(unique(data_all_frac$tube)) #42 tubes


# Checking the data 
# Looks at the dimension of your dataframe
# Take note of columns your taxa start and end. 
dim(data_all_frac) # 469 (rows) x 710 (columns)
head(data_all_frac[, 1:3]) # Looks at the first three columns
head(data_all_frac[, 702:710]) # Looks at the last 8 columns 

str(unique(data_all_frac$tube)) #should be 42 tubes (14 unlabeled, 14 Lab Inc, 14 Field Inc)

# Vector of metadata names (Sample IDs, Treatment(s), density, qpcr copies number)
# Make sure you type the name exactly as how it appeared in your df you read in (data_all_frac)
data_meta_col_names<-c("index", "tube", "management", "incubation",
                       "time_hrs", "density.g.ml", "copies")


#### 3. Calculate the relative abundance for each fractions ####
# Note: Every experiment will have different numbers of taxa. I have 703 taxa. 
# For this, the taxa are located on columns 2 - 704. 

# Add column for the sum abundance of all taxa by fraction 
# This will be used to calculate the relative abundance of taxa in each fraction
data_all_frac$sum.abundance<- rowSums(data_all_frac[, 2:704], na.rm = TRUE)


# I like to look at my data and check to see if everything is calculated correctly
dim(data_all_frac) # 469 x 711, you can see that the df is now larger
head(data_all_frac[, 708:711]) # There is now a column for the sum abundance of all*(KR change from each?) taxa across row (by fractions)


# We are making a new vector of metadata names by adding sum.abundance to the data_meta_col_names
# Not all our df requires the sum.abudance column, which is why we needed to create a new vector
data_meta_col_names_sum <- c(data_meta_col_names, "sum.abundance")

# Calcalate the relative abudance of each taxa in each fraction:
data.rel.fraction <- (1/data_all_frac$sum.abundance)*data_all_frac[,!names(data_all_frac) %in% data_meta_col_names_sum]

table(is.na(data.rel.fraction))

# Finds NA or missing data and replaces it with 0
#I don't have any NA this time so skip
#data.rel.fraction <- sapply(data.rel.fraction,  FUN = function(x) {
#  x <- coalesce(x, 0)
#})

# Add metadata column to data.rel
data.rel.fraction<-cbind(data_all_frac[, names(data_all_frac) %in% data_meta_col_names_sum], data.rel.fraction)


# To double check if the relative abundance is calculated correct, 
# Calculate the porportional abundance of all taxa by fraction
# Everything should equal to 1 
data.rel.fraction$sum.abundance <- rowSums(data.rel.fraction[,!names(data.rel.fraction) %in% data_meta_col_names_sum], na.rm = TRUE)
data.rel.fraction$sum.abundance[1:50] #should be all 1's. so we are good



# Calculate number of copies per uL, based on relative abundance and total number of copies per uL
# ncopies = is the copies number (from qpcr) * relative abundance (sequence-qiime data!)
# The ncopies will be used to calculate your WAD 

ncopies <- data.rel.fraction$copies*data.rel.fraction[,!names(data.rel.fraction) %in% data_meta_col_names_sum] 
dim(ncopies) #469 703
ncopies[1:4, c(1:8, 695:703)]




# Add metadata columns to ncopies df
ncopies <- cbind(data_all_frac[, names(data_all_frac) %in% data_meta_col_names], ncopies)
dim(ncopies) #469 x 710 #This should be the same size our original 



# Pivot the ncopies and data.rel df to long format instead of wide by tube, trt.code, density, and qpcr.copies (copies)
ncopies.long <- pivot_longer(ncopies,
                             cols = c(names(ncopies)[!names(ncopies) %in% data_meta_col_names]), 
                             names_to = "taxon",
                             values_to = "ncopies") 


dim(ncopies.long) #329707 x 9

data.rel.long <- pivot_longer(data.rel.fraction, 
                              cols = c(names(data.rel.fraction)[!names(data.rel.fraction) %in% data_meta_col_names_sum]),
                              names_to = "taxon",
                              values_to = "rel.abund.fract")


dim(data.rel.long) #329707 x 10





# Now merge data.rel.long and ncopies into one df
# Do not change the data.rel.ncopies.merge df name. 
# The WAD calculation in step 5 use this dataframe
merge_col_names<-names(ncopies.long)[!names(ncopies.long) %in% "ncopies"]
data.rel.ncopies.merge <- merge(data.rel.long, ncopies.long, by = merge_col_names)  


# Save your dataframe: 
# Optional: Create a new directory and save 
# I like to keep my files organized and separate. 
# For me, it is easier to find and I'm not overwhelmed by the appt of .csv files 

dir.create("rel.abund_ncopies")
write.csv(data.rel.ncopies.merge, "rel.abund_ncopies/df.rel.ncopies.csv", row.names = FALSE)
write.csv(ncopies, "rel.abund_ncopies/ncopies.csv", row.names = FALSE)
write.csv(data.rel.fraction, "rel.abund_ncopies/rel.abundance.fraction.csv", row.names = FALSE)



#### 4. Calculate the relative abundance of each ASV in tubes ####
# This step will be used for filtering
# Due to differences in sequence depth between fractions, we will used the qpcr values to account for this

# Calculate the ncopies sum of each tube (sum all the ncopies per taxa by tube)
ncopies.tube.sum <- summaryBy(ncopies ~ tube + management + incubation + time_hrs, data = ncopies.long, FUN = sum, fun.names = "tube.sum") 
ncopies.taxa.sum <- summaryBy(ncopies ~ tube + management + incubation + time_hrs + taxon, data = ncopies.long, FUN = sum, fun.names = "taxa.sum") 

# Merge df together
ncopies.tube.taxa <- merge(ncopies.taxa.sum, ncopies.tube.sum, merge = c("tube", "management" ,"incubation", "time_hrs"))

# Calculate the relative ncopies for each taxa in each tube 
# The rel.ncopies value will be used to calculate the substrate assimilation later on (Step 12)

ncopies.tube.taxa$rel.ncopies <- (ncopies.tube.taxa$ncopies.taxa.sum / ncopies.tube.taxa$ncopies.tube.sum)


# Lets double check that our calculation is correct! 
# The total sum of the rel.ncopies of each taxa in a fraction should equal the qPCR value of that fraction  


ck.ncopies.sum <- summaryBy(ncopies.taxa.sum ~ tube, data = ncopies.tube.taxa, FUN = sum)
ck.qpcr.sum <- summaryBy(copies ~ tube, data = data_all_frac, FUN = sum)

# These values should equal to one another
#might have to remove all() or not
if(all(ck.ncopies.sum$ncopies.taxa.sum.sum == ck.qpcr.sum$qpcr.copies.sum)){
  print("Values are equal")
} else {
  print("not equal")
}


# Optional: I like to check to see if everything is calculated correctly
# I subset a taxon and the total tube/IDs should equal your sample size
# For this dataset, N = 42 

dim(subset(ncopies.tube.taxa, taxon == "d__Archaea.p__Crenarchaeota.c__Nitrososphaeria.o__Nitrososphaerales.f__Nitrososphaeraceae.__")) #42 x 8

write.csv(ncopies.tube.taxa, "rel.abund_ncopies/ncopies.tube.taxa.csv", row.names = FALSE)

#### 5. Calculating the WAD using the long-form fraction data #### 
# Make sure you ran all the functions in step 1 before you start
# This step will calculate the weighted average density (WAD) of each taxa that is present in your data

# Make an output file (This is a placeholder for the WAD output and should have nothing in it) 
wad.output<- data.frame(matrix(ncol = 3 + length(unique(data.rel.ncopies.merge)), nrow = 0))

# Depending on your computer's processing power, this will take about 1 to 2 minute to run
# The start.time and end.time is used to determine how long this set took! 
start.time <- Sys.time()
for (i in unique(data.rel.ncopies.merge$taxon)) {
  # subset the dataframe
  X <- data.rel.ncopies.merge[data.rel.ncopies.merge$taxon %in% i,]
  X.boot.wad <- boot.WAD.func(X, vars=c("density.g.ml", "ncopies", "tube"), CI=0.90, draws=1, size=NULL)
  dmc1.test.df <- as.data.frame(X.boot.wad$obs.wads)
  dmc1.test.df.mean <- as.data.frame(X.boot.wad$obs.wad.mean)
  dmc1.test.df.mean$rep <- "mean"
  names(dmc1.test.df.mean) <- c("wad", "rep")
  
  dmc1.test.df.mean <- as.data.frame(X.boot.wad$obs.wad.mean)
  dmc1.test.df.mean$rep <- "mean"
  names(dmc1.test.df.mean) <- c("wad", "rep")
  
  dmc1.test.df.5CI <- as.data.frame(X.boot.wad$boot.wads.CI[1])
  dmc1.test.df.5CI$rep <- "CI_lower"
  names(dmc1.test.df.5CI) <- c("wad", "rep")
  
  dmc1.test.df.95CI <- as.data.frame(X.boot.wad$boot.wads.CI[2])
  dmc1.test.df.95CI$rep <- "CI_upper"
  names(dmc1.test.df.95CI) <- c("wad", "rep")
  
  dmc1.test <- rbind(dmc1.test.df.mean, dmc1.test.df.5CI, dmc1.test.df.95CI, dmc1.test.df)
  rownames(dmc1.test) <- dmc1.test$rep
  
  step_result <- as.data.frame(t(dmc1.test))
  step_result <- step_result[1,]
  rownames(step_result) <- i
  
  wad.output <- rbind(wad.output,step_result)
}

end.time <- Sys.time()
end.time - start.time # 12 sec
str(wad.output)
dim(wad.output) #703 x 45 


# Convert output to numeric 
wad.output.num <- as.data.frame(sapply(wad.output, as.numeric))
dim(wad.output.num) #703 x 45
rownames(wad.output.num) <- rownames(wad.output)
wad.output.num[1:10, 1:10]

# New column for taxa 
wad.output.num$taxon <- rownames(wad.output.num)
head(wad.output.num)
dim(wad.output.num) #703 x 46

# Save output 
dir.create("WAD")
write.csv(wad.output.num, "WAD/wad.output.csv", row.names = FALSE)
wad.output.num <- read.csv("WAD/wad.output.csv")



#### 6. Filtering ASV to remove low abundance taxa ####
#  I will be removing anything that is less than 0.0001. Low-abundance taxa are still relevant to my experimental design as long as they are still present in enough tubes.
# Depending on your data, you may need to go lower or higher
# Again, because of variability in the sequence depth, we will use the rel.ncopies value to filter
# This data is located on the ncopies.tube.taxa df
#also this section will remove taxa that show up in less than 5 tubes per farm

filtered.ncopies.tube.taxa <- subset(ncopies.tube.taxa, rel.ncopies > 0.0001) #was 0.0005. changed 29 Mar 2023 for graphing purposes
length(unique(filtered.ncopies.tube.taxa$taxon)) #692 -  We initially started with 703 taxa

# Double check to see if subsetting work
min(filtered.ncopies.tube.taxa$rel.ncopies) #0.0001000203
max(filtered.ncopies.tube.taxa$rel.ncopies) #0.1196862

# String of taxa names to include in fractions
filtered.taxa <- as.vector(filtered.ncopies.tube.taxa$taxon)

# Filtered out taxa using the filtered.taxa list 
wad.filtered.low <- subset(wad.output.num, taxon %in% filtered.taxa)
length(unique(wad.filtered.low$taxon))  #692 should be same as above 



# We will remove taxa that shows up in less than x amount of samples.


#change code to pull out unlabeled tubes
#KR has to sep by farm since soils are so different
wad.unlab.conv <- select(wad.filtered.low, matches(c("taxon","^A0[0-9]")))
wad.unlab.org <- select(wad.filtered.low, matches(c("taxon","^O0[0-9]")))

# Examine the n of the unlabeled taxa 
wad.unlab.conv$n_wad.conv <- rowSums(!is.na(wad.unlab.conv[2:8]))
wad.unlab.org$n_wad.org <- rowSums(!is.na(wad.unlab.org[2:8]))

#need to get them back together now
wad.unlab.2filt <- full_join(wad.unlab.conv, wad.unlab.org, by = "taxon")
#stats that end in .x is  conv and .y is org 


# Let's remove taxa that shows up in less than 5 samples per farm.
# This value is arbitrary and based on your experimental design
#kinsey: 2 farms x 7 reps = 14. Note that the incubations aren't included as a treatment here since don't have their own sep time 0.

wad.unlab.tube.filt <- subset(wad.unlab.2filt, n_wad.conv >= 5 & n_wad.org >= 5) ##5 TUBES
length(unique(wad.unlab.tube.filt$taxon)) #513 left (previously 692)

# Double check to see if filtering work
min(wad.unlab.tube.filt$n_wad.conv) #5
min(wad.unlab.tube.filt$n_wad.org) #5

taxa2filt <- as.vector(wad.unlab.tube.filt$taxon)  


wad.unlab <- select(wad.filtered.low, matches(c("taxon","^[A-Z]0[0-9]")))
length(unique(wad.unlab$taxon)) #original 692 check
wad.unlab <- subset(wad.filtered.low, taxon %in% taxa2filt)
length(unique(wad.unlab$taxon)) #now 513



# Add the relative ncopies of each taxon in the unlabeled samples
# The rel_abundance.sum is the sum of the relative abundance of taxon A in tube A (summing all the fractions)
# Lets subset the unlabeled samples from the ncopies.tube.taxa 
rel.ncopies.unlab <- subset(ncopies.tube.taxa, incubation = "None")

# Lets average out the relative abundance of each taxa 
ave.ncopies.unlab <- summaryBy(rel.ncopies ~ taxon, data = rel.ncopies.unlab, 
                               FUN = trt.sum.fun, fun.names = c("mean", "sd", "n"))


# Merge the wad.unlab with ave.ncopies.unlab
wad.unlab.abund <- merge(wad.unlab, ave.ncopies.unlab)


# String of filtered taxa names to subset fractions
contingency.taxa <- as.vector(wad.unlab.abund$taxon)   


# Optional: Determine the percentage of community that was retained after filtering
sum(wad.unlab.abund$rel.ncopies.ave) #98.6%


#### 7. Data organization and visualization #### 

# Transpose the df. Remove  the first 3 columns: mean and CIs 
wad.transposed <- as.data.frame(t(wad.output.num[, -c(1:3)]))

# The taxon is located on row 43. Convert this to column name
# Depending on your dataframe, this might be located in a different position
# Also, convert the tube row into row.names
colnames(wad.transposed) <- wad.transposed[43, ]
wad.transposed$tube <- row.names(wad.transposed)

# Add metadata to the df by merging the metadata df = metadata_tube with wad.transposed by tube ID
wad.transposed <- merge(metadata_tube, wad.transposed, by = "tube")

# Check the df to see if it looks correct. 
view(wad.transposed[1:3, 1:10])
dim(wad.transposed) # 469 x 709


# Convert df into long format for filtering proposes 
wad.long <- pivot_longer(wad.transposed, cols = c(names(wad.transposed[!names(wad.transposed) %in% data_meta_col_names])),
                         names_to = "taxon", values_to = "WAD")

# Convert WAD value from characters to numeric
wad.long$WAD <- as.numeric(wad.long$WAD)

# Filtered the wad.long df to keep desired taxa
#i used 5 tubes to start. filtered for each farm individually
wad.long.filtered <- subset(wad.long, taxon %in% contingency.taxa)



# Check to see if filtering work should be same as before
length(unique(wad.long.filtered$taxon)) #513 now which is correct

# Save 
write.csv(wad.long, "WAD/new.wad.long.csv", row.names = FALSE)
write.csv(wad.long.filtered, "WAD/new.wad.long.filtered.csv", row.names = FALSE)


# Subset the df by tree x treatment 
wad.OF <- subset(wad.long.filtered, management == "Organic" & incubation %in% c("Field", "None"))
wad.OL <- subset(wad.long.filtered, management == "Organic" & incubation %in% c("Lab", "None"))
wad.CF <- subset(wad.long.filtered, management == "Conventional" & incubation %in% c("Field", "None"))
wad.CL <- subset(wad.long.filtered, management == "Conventional" & incubation %in% c("Lab", "None"))



# Let's look at the mean WAD 
aggregate(wad.OF$WAD, list(wad.OF$time_hrs), FUN = mean, na.rm = TRUE) 
aggregate(wad.OL$WAD, list(wad.OL$time_hrs), FUN = mean, na.rm = TRUE) 
aggregate(wad.CF$WAD, list(wad.CF$time_hrs), FUN = mean, na.rm = TRUE) 
aggregate(wad.CL$WAD, list(wad.CL$time_hrs), FUN = mean, na.rm = TRUE) 



# Lets plot the WAD . 
# In this step, ideally you want the mean WAD of each treatment to be similar. This won't always happen

plot.wad.OF <- ggplot(wad.OF, aes(x = tube, y = WAD))+
  geom_boxplot()+
  geom_violin()+
  geom_hline(yintercept= 1.703774)+
  ggtitle("OF WAD")+
  theme_classic()+
  theme(legend.position = "top")
plot.wad.OF
#keep everything, though not super consistent
ggsave("WAD/new.wad.OF.png")

plot.wad.OL <- ggplot(wad.OL, aes(x = tube, y = WAD))+
  geom_boxplot()+
  geom_violin()+
  geom_hline(yintercept= 1.707330)+
  ggtitle("OL WAD")+
  theme_classic()+
  theme(legend.position = "top")
plot.wad.OL
ggsave("WAD/new.wad.OL.png")



plot.wad.CF <- ggplot(wad.CF, aes(x = tube, y = WAD))+
  geom_boxplot()+
  geom_violin()+
  geom_hline(yintercept= 1.706305)+
  ggtitle("CF WAD")+
  theme_classic()+
  theme(legend.position = "top")
plot.wad.CF
ggsave("WAD/new.wad.CF.png")


plot.wad.CL <- ggplot(wad.CL, aes(x = tube, y = WAD))+
  geom_boxplot()+
  geom_violin()+
  geom_hline(yintercept= 1.708841)+
  ggtitle("CL WAD")+
  theme_classic()+
  theme(legend.position = "top")
plot.wad.CL
ggsave("WAD/new.wad.CL.png")


# Optional: Save plots
library(ggpubr)
plot.all.wad <- ggarrange(plot.wad.OF, plot.wad.OL, plot.wad.CF, plot.wad.CL, ncol = 2, nrow = 2)
plot.all.wad
ggsave("WAD/plot.all.wad.png")




# Use the trt.sum.fun.CI function to calculate the WAD mean, std error, number of reps, and CI
# Depending on how you want to calculate the WAD shift, 
# you can average all unlabeled WAD (regardless of treatment) and use that value or 
# use the unlabeled mean WAD to calculate shift from that treatment: 
# So the unlabeled from Organic is only use to calculate shift from Organic group and 
# Conventional unlabeled is only use to calculate shift from Conv. group

# For this example, I am going to average all the WAD unlabeled together
# To do so, I need to subset the unlabeled (time_hrs == 0) and row bind together 

OF.WAD.unlab <- subset(wad.OF, time_hrs == "0")
OL.WAD.unlab <- subset(wad.OL, time_hrs == "0")
CF.WAD.unlab <- subset(wad.CF, time_hrs == "0")
CL.WAD.unlab <- subset(wad.CL, time_hrs == "0")

wad.unlab.all <- as.data.frame(rbind(OF.WAD.unlab, OL.WAD.unlab, CF.WAD.unlab, CL.WAD.unlab))


# Removing the unlab from the df 
OF.WAD.lab <- subset(wad.OF, time_hrs != "0")
OL.WAD.lab <- subset(wad.OL, time_hrs != "0")
CF.WAD.lab <- subset(wad.CF, time_hrs != "0")
CL.WAD.lab <- subset(wad.CL, time_hrs != "0")


# Lets calculate the WAD mean for each treatment
wad.all.sum <- summaryBy(WAD ~ taxon, data = wad.unlab.all, FUN = trt.sum.fun.CI)
wad.OF.sum <- summaryBy(WAD ~ taxon, data = OF.WAD.lab, FUN = trt.sum.fun.CI)
wad.OL.sum <- summaryBy(WAD ~ taxon, data = OL.WAD.lab, FUN = trt.sum.fun.CI)
wad.CF.sum <- summaryBy(WAD ~ taxon, data = CF.WAD.lab, FUN = trt.sum.fun.CI)
wad.CL.sum <- summaryBy(WAD ~ taxon, data = CL.WAD.lab, FUN = trt.sum.fun.CI)


# Rename the column header 
# Do not change the 'mean_WAD' name. The EAF calculation uses this column name!  
wad.sum.name <- c("taxon", "mean_WAD", "std_err", "n_reps", "CI")
names(wad.OF.sum) <- wad.sum.name
names(wad.OL.sum) <- wad.sum.name
names(wad.CF.sum) <- wad.sum.name
names(wad.CL.sum) <- wad.sum.name

# Rename the unlabeled columns to distinguish from the labeled 
# Do not change the mean_WAD_unlab column name. This is used in the EAF calculation function
wad.unlab.name <- c("taxon", "mean_WAD_unlab", "std_err_unlab", "n_reps_unlab", "CI_unlab")
names(wad.all.sum) <- wad.unlab.name


# Lets merge the df wad.all.sum with each treatment WAD df 
wad.OF.shift <- merge(wad.all.sum, wad.OF.sum, by = "taxon")
wad.OL.shift <- merge(wad.all.sum, wad.OL.sum, by = "taxon")
wad.CF.shift <- merge(wad.all.sum, wad.CF.sum, by = "taxon")
wad.CL.shift <- merge(wad.all.sum, wad.CL.sum, by = "taxon")

# Plot the WAD ordered by taxa 
plot.wad.OF.taxa <- ggplot(wad.OF.shift) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD), ymin = mean_WAD - std_err, ymax = mean_WAD + std_err), width=.2) +
  geom_point(aes(x = reorder(taxon, mean_WAD, mean), y = mean_WAD), size = 4, color = "black", shape = 21) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD_unlab), ymin = mean_WAD_unlab - std_err_unlab, ymax = mean_WAD_unlab + std_err_unlab, color = "Unlabeled"))+
  geom_point(aes(x = reorder(taxon, mean_WAD_unlab), y = mean_WAD_unlab, fill = "Unlabeled"), size = 3, color = "black", shape = 23, alpha = 0.5) +
  theme_classic2()+
  theme(legend.position = "top")+
  ggtitle("OF WAD")
plot.wad.OF.taxa
#note shift really small here, probably because had poor uptake of 15N due to tillage and increases drainage/evaporation
ggsave("WAD/wad.OF.taxa.png")

plot.wad.OL.taxa <- ggplot(wad.OL.shift) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD), ymin = mean_WAD - std_err, ymax = mean_WAD + std_err), width=.2) +
  geom_point(aes(x = reorder(taxon, mean_WAD, mean), y = mean_WAD), size = 4, color = "black", shape = 21) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD_unlab), ymin = mean_WAD_unlab - std_err_unlab, ymax = mean_WAD_unlab + std_err_unlab, color = "Unlabeled"))+
  geom_point(aes(x = reorder(taxon, mean_WAD_unlab), y = mean_WAD_unlab, fill = "Unlabeled"), size = 3, color = "black", shape = 23, alpha = 0.5) +
  theme_classic2()+
  theme(legend.position = "top")+
  ggtitle("OL WAD")
plot.wad.OL.taxa
#good shift
ggsave("WAD/wad.OL.taxa.png")

plot.wad.CF.taxa <- ggplot(wad.CF.shift) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD), ymin = mean_WAD - std_err, ymax = mean_WAD + std_err), width=.2) +
  geom_point(aes(x = reorder(taxon, mean_WAD, mean), y = mean_WAD), size = 4, color = "black", shape = 21) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD_unlab), ymin = mean_WAD_unlab - std_err_unlab, ymax = mean_WAD_unlab + std_err_unlab, color = "Unlabeled"))+
  geom_point(aes(x = reorder(taxon, mean_WAD_unlab), y = mean_WAD_unlab, fill = "Unlabeled"), size = 3, color = "black", shape = 23, alpha = 0.5) +
  theme_classic2()+
  theme(legend.position = "top")+
  ggtitle("CF WAD")
plot.wad.CF.taxa
#pretty good shift
ggsave("WAD/wad.CF.taxa.png")

plot.wad.CL.taxa <- ggplot(wad.CL.shift) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD), ymin = mean_WAD - std_err, ymax = mean_WAD + std_err), width=.2) +
  geom_point(aes(x = reorder(taxon, mean_WAD, mean), y = mean_WAD), size = 4, color = "black", shape = 21) +
  geom_errorbar(aes(x = reorder(taxon, mean_WAD_unlab), ymin = mean_WAD_unlab - std_err_unlab, ymax = mean_WAD_unlab + std_err_unlab, color = "Unlabeled"))+
  geom_point(aes(x = reorder(taxon, mean_WAD_unlab), y = mean_WAD_unlab, fill = "Unlabeled"), size = 3, color = "black", shape = 23, alpha = 0.5) +
  theme_classic2()+
  theme(legend.position = "top")+
  ggtitle("CL WAD")
plot.wad.CL.taxa
#great
ggsave("WAD/wad.CL.taxa.png")


# Optional: Save plots 
plot.all.wad.taxa <- ggarrange(plot.wad.OF.taxa, plot.wad.OL.taxa, plot.wad.CF.taxa, plot.wad.CL.taxa, ncol = 2, nrow = 2)
ggsave("WAD/new.plot.all.wad.taxa.png", plot.all.wad.taxa)



#### 9. Calculate WAD shift #### 
wad.OF.shift <- mutate(wad.OF.shift, shift = mean_WAD - mean_WAD_unlab)
wad.OL.shift <- mutate(wad.OL.shift, shift = mean_WAD - mean_WAD_unlab)
wad.CF.shift <- mutate(wad.CF.shift, shift = mean_WAD - mean_WAD_unlab)
wad.CL.shift <- mutate(wad.CL.shift, shift = mean_WAD - mean_WAD_unlab)

# Order shift (small to large) and add observation and NAs count column 
wad.OF.shift <- wad.OF.shift[order(wad.OF.shift$shift), ]
wad.OF.shift$observation <- 1:nrow(wad.OF.shift)
wad.OF.shift$nas <- rowSums(is.na(wad.OF.shift))

wad.OL.shift <- wad.OL.shift[order(wad.OL.shift$shift), ]
wad.OL.shift$observation <- 1:nrow(wad.OL.shift)
wad.OL.shift$nas <- rowSums(is.na(wad.OL.shift))


wad.CF.shift <- wad.CF.shift[order(wad.CF.shift$shift), ]
wad.CF.shift$observation <- 1:nrow(wad.CF.shift)
wad.CF.shift$nas <- rowSums(is.na(wad.CF.shift))


wad.CL.shift <- wad.CL.shift[order(wad.CL.shift$shift), ]
wad.CL.shift$observation <- 1:nrow(wad.CL.shift)
wad.CL.shift$nas <- rowSums(is.na(wad.CL.shift))


# Save 
write.csv(wad.OF.shift, "WAD/new.wad.OF.shift.csv", row.names = FALSE)
write.csv(wad.OL.shift, "WAD/new.wad.OL.shift.csv", row.names = FALSE)
write.csv(wad.CF.shift, "WAD/new.wad.CF.shift.csv", row.names = FALSE)
write.csv(wad.CL.shift, "WAD/new.wad.CL.shift.csv", row.names = FALSE)


# qplot of averaged WAD 
# First plot: Ideally, you are looking for a linear relationship between the unlab vs lab WAD. 
# Also, you want most of your points to fall below the linear regression line
# Second plot: Distribution of the shifts 
# Third plot: A second a way of looking at the distribution of shift
# Ideally, you want most of your observation to be above zero
qplot(mean_WAD, mean_WAD_unlab, data = wad.OF.shift) +
  geom_abline(reg=lm(mean_WAD_unlab~mean_WAD, data = wad.OF.shift))
ggsave("WAD/qplot1_OF.png")
qplot(shift, data = wad.OF.shift) 
ggsave("WAD/qplot2shift_OF.png")
qplot(observation, shift, data = wad.OF.shift)+
  geom_abline(intercept = 0, slope = 0)
ggsave("WAD/qplot3shift_OF.png")

qplot(mean_WAD, mean_WAD_unlab, data = wad.OL.shift)+
  geom_abline(reg=lm(mean_WAD_unlab~mean_WAD, data = wad.OL.shift))
ggsave("WAD/qplot1_OL.png")
qplot(shift, data = wad.OL.shift)
ggsave("WAD/qplot2shift_OL.png")
qplot(observation, shift, data = wad.OL.shift)+
  geom_abline(intercept = 0, slope = 0)
ggsave("WAD/qplot3shift_OL.png")

qplot(mean_WAD, mean_WAD_unlab, data = wad.CF.shift)+
  geom_abline(reg=lm(mean_WAD_unlab~mean_WAD, data = wad.CF.shift))
ggsave("WAD/qplot1_CF.png")
qplot(shift, data = wad.CF.shift)
ggsave("WAD/qplot2shift_CF.png")
qplot(observation, shift, data = wad.CF.shift)+
  geom_abline(intercept = 0, slope = 0)
ggsave("WAD/qplot3shift_CF.png")

qplot(mean_WAD, mean_WAD_unlab, data = wad.CL.shift)+
  geom_abline(reg=lm(mean_WAD_unlab~mean_WAD, data = wad.CL.shift))
ggsave("WAD/qplot1_CL.png")
qplot(shift, data = wad.CL.shift)
ggsave("WAD/qplot2shift_CL.png")
qplot(observation, shift, data = wad.CL.shift)+
  geom_abline(intercept = 0, slope = 0)
ggsave("WAD/qplot3shift_CL.png")


# Calculate the average WAD and shift of each taxon by trt (tube replicates)
# Merge WAD.XX df with the XX.unlab df by taxon 
wad.OF.tube <- merge(wad.OF, wad.OF.shift, by = "taxon")
wad.OL.tube <- merge(wad.OL, wad.OL.shift, by = "taxon")
wad.CF.tube <- merge(wad.CF, wad.CF.shift, by = "taxon")
wad.CL.tube <- merge(wad.CL, wad.CL.shift, by = "taxon")



# Convert to wide format with tube ids as column headers 
### without "values_fn = median" makes a list of the same value for each taxon

wad.OF.wide <- pivot_wider(wad.OF.tube, id_cols = c("taxon", "mean_WAD_unlab", "std_err_unlab", "n_reps_unlab"),
                           names_from = "tube", values_from = "WAD", values_fn = median)

dim(wad.OF.wide) #513 x 18



wad.OL.wide <- pivot_wider(wad.OL.tube, id_cols = c("taxon", "mean_WAD_unlab", "std_err_unlab", "n_reps_unlab"),
                           names_from = "tube", values_from = "WAD", values_fn = median)

wad.CF.wide <- pivot_wider(wad.CF.tube, id_cols = c("taxon", "mean_WAD_unlab", "std_err_unlab", "n_reps_unlab"),
                           names_from = "tube", values_from = "WAD", values_fn = median)

wad.CL.wide <- pivot_wider(wad.CL.tube, id_cols = c("taxon", "mean_WAD_unlab", "std_err_unlab", "n_reps_unlab"),
                           names_from = "tube", values_from = "WAD", values_fn = median)



#get rid of time 0 tubes

wad.OF.B <- select(wad.OF.wide, contains(c("taxon", "unlab", "O5")))
dim(wad.OF.B) #513 x 11 col, 7 reps 4 other
wad.OL.B <- select(wad.OL.wide, contains(c("taxon", "unlab", "O5")))
wad.CF.B <- select(wad.CF.wide, contains(c("taxon", "unlab", "A5")))
wad.CL.B <- select(wad.CL.wide, contains(c("taxon", "unlab", "A5")))



# Calculate the average WAD of each taxa by tube 
# Row 5:11 is based on this dataset with 7 reps. Change the row values to suit your df 

wad.OF.B$avetrt <- rowMeans(wad.OF.B[5:11], na.rm = TRUE)
wad.OL.B$avetrt <- rowMeans(wad.OL.B[5:11], na.rm = TRUE)
wad.CF.B$avetrt <- rowMeans(wad.CF.B[5:11], na.rm = TRUE)
wad.CL.B$avetrt <- rowMeans(wad.CL.B[5:11], na.rm = TRUE)


# Calculate the shift for each taxa 
wad.OF.B <- mutate(wad.OF.B, shift = avetrt - mean_WAD_unlab)
wad.OL.B <- mutate(wad.OL.B, shift = avetrt - mean_WAD_unlab)
wad.CF.B <- mutate(wad.CF.B, shift = avetrt - mean_WAD_unlab)
wad.CL.B <- mutate(wad.CL.B, shift = avetrt - mean_WAD_unlab)

# Ordered df by shift (small to large) and add observation columns
wad.OF.B <- wad.OF.B[order(wad.OF.B$shift), ]
wad.OF.B$observation <- 1:nrow(wad.OF.B)
wad.OL.B <- wad.OL.B[order(wad.OL.B$shift), ]
wad.OL.B$observation <- 1:nrow(wad.OL.B)
wad.CF.B <- wad.CF.B[order(wad.CF.B$shift), ]
wad.CF.B$observation <- 1:nrow(wad.CF.B)
wad.CL.B <- wad.CL.B[order(wad.CL.B$shift), ]
wad.CL.B$observation <- 1:nrow(wad.CL.B)


# Save
write.csv(wad.OF.B, "WAD/new.wad.OF.B.csv", row.names = FALSE)
write.csv(wad.OL.B, "WAD/new.wad.OL.B.csv", row.names = FALSE)
write.csv(wad.CF.B, "WAD/new.wad.CF.B.csv", row.names = FALSE)
write.csv(wad.CL.B, "WAD/new.wad.CL.B.csv", row.names = FALSE)




#### 10. Bottom tube correction ####
# For this data set, I am going to use the bottom 50 based on shift (low to high) 
# This value is completely arbitrary and based on how your data 
# Convert the observation column into numeric value 

wad.OF.B$observation <- as.numeric(wad.OF.B$observation) 
wad.OL.B$observation <- as.numeric(wad.OL.B$observation) 
wad.CF.B$observation <- as.numeric(wad.CF.B$observation) 
wad.CL.B$observation <- as.numeric(wad.CL.B$observation) 

# Pull out the bottom 50 - This will be used to generate the offset value for each tube correction
#also tested 25, 50, and 100. helped some but according to data the shifts stop being negative around row #17.
wad.OF.B.bottom <- wad.OF.B[1:50, ]
wad.OL.B.bottom <- wad.OL.B[1:50, ]
wad.CF.B.bottom <- wad.CF.B[1:50, ]
wad.CL.B.bottom <- wad.CL.B[1:50, ]


# Calculating the difference between WAD of each taxa x tube
# and the mean unlabeled WAD of each taxa 
# To get the offset values for each tube, calculate the median of the difference
# Convert to df 

# The WAD of each tube is located on columns 5-8. This will vary with every experiment/dataset
diff.OF.B <- wad.OF.B.bottom[, 5:11] - wad.OF.B.bottom$mean_WAD_unlab
diff.OL.B <- wad.OL.B.bottom[, 5:11] - wad.OL.B.bottom$mean_WAD_unlab
diff.CF.B <- wad.CF.B.bottom[, 5:11] - wad.CF.B.bottom$mean_WAD_unlab
diff.CL.B <- wad.CL.B.bottom[, 5:11] - wad.CL.B.bottom$mean_WAD_unlab


off.OF.B <- as.data.frame(rbind(apply(diff.OF.B, 2, median, na.rm = TRUE))) #2 indicates the function will be applied over the columns
off.OL.B <- as.data.frame(rbind(apply(diff.OL.B, 2, median, na.rm = TRUE)))
off.CF.B <- as.data.frame(rbind(apply(diff.CF.B, 2, median, na.rm = TRUE)))
off.CL.B <- as.data.frame(rbind(apply(diff.CL.B, 2, median, na.rm = TRUE)))

# Use the offset values to do the tube correction. This is just for the bottom 50. 

#created a df with the offset values repeated so that it equals the same size as the WAD df
#This way, I can subtract the two df from each other and not have to do each tube individually


off.OF.B <- sapply(off.OF.B, rep.int, times = nrow(wad.OF.B))
off.OL.B <- sapply(off.OL.B, rep.int, times = nrow(wad.OL.B))
off.CF.B <- sapply(off.CF.B, rep.int, times = nrow(wad.CF.B))
off.CL.B <- sapply(off.CL.B, rep.int, times = nrow(wad.CL.B))

# Tube correction: Tube WAD of each taxa - tube offset value

cor.OF.B <- wad.OF.B[, 5:11] - off.OF.B 
cor.OL.B <- wad.OL.B[, 5:11] - off.OL.B 
cor.CF.B <- wad.CF.B[, 5:11] - off.CF.B 
cor.CL.B <- wad.CL.B[, 5:11] - off.CL.B 

# Change the tube correction column names 
colnames(cor.OF.B) <- paste("x", colnames(cor.OF.B), sep = "")
colnames(cor.OL.B) <- paste("x", colnames(cor.OL.B), sep = "")
colnames(cor.CF.B) <- paste("x", colnames(cor.CF.B), sep = "")
colnames(cor.CL.B) <- paste("x", colnames(cor.CL.B), sep = "")


# cbind the tube corrected column with the uncorrected tube 
wad.OF.C <- cbind(wad.OF.B, cor.OF.B)
wad.OL.C <- cbind(wad.OL.B, cor.OL.B)
wad.CF.C <- cbind(wad.CF.B, cor.CF.B)
wad.CL.C <- cbind(wad.CL.B, cor.CL.B)

# Add summary for all corrected reps
# This will use the last four columns, containing the corrected WAD
dim(wad.OF.C) #513:21, so target columns will be be 15:21 (there are 7 reps)

wad.OF.C$avetrt.c <- rowMeans(wad.OF.C[15:21], na.rm = TRUE)
wad.OL.C$avetrt.c <- rowMeans(wad.OL.C[15:21], na.rm = TRUE)
wad.CF.C$avetrt.c <- rowMeans(wad.CF.C[15:21], na.rm = TRUE)
wad.CL.C$avetrt.c <- rowMeans(wad.CL.C[15:21], na.rm = TRUE)

# Save: 
write.csv(wad.OF.C, "WAD/new.corrected.wad.OF.C.csv", row.names = FALSE)
write.csv(wad.OL.C, "WAD/new.corrected.wad.OL.C.csv", row.names = FALSE)
write.csv(wad.CF.C, "WAD/new.corrected.wad.CF.C.csv", row.names = FALSE)
write.csv(wad.CL.C, "WAD/new.corrected.wad.CL.C.csv", row.names = FALSE)

# Plot the WAD of each tube vs mean_WAD_unlab. 
# This is optional: But in this example, we will plot both the uncorrected vs the corrected. 
# Ideally, the corrected WAD plot should be tighter than the uncorrected.
# What you are looking: Reps to be on top of (or as close to) one another
# Lower mean_WAD_unlab in comparison to your labeled WAD
# You want to be able to see pattern similar to your DNA curves 

plot.OF.C <- ggplot(wad.OF.C, aes(x = observation)) +
  geom_line(aes(y = O51F, color = "O51F"), size = 0.25)+
  geom_line(aes(y = O52F, color = "O52F"), size = 0.25)+
  geom_line(aes(y = O53F, color = "O53F"), size = 0.25)+
  geom_line(aes(y = O54F, color = "O54F"), size = 0.25)+
  geom_line(aes(y = O55F, color = "O55F"), size = 0.25)+
  geom_line(aes(y = O56F, color = "O56F"), size = 0.25)+
  geom_line(aes(y = O57F, color = "O57F"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.5)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("original OF")
plot.OF.C

plot.OF.C.cor <- ggplot(wad.OF.C, aes(x = observation)) +
  geom_line(aes(y = xO51F, color = "xO51F"), size = 0.25)+
  geom_line(aes(y = xO52F, color = "xO52F"), size = 0.25)+
  geom_line(aes(y = xO53F, color = "xO53F"), size = 0.25)+
  geom_line(aes(y = xO54F, color = "xO54F"), size = 0.25)+
  geom_line(aes(y = xO55F, color = "xO55F"), size = 0.25)+
  geom_line(aes(y = xO56F, color = "xO56F"), size = 0.25)+
  geom_line(aes(y = xO57F, color = "xO57F"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.5)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("corrected OF - 50")
plot.OF.C.cor

plot.OL.C <- ggplot(wad.OL.C, aes(x = observation)) +
  geom_line(aes(y = O51L, color = "O51L"), size = 0.25)+
  geom_line(aes(y = O52L, color = "O52L"), size = 0.25)+
  geom_line(aes(y = O53L, color = "O53L"), size = 0.25)+
  geom_line(aes(y = O54L, color = "O54L"), size = 0.25)+
  geom_line(aes(y = O55L, color = "O55L"), size = 0.25)+
  geom_line(aes(y = O56L, color = "O56L"), size = 0.25)+
  geom_line(aes(y = O57L, color = "O57L"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.75)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("original OL")
plot.OL.C

plot.OL.C.cor <- ggplot(wad.OL.C, aes(x = observation)) +
  geom_line(aes(y = xO51L, color = "xO51L"), size = 0.25)+
  geom_line(aes(y = xO52L, color = "xO52L"), size = 0.25)+
  geom_line(aes(y = xO53L, color = "xO53L"), size = 0.25)+
  geom_line(aes(y = xO54L, color = "xO54L"), size = 0.25)+
  geom_line(aes(y = xO55L, color = "xO55L"), size = 0.25)+
  geom_line(aes(y = xO56L, color = "xO56L"), size = 0.25)+
  geom_line(aes(y = xO57L, color = "xO57L"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.75)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("corrected OL - 50")
plot.OL.C.cor
##O51L is a problem?

plot.CF.C <- ggplot(wad.CF.C, aes(x = observation)) +
  geom_line(aes(y = A51F, color = "A51F"), size = 0.25)+
  geom_line(aes(y = A52F, color = "A52F"), size = 0.25)+
  geom_line(aes(y = A53F, color = "A53F"), size = 0.25)+
  geom_line(aes(y = A54F, color = "A54F"), size = 0.25)+
  geom_line(aes(y = A55F, color = "A55F"), size = 0.25)+
  geom_line(aes(y = A56F, color = "A56F"), size = 0.25)+
  geom_line(aes(y = A57F, color = "A57F"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.75)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("original CF")
plot.CF.C

plot.CF.C.cor <- ggplot(wad.CF.C, aes(x = observation)) +
  geom_line(aes(y = xA51F, color = "xA51F"), size = 0.25)+
  geom_line(aes(y = xA52F, color = "xA52F"), size = 0.25)+
  geom_line(aes(y = xA53F, color = "xA53F"), size = 0.25)+
  geom_line(aes(y = xA54F, color = "xA54F"), size = 0.25)+
  geom_line(aes(y = xA55F, color = "xA55F"), size = 0.25)+
  geom_line(aes(y = xA56F, color = "xA56F"), size = 0.25)+
  geom_line(aes(y = xA57F, color = "xA57F"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.75)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("corrected CF - 50")
plot.CF.C.cor

plot.CL.C <- ggplot(wad.CL.C, aes(x = observation)) +
  geom_line(aes(y = A51L, color = "A51L"), size = 0.25)+
  geom_line(aes(y = A52L, color = "A52L"), size = 0.25)+
  geom_line(aes(y = A53L, color = "A53L"), size = 0.25)+
  geom_line(aes(y = A54L, color = "A54L"), size = 0.25)+
  geom_line(aes(y = A55L, color = "A55L"), size = 0.25)+
  geom_line(aes(y = A56L, color = "A56L"), size = 0.25)+
  geom_line(aes(y = A57L, color = "A57L"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.75)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("original CL")
plot.CL.C

plot.CL.C.cor <- ggplot(wad.CL.C, aes(x = observation)) +
  geom_line(aes(y = xA51L, color = "xA51L"), size = 0.25)+
  geom_line(aes(y = xA52L, color = "xA52L"), size = 0.25)+
  geom_line(aes(y = xA53L, color = "xA53L"), size = 0.25)+
  geom_line(aes(y = xA54L, color = "xA54L"), size = 0.25)+
  geom_line(aes(y = xA55L, color = "xA55L"), size = 0.25)+
  geom_line(aes(y = xA56L, color = "xA56L"), size = 0.25)+
  geom_line(aes(y = xA57L, color = "xA57L"), size = 0.25)+
  geom_line(aes(y = mean_WAD_unlab), linetype = "dashed", size = 0.75)+
  theme_classic2()+
  theme(legend.position = "top") + 
  ylab("WAD")+
  ggtitle("corrected C-50")
plot.CL.C.cor



# Optional: Save plot 
plot.OF.comp <- ggarrange(plot.OF.C, plot.OF.C.cor, nrow = 2, ncol = 1)
plot.OL.comp <- ggarrange(plot.OL.C, plot.OL.C.cor, nrow = 2, ncol = 1)
plot.CF.comp <- ggarrange(plot.CF.C, plot.CF.C.cor, nrow = 2, ncol = 1)
plot.CL.comp <- ggarrange(plot.CL.C, plot.CL.C.cor, nrow = 2, ncol = 1)


ggsave("WAD/new.plot.OF.comp.50.png", plot.OF.comp)
ggsave("WAD/new.plot.OL.comp.50.png", plot.OL.comp)
ggsave("WAD/new.plot.CF.comp.50.png", plot.CF.comp)
ggsave("WAD/new.plot.CL.comp.50.png", plot.CL.comp)


#### 11. Calculating EAF values #### 
# The columns that you need for the eaf function: mean_WAD_unlab, and corrected tube WAD (change column number to match yours in the eaf function)
# Make sure you grab the right eqn based on the isotopes you used!

eaf_eqn <-function (df){
  df.func <- df  
  # Gi = Calculating the GC content in taxon i   
  df.func$GC <- (1/0.083506) * (df.func$mean_WAD_unlab - 1.646057)
  # M. Light i = observed mol wt of the DNA fragment containing the 16S RNA gene for taxon i in the unlab treatment (natural abundance) (g*mol-1)
  df.func$MW_unlabeled <- (0.496 * df.func$GC) + 307.691
  # M. Heavy Max i = Theoretical mol wt of the DNA fragment containing the 16S RNA gene for taxon i assuming maximum labeling by the heavy isotopes (g * mol-1)
  df.func$MW_max<- (0.5024851 * df.func$GC) + 3.517396 + df.func$MW_unlabeled
  #0.5 changed for 15N... is it supposed to be negative? Chansos had a - in front of the C/Ovalue
  
  # Get the difference of the labeled and unlabeled WAD 
  # Note: df will vary with each experiment. the first df is the tube corrected columns! 
  # The second df is the mean_WAD_unlab
  df.diff <- as.data.frame(sapply (df[,15:21], function (x) x - df[, 2]))
  # Change the column names of the diff columns 
  colnames(df.diff) <- paste(colnames(df)[c(15:21)], "diff", sep = "_")
  
  # M. Lab i = Observed mol wt of the DNA fragment containing the 16S RNA gene for taxon i in the labeled treatment (g * mol−1)
  # This is calculating the MW_labled = ((WAD_Differences/WAD_unlabeled) + 1) * MW_unlabled 
  df.MW_lab <- as.data.frame(sapply(df.diff, function(x) (((x/df.func$mean_WAD_unlab) + 1) * df.func$MW_unlabeled)))
  # Change the column names of MW_labeled columns 
  colnames(df.MW_lab) <- paste(colnames(df)[c(15:21)], "MW", sep = "_")
  
  # A. carbon i = Atom fraction excess of 15N in the labeled vs unlabeled treatment for taxon i (unitness)
  df.eaf <- as.data.frame(sapply(df.MW_lab, function(x) ((x - df.func$MW_unlabeled)/(df.func$MW_max - df.func$MW_unlabeled)) * (1 - 0.003663004)))
  # Change the column names of the EAF columns 
  colnames(df.eaf) <- paste(colnames(df)[c(15:21)], "EAF", sep = "_")
  
  # Column bind all the dataframes into one output 
  df.out <- cbind(df.func, df.diff, df.MW_lab, df.eaf)
  # Mean EAF 
  df.out$mean_eaf <- rowMeans(df.eaf, na.rm = TRUE)
  
  ## Each calculation are based on functions that should already be ran.
  # Get the number of treatment tubes each taxon in obs in 
  df.out$n_trt_tubes <- apply(df.eaf, 1, len.sumfun)
  # Standard deviation of EAF
  df.out$sd_eaf <- apply(df.eaf, 1, sd.sumfun)
  # Standard error of EAF
  df.out$se_eaf <- apply(df.eaf, 1, se.sumfun)
  # Calculate the 90% confidence interbal of EAF
  df.out$CI90_eaf <-apply(df.eaf, 1, CI.sumfun)
  
  df.out
}

# EAF calculation 
eaf.OF.C <- eaf_eqn(wad.OF.C)
eaf.OL.C <- eaf_eqn(wad.OL.C)
eaf.CF.C <- eaf_eqn(wad.CF.C)
eaf.CL.C <- eaf_eqn(wad.CL.C)


# Plot the EAF 

Plot.EAF.OF.C <- ggplot(eaf.OF.C, aes(x = reorder(taxon, mean_eaf, mean), y = mean_eaf, color = n_reps_unlab)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_eaf - se_eaf, ymax = mean_eaf + se_eaf)) + 
  geom_hline(yintercept = 0) + 
  theme_classic2() + 
  ggtitle("OF EAF")
Plot.EAF.OF.C

Plot.EAF.OL.C <- ggplot(eaf.OL.C, aes(x = reorder(taxon, mean_eaf, mean), y = mean_eaf, color = n_reps_unlab)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_eaf - se_eaf, ymax = mean_eaf + se_eaf)) + 
  geom_hline(yintercept = 0) + 
  theme_classic2() + 
  ggtitle("OL EAF")
Plot.EAF.OL.C

Plot.EAF.CF.C <- ggplot(eaf.CF.C, aes(x = reorder(taxon, mean_eaf, mean), y = mean_eaf, color = n_reps_unlab)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_eaf - se_eaf, ymax = mean_eaf + se_eaf)) + 
  geom_hline(yintercept = 0) + 
  theme_classic2() + 
  ggtitle("CF EAF")
Plot.EAF.CF.C

Plot.EAF.CL.C <- ggplot(eaf.CL.C, aes(x = reorder(taxon, mean_eaf, mean), y = mean_eaf, color = n_reps_unlab)) +
  geom_point() +
  geom_errorbar(aes(ymin = mean_eaf - se_eaf, ymax = mean_eaf + se_eaf)) + 
  geom_hline(yintercept = 0) + 
  theme_classic2() + 
  ggtitle("CL EAF")
Plot.EAF.CL.C


# Create new directory for APE dfs and plots 
dir.create("APE")

Plot.EAF.Organic <- ggarrange(Plot.EAF.OF.C, Plot.EAF.OL.C, nrow = 1, ncol = 2)
Plot.EAF.Organic 
Plot.EAF.Conventional <- ggarrange(Plot.EAF.CF.C, Plot.EAF.CL.C, nrow = 1, ncol = 2)
Plot.EAF.Conventional 

ggsave("APE/new.Plot.EAF.Organic.png", Plot.EAF.Organic)
ggsave("APE/new.Plot.EAF.Conventional.png", Plot.EAF.Conventional)


#### 12. Calculating the nitrogen assimilation ####
# We only need to include columns: _EAF, taxon id, and observation, and stats  

eaf1.OF.C <- select(eaf.OF.C, matches(c("taxon", "^xO5[0-9]F+_EAF", "mean_eaf", "n_trt_tubes", "sd_eaf", "se_eaf", "CI90_eaf" )))
eaf1.OL.C <- select(eaf.OL.C, matches(c("taxon", "^xO5[0-9]L+_EAF", "mean_eaf", "n_trt_tubes", "sd_eaf", "se_eaf", "CI90_eaf" )))
eaf1.CF.C <- select(eaf.CF.C, matches(c("taxon", "^xA5[0-9]F+_EAF", "mean_eaf", "n_trt_tubes", "sd_eaf", "se_eaf", "CI90_eaf" )))
eaf1.CL.C <- select(eaf.CL.C, matches(c("taxon", "^xA5[0-9]L+_EAF", "mean_eaf", "n_trt_tubes", "sd_eaf", "se_eaf", "CI90_eaf" )))


# Order by mean_EAF value (high to low)
eaf1.OF.C <- eaf1.OF.C[order(eaf1.OF.C$mean_eaf, decreasing = TRUE), ]
eaf1.OL.C <- eaf1.OL.C[order(eaf1.OL.C$mean_eaf, decreasing = TRUE), ]
eaf1.CF.C <- eaf1.CF.C[order(eaf1.CF.C$mean_eaf, decreasing = TRUE), ]
eaf1.CL.C <- eaf1.CL.C[order(eaf1.CL.C$mean_eaf, decreasing = TRUE), ]


# There is alot of NaN in my df. I am now going to remove it
#5/3 skip this for a sec and just put na.omit into functions
dim(eaf1.OF.C) #513 x 13
#all would have started with 513

#eaf1.OF.C <- na.omit(eaf1.OF.C)
#dim(eaf1.OF.C) #419 x 13 

#eaf1.OL.C <- na.omit(eaf1.OL.C)
#dim(eaf1.OL.C) #438 x 13 

#eaf1.CF.C <- na.omit(eaf1.CF.C)
#dim(eaf1.CF.C) #468 x 13 

#eaf1.CL.C <- na.omit(eaf1.CL.C)
#dim(eaf1.CL.C) #460 x 13 



 # Make eaf1 df into long format 

eaf1.long.OF.C <- pivot_longer(eaf1.OF.C, cols = names(eaf1.OF.C[2:8]), names_to = "tube", values_to = "EAF")
eaf1.long.OL.C <- pivot_longer(eaf1.OL.C, cols = names(eaf1.OL.C[2:8]), names_to = "tube", values_to = "EAF")
eaf1.long.CF.C <- pivot_longer(eaf1.CF.C, cols = names(eaf1.CF.C[2:8]), names_to = "tube", values_to = "EAF")
eaf1.long.CL.C <- pivot_longer(eaf1.CL.C, cols = names(eaf1.CL.C[2:8]), names_to = "tube", values_to = "EAF")


# Rename tube columns 
eaf1.long.OF.C$tube <-gsub("x", "", gsub("_EAF", "", eaf1.long.OF.C$tube))
eaf1.long.OL.C$tube <-gsub("x", "", gsub("_EAF", "", eaf1.long.OL.C$tube))
eaf1.long.CF.C$tube <-gsub("x", "", gsub("_EAF", "", eaf1.long.CF.C$tube))
eaf1.long.CL.C$tube <-gsub("x", "", gsub("_EAF", "", eaf1.long.CL.C$tube))

# Make an negative EAF values = NA
# Impossible to have a negative ort above 1 EAF: Either took up the labeled substrates or it didn't 
# Depending on how you want to analyze your data, you might want to leave the negative EAF
# March 21 2024, putting this back in but changing from 0 to NA and adding over 1 as NA


Ridge <- subset(median.all.asv.unfilt, Management == "Ridge")


eaf1.long.OF.C$EAF[eaf1.long.OF.C$EAF > 1] <- NA
eaf1.long.OL.C$EAF[eaf1.long.OL.C$EAF > 1] <- NA
eaf1.long.CF.C$EAF[eaf1.long.CF.C$EAF > 1] <- NA
eaf1.long.CL.C$EAF[eaf1.long.CL.C$EAF > 1] <- NA

eaf1.long.OF.C$EAF[eaf1.long.OF.C$EAF < 0] <- NA
eaf1.long.OL.C$EAF[eaf1.long.OL.C$EAF < 0] <- NA
eaf1.long.CF.C$EAF[eaf1.long.CF.C$EAF < 0] <- NA
eaf1.long.CL.C$EAF[eaf1.long.CL.C$EAF < 0] <- NA

# Bring in the ncopies.tube.taxa df and merge together with eaf df 

eaf.rel.OF.C <- merge(eaf1.long.OF.C, ncopies.tube.taxa, all.x = TRUE)
eaf.rel.OL.C <- merge(eaf1.long.OL.C, ncopies.tube.taxa, all.x = TRUE)
eaf.rel.CF.C <- merge(eaf1.long.CF.C, ncopies.tube.taxa, all.x = TRUE)
eaf.rel.CL.C <- merge(eaf1.long.CL.C, ncopies.tube.taxa, all.x = TRUE)

# double check to make sure that the correct tube is merge together 
unique(eaf.rel.OF.C$tube) #Yes, it worked 
unique(eaf.rel.CF.C$tube)

# Run function to calculate substrate assimilation
# There will be two columns for the uncorrected assimilation 
# and corrected assimilation. Because, relative abundance 
# data is use to calculate assimilation, the corrected 
# assimilation accounts for all the filtering that was done 

assimilation_eqn <- function (df){
  df1 <- df
  # Sum of filtered ncopies.taxa.sum (this value should be lower than the ncopies.tube.sum)
  for (i in 1:length(df1$tube)) {
    df1$rel.ncopies.sum.filtered[i] <- sum(df1$rel.ncopies[df1$tube %in% df1$tube[i]])
  }
  
  # piAi = Proportion of rel.ncopies (pi) * EAF (Ai)
  df1$piAi <- df1$rel.ncopies * df1$EAF
  
  # Sum of piAi by tube
  for(i in 1:length(df1$tube)){
    df1$piAi.sum[i] <- sum(df1$piAi[df1$tube %in% df1$tube[i]], na.rm = TRUE)
  }
  
  # Corrected piAi.sum -- Accounts for filtering (value less than 1)
  df1$piAi.sum.cor <- df1$piAi.sum/df1$rel.ncopies.sum.filtered
  
  # Calculate the percent substrate assimilation and corrected assimilation
  # Optional: You can change iso to match the isotope you used (13C, 15N, or 18O)
  df1$iso.Ai <- (df1$piAi/df1$piAi.sum)*100
  df1$iso.Ai.cor <- (df1$piAi/df1$piAi.sum.cor)*100
  
  # Double check to see if your calculation is correct
  # For the uncorrected iso.Ai, sum = 100
  # For the corrected iso.Ai, sum = rel.ncopies.sum.filtered * 100
  for(i in 1:length(df1$tube)){
    df1$iso.Ai.sum[i] <- sum(df1$iso.Ai[df1$tube %in% df1$tube[i]], na.rm = TRUE)
  }
  
  for(i in 1:length(df1$tube)){
    df1$iso.Ai.sum.corr[i] <- sum(df1$iso.Ai.cor[df1$tube %in% df1$tube[i]], na.rm = TRUE)
  }
  
  return (df1)
}


eaf.rel.OF.C <- assimilation_eqn(eaf.rel.OF.C)
eaf.rel.OL.C <- assimilation_eqn(eaf.rel.OL.C)
eaf.rel.CF.C <- assimilation_eqn(eaf.rel.CF.C)
eaf.rel.CL.C <- assimilation_eqn(eaf.rel.CL.C)


# Save 
write.csv(eaf.rel.OF.C, "APE/new.eaf.rel.OF.C.csv", row.names = FALSE)
write.csv(eaf.rel.OL.C, "APE/new.eaf.rel.OL.C.csv", row.names = FALSE)
write.csv(eaf.rel.CF.C, "APE/new.eaf.rel.CF.C.csv", row.names = FALSE)
write.csv(eaf.rel.CL.C, "APE/new.eaf.rel.CL.C.csv", row.names = FALSE)



Plotting
#####


#change relevant column names to match code
colnames(eaf.rel.OF.C)[which(names(eaf.rel.OF.C) == "iso.Ai.cor")] <- "NAi"
colnames(eaf.rel.OL.C)[which(names(eaf.rel.OL.C) == "iso.Ai.cor")] <- "NAi"
colnames(eaf.rel.CF.C)[which(names(eaf.rel.CF.C) == "iso.Ai.cor")] <- "NAi"
colnames(eaf.rel.CL.C)[which(names(eaf.rel.CL.C) == "iso.Ai.cor")] <- "NAi"


#merge farm files together
genus.all.conv <- rbind(eaf.rel.CF.C, eaf.rel.CL.C)
genus.all.org <- rbind(eaf.rel.OF.C, eaf.rel.OL.C)


length(unique(genus.all.conv$taxon)) #both 513... (Mar 21 2024; was 481 in code before...)
length(unique(genus.all.org$taxon)) #now 513 was 453 ... unclear if intentional or not

#make sure to keep taxon column for ID

genus.all.conv <- genus.all.conv %>% separate(taxon, 
                                              c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
                                              ".[a-z]__", remove = FALSE)

genus.all.org <- genus.all.org %>% separate(taxon, 
                                            c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
                                            ".[a-z]__", remove = FALSE)



#Cleaning Up
#remove .__ at end of taxa name on some 
genus.all.conv$Family <- gsub('.__', '', genus.all.conv$Family)
genus.all.conv$Order <- gsub('.__.__', '', genus.all.conv$Order)
genus.all.conv$Class <- gsub('.__.__.__', '', genus.all.conv$Class)
genus.all.conv$Domain <- gsub('d__', '', genus.all.conv$Domain)
genus.all.org$Family <- gsub('.__', '', genus.all.org$Family)
genus.all.org$Order <- gsub('.__.__', '', genus.all.org$Order)
genus.all.org$Class <- gsub('.__.__.__', '', genus.all.org$Class)
genus.all.org$Domain <- gsub('d__', '', genus.all.org$Domain)
genus.all.org$Domain <- gsub('.__.__.__.__.__', '', genus.all.org$Domain)
genus.all.conv$Domain <- gsub('.__.__.__.__.__', '', genus.all.conv$Domain)
genus.all.org$Phylum <- gsub('.__.__.__.__', '', genus.all.org$Phylum)
genus.all.conv$Phylum <- gsub('.__.__.__.__', '', genus.all.conv$Phylum)


#remove "chloroplasts and mitochondria"

genus.all.conv <- genus.all.conv %>%
  filter(!grepl('Chloroplast', Genus)) %>%
  filter(!grepl('Mitochondria', Genus))
  
genus.all.org <- genus.all.org %>%
  filter(!grepl('Chloroplast', Genus)) %>%
  filter(!grepl('Mitochondria', Genus))


#replace NA with highest given taxonomy. Have to start high and go low ****
genus.all.conv$Phylum[is.na(genus.all.conv$Phylum)] <- genus.all.conv$Domain[is.na(genus.all.conv$Phylum)]
genus.all.conv$Class[is.na(genus.all.conv$Class)] <- genus.all.conv$Phylum[is.na(genus.all.conv$Class)]
genus.all.conv$Order[is.na(genus.all.conv$Order)] <- genus.all.conv$Class[is.na(genus.all.conv$Order)]
genus.all.conv$Family[is.na(genus.all.conv$Family)] <- genus.all.conv$Order[is.na(genus.all.conv$Family)]

genus.all.org$Phylum[is.na(genus.all.org$Phylum)] <- genus.all.org$Domain[is.na(genus.all.org$Phylum)]
genus.all.org$Class[is.na(genus.all.org$Class)] <- genus.all.org$Phylum[is.na(genus.all.org$Class)]
genus.all.org$Order[is.na(genus.all.org$Order)] <- genus.all.org$Class[is.na(genus.all.org$Order)]
genus.all.org$Family[is.na(genus.all.org$Family)] <- genus.all.org$Order[is.na(genus.all.org$Family)]

#make NAs in Genus column "unknown"
genus.all.org$Genus[is.na(genus.all.org$Genus)] <- paste(genus.all.org$Family[is.na(genus.all.org$Genus)],c("Unknown"))
genus.all.conv$Genus[is.na(genus.all.conv$Genus)] <- paste(genus.all.conv$Family[is.na(genus.all.conv$Genus)],c("Unknown"))

#make "uncultured" = NA now
genus.all.conv[genus.all.conv == "uncultured" ] <- NA
genus.all.conv$Class[is.na(genus.all.conv$Class)] <- paste(genus.all.conv$Phylum[is.na(genus.all.conv$Class)])
genus.all.conv$Order[is.na(genus.all.conv$Order)] <- paste(genus.all.conv$Class[is.na(genus.all.conv$Order)])
genus.all.conv$Family[is.na(genus.all.conv$Family)] <- paste(genus.all.conv$Order[is.na(genus.all.conv$Family)])
genus.all.conv$Genus[is.na(genus.all.conv$Genus)] <- paste(genus.all.conv$Family[is.na(genus.all.conv$Genus)],c("uncultured"))

genus.all.org[genus.all.org == "uncultured" ] <- NA
genus.all.org$Class[is.na(genus.all.org$Class)] <- paste(genus.all.org$Phylum[is.na(genus.all.org$Class)])
genus.all.org$Order[is.na(genus.all.org$Order)] <- paste(genus.all.org$Class[is.na(genus.all.org$Order)])
genus.all.org$Family[is.na(genus.all.org$Family)] <- paste(genus.all.org$Order[is.na(genus.all.org$Family)])
genus.all.org$Genus[is.na(genus.all.org$Genus)] <- paste(genus.all.org$Family[is.na(genus.all.org$Genus)],c("uncultured"))

#double check to make sure didn't get rid of any ASVs except mitochondria and chloroplast
length(unique(genus.all.org$Genus)) 
#now 511 (Mar 21 2024; was 451)
length(unique(genus.all.conv$Genus))
#now 511 was 479

#want to remove NAs so 
median.na_rm <- function(x, ...){
  median=median(x, na.rm=TRUE, ...)
}

sum.na_rm <- function(x, ...){
  sum=sum(x, na.rm=TRUE, ...)
}

#I don't think matters for sum?
#sum for NAi
genus.all.conv <- summaryBy(NAi ~ taxon + incubation + tube, data = genus.all.conv, FUN = sum.na_rm, id = c("Genus", "Phylum", "Class", "Order", "Family", "EAF"))
genus.all.org <- summaryBy(NAi ~ taxon + incubation + tube, data = genus.all.org, FUN = sum.na_rm, id = c("Genus", "Phylum", "Class", "Order", "Family", "EAF"))

#median for EAF
genus.all.conv <- summaryBy(EAF ~ taxon + incubation + tube, data = genus.all.conv, FUN = median.na_rm, id = c("Genus", "Phylum", "Class", "Order", "Family", "NAi.sum.na_rm"))
genus.all.org <- summaryBy(EAF ~ taxon + incubation + tube, data = genus.all.org, FUN = median.na_rm, id = c("Genus", "Phylum", "Class", "Order", "Family", "NAi.sum.na_rm"))


#remove .sum.na_rm and .median.na_rm
names(genus.all.org) <- gsub('.sum.na_rm', '', names(genus.all.org))
names(genus.all.org) <- gsub('.median.na_rm', '', names(genus.all.org))
names(genus.all.conv) <- gsub('.sum.na_rm', '', names(genus.all.conv))
names(genus.all.conv) <- gsub('.median.na_rm', '', names(genus.all.conv))

write.csv(genus.all.conv, "genus.all.conv.csv")
write.csv(genus.all.org, "genus.all.org.csv")




#get "Lab_EAF", "Field_EAF", "Lab_NAi", and "Field_NAi" as columns 
#remove F and L from tube names?
genus.all.conv$tube <- gsub('F', '', genus.all.conv$tube)
genus.all.conv$tube <- gsub('L', '', genus.all.conv$tube)
genus.all.org$tube <- gsub('F', '', genus.all.org$tube)
genus.all.org$tube <- gsub('L', '', genus.all.org$tube)



conv.df <- pivot_wider(genus.all.conv, names_from = "incubation", values_from = c("NAi", "EAF"), id_cols = c("tube", "Genus", "Phylum", "Class", "Order", "Family", "taxon"))
org.df <- pivot_wider(genus.all.org, names_from = "incubation", values_from = c("NAi", "EAF"), id_cols = c("tube", "Genus", "Phylum", "Class", "Order", "Family", "taxon"))

write.csv(conv.df, "conv.df.csv")
write.csv(org.df, "org.df.csv")


#get both farms together

#add management column
conv.df <- conv.df %>%
  mutate(Management = "Valley") 
org.df <- org.df %>%
  mutate(Management = "Ridge") 


both.farms <- bind_rows(conv.df, org.df)
#7151 (511 *14 tubes)


write.csv(both.farms, "both.farms.csv")


median.all.asv.unfilt <- summaryBy(NAi_Field + NAi_Lab + EAF_Field + EAF_Lab ~ taxon + Management, data = both.farms, FUN = median.na_rm, id = c("Phylum", "Family", "Order", "Class", "Genus"), keep.names = TRUE)
#


#Note, some of the stuff immediately above this won't be used... edit before pub.

#Agreement Analyses
#####

#Fig 4a. is a Equivariant Passing-Bablock Regression over a scatter plot
#Passing Bablok regression like Deming, but is a non-parametric method that does not make any assumptions about the distributions of the samples or their measurement errors (Passing and Bablok, 1983). 
#The method does, however, assume the two variables are highly correlated and have a linear relationship.
#For Passing-Bablok regression a nested bootstrap (bias corrected and accelerated) interval is preferred, but this can take a very long time with large datasets, where an approximate Passing-Bablok procedure may be more practical.
#If the confidence interval for the slope does not contain the value 1 then you reject the null hypothesis that the slope is equal to 1 - in other words there is statistically significant evidence of at least a proportional difference between the two methods.
#If the confidence interval for the intercept does not contain the value 0 then you reject the null hypothesis that the intercept is equal to 0 - in other words there is statistically significant evidence that the methods differ by at least a constant amount.

#scatter plot with agreement line
# You can select a residual plot, which is useful for spotting outliers and non linear patterns, and for checking how the agreement varies over the range of measurement.


Ridge <- subset(median.all.asv.unfilt, Management == "Ridge")
Valley <- subset(median.all.asv.unfilt, Management == "Valley")


NAi.PaBa.Ridge <- mcreg(x = Ridge$NAi_Lab,
                        y = Ridge$NAi_Field,
                        method.reg = "PBequi",
                        mref.name = "Lab",
                        mtest.name = "Field",
                        sample.names = NULL,
                        na.rm = TRUE
) #default alpha = 0.05


EAF.PaBa.Ridge <- mcreg(x = Ridge$EAF_Lab,
                        y = Ridge$EAF_Field,
                        method.reg = "PBequi",
                        mref.name = "Lab",
                        mtest.name = "Field",
                        sample.names = NULL,
                        na.rm = TRUE
) #default alpha = 0.05


print(NAi.PaBa.Ridge@para)

#               EST SE           LCI         UCI
#Intercept 0.0000000 NA -7.466146e-06 0.000114312
#Slope     0.9513317 NA  9.040492e-01 0.995927631

#CI for Intercept does contain the value 0, accept null hypothesis, NO evidence that methods differ by at least a constant amount
#CI for Slope DOES NOT contain the value 1, so statistically significant evidence of even proportional differences between the two methods.


print(EAF.PaBa.Ridge@para)

#Please note: 
 # 4 of 511 observations contain missing values and have been removed.
# Number of data points in analysis is 507.

#                 EST SE         LCI         UCI
#Intercept -0.01488444 NA -0.03059966 -0.00205814
#Slope      0.86630679 NA  0.78776927  0.94927549

#CI for Intercept does NOT contain the value 0, reject null hypothesis, evidence that methods differ by at least a constant amount
#CI for Slope DOES not contain the value 1, so statistically significant evidence of a proportional differences between the two methods.

NAi.PaBa.Valley <- mcreg(x = Valley$NAi_Lab,
                         y = Valley$NAi_Field,
                         method.reg = "PBequi",
                         mref.name = "Lab",
                         mtest.name = "Field",
                         sample.names = NULL,
                         na.rm = TRUE
) #default alpha = 0.05

#252 of 511 observations contain missing values and have been removed.
#Number of data points in analysis is 259.


EAF.PaBa.Valley <- mcreg(x = Valley$EAF_Lab,
                         y = Valley$EAF_Field,
                         method.reg = "PBequi",
                         mref.name = "Lab",
                         mtest.name = "Field",
                         sample.names = NULL,
                         na.rm = TRUE
                         ) #default alpha = 0.05

#Please note: 
#5 of 511 observations contain missing values and have been removed.
#Number of data points in analysis is 506.

print(NAi.PaBa.Valley@para)

#                   EST SE          LCI      UCI
#Intercept -2.229398e-05 NA -0.001274343 0.000000
#Slope      1.018069e+00 NA  0.978567933 1.062565

#CI for Intercept does contain the value 0, accept null hypothesis, no evidence that methods differ by at least a constant amount
#CI for Slope DOES contain the value 1, so no statistically significant evidence of even proportional differences between the two methods.


print(EAF.PaBa.Valley@para)
#                  EST SE         LCI        UCI
#Intercept 0.01535646 NA 0.003636968 0.02618975
#Slope     0.87354369 NA 0.793714878 0.95367251
#CI for Intercept does NOT contain the value 0, reject null hypothesis, evidence that methods differ by at least a constant amount
#CI for Slope DOES not contain the value 1, so statistically significant evidence of a proportional differences between the two methods.
#^^^ intercept conclusion is different than Ridge



colors <- sample(c("red", "blue"), replace = T, size = length(median.all.asv.unfilt))

MCResult.plot(NAi.PaBa.Ridge,
              equal.axis = FALSE,
              x.lab = "Lab",
              y.lab = "Field",
              alpha = 0.05,
              points.col = "black",
              points.pch = 21,
              ci.area = TRUE,
              ci.area.col = "#0000FF50",
              main = expression('%'^~'15'~N~Assimilated),
              sub = "",
              add.grid = FALSE,
              points.cex = 1,
              legend.place = "topleft",
              cor.method = "pearson",
              digits = list(coef = 4, cor = 3)
)


printSummary(NAi.PaBa.Ridge)


MCResult.plot(EAF.PaBa.Ridge,
              equal.axis = FALSE,
              x.lab = "Lab",
              y.lab = "Field",
              alpha = 0.05,
              points.col = "black",
              points.pch = 21,
              ci.area = TRUE,
              ci.area.col = "#0000FF50",
              main =  expression(Relative^~'15'~N~Assimilation~Rate~'/'~5~Days),
              sub = "",
              add.grid = FALSE,
              points.cex = 1,
              legend.place = "topleft",
              cor.method = "spearman",
              digits = list(coef = 4, cor = 3)
)

printSummary(EAF.PaBa.Ridge)

MCResult.plot(EAF.PaBa.Valley,
              equal.axis = FALSE,
              x.lab = "Lab",
              y.lab = "Field",
              alpha = 0.05,
              points.col = "black",
              points.pch = 21,
              ci.area = TRUE,
              ci.area.col = "#0000FF50",
              main =  expression(Relative^~'15'~N~Assimilation~Rate~'/'~5~Days),
              sub = "",
              add.grid = FALSE,
              points.cex = 1,
              legend.place = "topleft",
              cor.method = "spearman",
              digits = list(coef = 4, cor = 3)
)

printSummary(EAF.PaBa.Valley)


MCResult.plot(NAi.PaBa.Valley,
              equal.axis = FALSE,
              x.lab = "Lab",
              y.lab = "Field",
              alpha = 0.05,
              points.col = "black",
              points.pch = 21,
              ci.area = TRUE,
              ci.area.col = "#0000FF50",
              main = expression('%'^~'15'~N~Assimilated),
              sub = "",
              add.grid = FALSE,
              points.cex = 1,
              legend.place = "topleft",
              cor.method = "spearman",
              digits = list(coef = 4, cor = 3)
)

printSummary(NAi.PaBa.Valley)





#reg line for both together

NAi.PaBa.unfilt <- mcreg(x = median.all.asv.unfilt$NAi_Lab,
                         y = median.all.asv.unfilt$NAi_Field,
                         method.reg = "PBequi",
                         mref.name = "Lab",
                         mtest.name = "Field",
                         sample.names = NULL,
                         na.rm = TRUE
) #default alpha = 0.05



print(NAi.PaBa.unfilt@para)

#                EST SE          LCI      UCI
#Intercept 0.0000000 NA -0.000181743 0.000000
#Slope     0.9854867 NA  0.952124175 1.019613

#Intercept does include 0 so no constant diff
#slope does include 1 so no prop diff


EAF.PaBa.unfilt <- mcreg(x = median.all.asv.unfilt$EAF_Lab,
                         y = median.all.asv.unfilt$EAF_Field,
                         method.reg = "PBequi",
                         mref.name = "Lab",
                         mtest.name = "Field",
                         sample.names = NULL,
                         na.rm = TRUE
) #default alpha = 0.05

#9 of 1022 observations contain missing values and have been removed.
#Number of data points in analysis is 1013.

print(EAF.PaBa.unfilt@para)

#                  EST SE        LCI        UCI
#Intercept 0.01038293 NA 0.00109323 0.01882823
#Slope     0.81225241 NA 0.75517484 0.87095758

#Intercept CI does not include 0, so constant diff
#Slope does not include 1, so proportional diff


#possible supplemental
compareFit(NAi.PaBa.Ridge, NAi.PaBa.Valley,NAi.PaBa.unfilt, EAF.PaBa.Ridge, EAF.PaBa.Valley, EAF.PaBa.unfilt)



#FOR PUBLICATION (i.e. WITH COLOR)
plot(NAi.PaBa.Valley, add.legend = FALSE, 
     points.col = rgb(red = 147, green = 112, blue = 219,maxColorValue = 255), points.pch = 21,
     reg.col = rgb(red = 147, green = 112, blue = 219,  maxColorValue = 255),
     ci.area = TRUE, ci.area.col = rgb(red = 147, green = 112, blue = 219, alpha = 150,  maxColorValue = 255),
     identity = TRUE, identity.col = "black", identity.lwd = 2, 
     add.grid = FALSE,
     main = expression('%'~N~Assimilated),
     add.cor = FALSE)

plot(NAi.PaBa.Ridge,add.legend = FALSE, add = TRUE,
     points.col = rgb(red = 60, green = 179, blue = 113, maxColorValue = 255), points.pch = 21,
     reg.col =  rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),
     main = expression('%'~N~Assimilated),
     ci.area = TRUE, ci.area.col =  rgb(red = 60, green = 179, blue = 113, alpha = 150, maxColorValue = 255),
     identity = FALSE, add.cor = FALSE, 
     add.grid = FALSE)


plot(NAi.PaBa.unfilt,add.legend = FALSE, add = TRUE,
     draw.points = FALSE, 
     reg.col =  rgb(red = 255, green = 255, blue = 51, maxColorValue = 255),
     ci.area = TRUE, ci.area.col =  rgb(red = 255, green = 255, blue = 51, alpha = 150, maxColorValue = 255),
     identity = FALSE,
     add.grid = FALSE)

includeLegend(place="topleft",models=list(NAi.PaBa.Ridge,NAi.PaBa.Valley,NAi.PaBa.unfilt ),
              colors=c(rgb(red = 147, green = 112, blue = 219,maxColorValue = 255),rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),rgb(red = 255, green = 255, blue = 51, maxColorValue = 255)),
              design="1", digits=4, cex = 1.2,
              model.names = c("Valley", "Ridge", "All Sites"))

NAi_PaBa.plot <- recordPlot()

#saved
dev.off()

plot(EAF.PaBa.Valley, add.legend = FALSE, 
     points.col = rgb(red = 147, green = 112, blue = 219,maxColorValue = 255), points.pch = 21,
     reg.col = rgb(red = 147, green = 112, blue = 219,  maxColorValue = 255),
     ci.area = TRUE, ci.area.col = rgb(red = 147, green = 112, blue = 219, alpha = 150,  maxColorValue = 255),    identity = TRUE, identity.col = "black", identity.lwd = 2, 
     add.grid = FALSE,
     main = expression(Relative~N~Assimilation~Rate~'/'~5~Days),
     add.cor = FALSE)

plot(EAF.PaBa.Ridge,add.legend = FALSE, add = TRUE,
     points.col = rgb(red = 60, green = 179, blue = 113, maxColorValue = 255), points.pch = 21,
     reg.col =  rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),
     ci.area = TRUE, ci.area.col =  rgb(red = 60, green = 179, blue = 113, alpha = 150, maxColorValue = 255),
     main = expression(Relative~N~Assimilation~Rate~'/'~5~Days),
     identity = FALSE, add.cor = FALSE, cor.method = 'pearson',
     add.grid = FALSE)

plot(EAF.PaBa.unfilt,add.legend = FALSE, add = TRUE,
     draw.points = FALSE, 
     reg.col =  rgb(red = 255, green = 255, blue = 51, maxColorValue = 255),
     ci.area = TRUE, ci.area.col =  rgb(red = 255, green = 255, blue = 51, alpha = 150, maxColorValue = 255),
     identity = FALSE,
     add.cor = TRUE,
     add.grid = FALSE)

includeLegend(place="topleft",models=list(EAF.PaBa.Valley,EAF.PaBa.Ridge,EAF.PaBa.unfilt ),
              colors=c(rgb(red = 147, green = 112, blue = 219,maxColorValue = 255),rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),rgb(red = 255, green = 255, blue = 51, maxColorValue = 255)),
              design="1", digits=4, cex = 1.2,
              model.names = c("Valley", "Ridge", "All Sites"))

EAF_PaBa.plot <- recordPlot()


#for combining
dev.off()
plot(NAi.PaBa.Valley, add.legend = FALSE, 
     points.col = rgb(red = 147, green = 112, blue = 219,maxColorValue = 255), points.pch = 21,
     reg.col = rgb(red = 147, green = 112, blue = 219,  maxColorValue = 255),
     ci.area = TRUE, ci.area.col = rgb(red = 147, green = 112, blue = 219, alpha = 150,  maxColorValue = 255),
     identity = TRUE, identity.col = "black", identity.lwd = 2, 
     add.grid = FALSE,
     main = "", sub = "",
     add.cor = FALSE)

plot(NAi.PaBa.Ridge, add.legend = FALSE, add = TRUE,
     points.col = rgb(red = 60, green = 179, blue = 113, maxColorValue = 255), points.pch = 21,
     reg.col =  rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),
     main = expression('%'~N~Assimilated),
     ci.area = TRUE, ci.area.col =  rgb(red = 60, green = 179, blue = 113, alpha = 150, maxColorValue = 255),
     identity = FALSE, add.cor = FALSE,
     add.grid = FALSE)


plot(NAi.PaBa.unfilt,add.legend = FALSE, add = TRUE,
     draw.points = FALSE, 
     reg.col =  rgb(red = 255, green = 255, blue = 51, maxColorValue = 255),
     ci.area = TRUE, ci.area.col =  rgb(red = 255, green = 255, blue = 51, alpha = 150, maxColorValue = 255),
     identity = FALSE,
     add.grid = FALSE)

includeLegend(place="topleft",models=list(NAi.PaBa.Ridge,NAi.PaBa.Valley,NAi.PaBa.unfilt ),
              colors=c(rgb(red = 147, green = 112, blue = 219,maxColorValue = 255),rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),rgb(red = 255, green = 255, blue = 51, maxColorValue = 255)),
              design="1", digits=4, cex = 0.9,
              model.names = c("Valley", "Ridge", "All Sites"))

NAi_PaBa.plot.2 <- recordPlot()
#saved as NAi_PaBa_Nocapt


dev.off()

plot(EAF.PaBa.Valley, add.legend = FALSE, 
     points.col = rgb(red = 147, green = 112, blue = 219,maxColorValue = 255), points.pch = 21,
     reg.col = rgb(red = 147, green = 112, blue = 219,  maxColorValue = 255),
     ci.area = TRUE, ci.area.col = rgb(red = 147, green = 112, blue = 219, alpha = 150,  maxColorValue = 255),    identity = TRUE, identity.col = "black", identity.lwd = 2, 
     add.grid = FALSE,
     main = "", sub = "",
     add.cor = FALSE)

plot(EAF.PaBa.Ridge,add.legend = FALSE, add = TRUE,
     points.col = rgb(red = 60, green = 179, blue = 113, maxColorValue = 255), points.pch = 21,
     reg.col =  rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),
     ci.area = TRUE, ci.area.col =  rgb(red = 60, green = 179, blue = 113, alpha = 150, maxColorValue = 255),
     main = expression(Relative~N~Assimilation~Rate~'/'~5~Days),
     identity = FALSE, add.cor = FALSE, cor.method = 'pearson',
     add.grid = FALSE)

plot(EAF.PaBa.unfilt,add.legend = FALSE, add = TRUE,
     draw.points = FALSE, 
     reg.col =  rgb(red = 255, green = 255, blue = 51, maxColorValue = 255),
     ci.area = TRUE, ci.area.col =  rgb(red = 255, green = 255, blue = 51, alpha = 150, maxColorValue = 255),
     identity = FALSE,
     add.cor = TRUE,
     add.grid = FALSE)

includeLegend(place="topleft",models=list(EAF.PaBa.Ridge,EAF.PaBa.Valley,EAF.PaBa.unfilt ),
              colors=c(rgb(red = 147, green = 112, blue = 219,maxColorValue = 255),rgb(red = 60, green = 179, blue = 113, maxColorValue = 255),rgb(red = 255, green = 255, blue = 51, maxColorValue = 255)),
              design="1", digits=4, cex = 0.9,
              model.names = c("Valley", "Ridge", "All Sites"))

EAF_PaBa.plot2 <- recordPlot()
#saved as EAF_PaBa_NoCapt

#These will be Fig 2 C and D 


#Density Plot
##### 
#Fig 2 Density Plot



NAi.long <- pivot_longer(both.farms,
                         cols = c("NAi_Lab", "NAi_Field"),
                         names_to = "Incubation",
                         names_prefix = "NAi_",
                         values_to = "NAi")

EAF.only <- pivot_longer(both.farms,
                         cols = c("EAF_Lab", "EAF_Field"),
                         names_to = "Incubation",
                         names_prefix = "EAF_",
                         values_to = "EAF")

#remove extra (I.e. either EAF or NAi) columns. CHECK THIS IS RIGHT EACH TIME BECAUSE ALWAYS GETS THROWN OFF
NAi.long <- NAi.long[,-c(8:9)]
EAF.only <- EAF.only[,-c(8:9)]

#now merge EAF and NAi
bf.bi.long <- left_join(NAi.long, EAF.only, 
                        by = c("tube", "Genus", "Class", "Phylum", "Family", "Order", "Management", "Incubation", "taxon"),
                        keep = FALSE)

bf.bi.long$Incubation <- as.factor(bf.bi.long$Incubation)
bf.bi.long$Management <- as.factor(bf.bi.long$Management)

#Note Incubation = Method
bf.bi.long$Management_Method <- paste(bf.bi.long$Management,bf.bi.long$Incubation)
bf.bi.long$Management_Method <- as.factor(bf.bi.long$Management_Method)

write.csv(bf.bi.long, "bf.bi.long.csv")


NAicaption = str_wrap("")


ggsave("PubGraphs/NAi_density.png")
ggsave("PubGraphs/NAi_density.svg")


bf.bi.long <- read.csv("bf.bi.long.csv")

#use median values

median.bf.bi.long <- bf.bi.long %>%
  #make NAi 0's = NA
  mutate_at(c('NAi'), ~na_if(., 0)) %>%
  group_by(taxon, Genus, Phylum, Order, Family, Class, Management_Method) %>%
  summarise(medianNAi = median.na_rm(NAi), medianEAF = median.na_rm(EAF))
  
write.csv(median.bf.bi.long, "median.bf.bi.long.csv")

EAF_dens_median <- median.bf.bi.long %>%
  ggplot( aes(x=medianEAF, color=Management_Method)) +
  geom_density(size = 1.3, alpha= 0.7 , na.rm = TRUE) +
  scale_color_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  labs(title = expression(Relative~N~Assimilation~Rate), y = 'Density', x = '(per 5 Days)') +
  theme_classic() +
  labs(color = "Site - Method") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12),
      #  axis.title.x = element_blank(),
        axis.text = element_text(size = 12),
        plot.caption = element_text(hjust = 0),
        plot.title = element_text(face = 'bold', hjust = 0.5, size = 14),
        axis.ticks = element_blank(),
        plot.background = element_blank()) 

EAF_dens_median
ggsave("PubGraphs/EAF_Fig2A.png")
ggsave("PubGraphs/EAF_Fig2A.svg")


NAi_dens.noyaxislab_median <-
  median.bf.bi.long %>%
  ggplot(aes(x=medianNAi, color=Management_Method)) +
  geom_density(size = 1.3, alpha= 0.7, na.rm = TRUE ) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) +
 # scale_color_viridis(discrete=TRUE, option = "H") +
  scale_color_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  labs(title = expression('%'~N~Assimilated), y = 'Density', x = 'log(10)') +
  theme_classic() +
  labs(color = "Site - Method") +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 12),
       # axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        plot.caption = element_text(hjust = 0),
       plot.title = element_text(face = 'bold', hjust = 0.5, size = 14),
        axis.ticks = element_blank(),
        plot.background = element_blank()) 



NAi_dens.noyaxislab_median

ggsave("PubGraphs/NAi_density_Fig2B.png")
ggsave("PubGraphs/NAi_density_Fig2B.svg")

ggarrange(EAF_dens_median,NAi_dens.noyaxislab_median, ncol = 2, common.legend = TRUE, legend = "bottom", labels = c("A", "B"), font.label = list(size = 16))

#Used for Fig 2 A and B
ggsave("PubGraphs/Combined_density_median.png")
ggsave("PubGraphs/Combined_density.median.svg")

  
#####
#Ridgeline Plot Fig 3


tokeep <- bf.bi.long %>%
  group_by(taxon, Management_Method) %>%
  summarise(NAi.n = n_distinct(NAi, na.rm = FALSE)) %>%
  filter(NAi.n > 3) %>% #change to <2 and top 3 lines to see how many and which will be removed... so 41 total
  #at least 3 replicates req
  summarise(trt.n = n_distinct(Management_Method)) %>%
  filter(trt.n == 4)

#these are ones to keep

bf.bi.long.filt <- merge(bf.bi.long, tokeep, by.y = 'taxon')

#make sure pkg 'plyr' was loaded before 'dplyr' or else won't group 
bf.bi.long.filt$Genus <- as_factor(bf.bi.long.filt$Genus)


NAi.t20.LabIncs<- bf.bi.long.filt %>%  
  filter(Incubation == "Lab") %>%
  group_by(taxon) %>%
  summarise(median = median(NAi, na.rm = TRUE)) %>% #?na.rm = FALSE otherwise deletes whole taxon for 1 rep NA
  slice_max(n = 20, order_by = median).             #now not sure true because deletes bacillus?

NAi.t20.FieldIncs <- bf.bi.long.filt %>%
  filter(Incubation == "Field") %>%
  group_by(taxon) %>%
  summarise(median = median(NAi, na.rm = TRUE)) %>%
  slice_max(n = 20, order_by = median)



combined.NAi <- full_join(NAi.t20.FieldIncs, NAi.t20.LabIncs, by = "taxon")


taxon <- c("taxon")

NAi.names <- combined.NAi[, taxon]


genus.NAi.t20 <- merge(bf.bi.long, NAi.names, by.y = taxon)

length(unique(genus.NAi.t20$Genus)) #21 #23


#Get actual highest median values for discussion 
winnersNAi <- genus.NAi.t20 %>%
  group_by(taxon, Genus) %>%
  summarise(median.EAF = median.na_rm(EAF), median.NAi = median.na_rm(NAi))
# is Candidatus_Udaeobacter EAF0.1490902 NAi 1.4272290
#after that close but Vicinamibacterales uncultured #0.1335499 #3.8814699

#5/14 switched to one-way ANOVA.... since sites aren't comparable .  one site at a time
aov_Valley <- genus.NAi.t20 %>%
  filter(Management == "Valley") %>%
  group_by(Genus) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(NAi ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay

aov_Valley$aov_results[[2]]#change ref for each taxa

#copy and pasted results to .doI c titled "NAi_genus_ridgeline_stats"

#NEWDATA <- data.frame(matrix(unlist(DATA), nrow=length(DATA), byrow=T))

unlist.NAi.Valley <- data.frame(matrix(unlist(aov_Valley$aov_results), nrow = length(aov_Valley$aov_results), byrow = T))
write.csv(unlist.NAi.Valley, "NAi.valley.1wayanova.csv")


aov_Ridge <- genus.NAi.t20 %>%
  filter(Management == "Ridge") %>%
  group_by(Genus) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(NAi ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay

aov_Ridge$aov_results[[1]]#change ref for each taxa


#NEWDATA <- data.frame(matrix(unlist(DATA), nrow=length(DATA), byrow=T))

unlist.NAi.Ridge <- data.frame(matrix(unlist(aov_Ridge$aov_results), nrow = length(aov_Ridge$aov_results), byrow = T))
write.csv(unlist.NAi.Ridge, "NAi.ridge.1wayanova.csv")


#see google sheet for 1 way aov results



genus.NAi.t20$Genus <- as.factor(genus.NAi.t20$Genus)

#library(ggridges)
genus.ridgeplot <- genus.NAi.t20 %>%
  mutate(Genus = fct_reorder(Genus, NAi, .fun = 'median', .na_rm = TRUE)) %>%
  ggplot(aes(x = NAi, y = Genus, fill = Management_Method, color = Management_Method)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.0001, alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  #scale_color_viridis(discrete=TRUE) +
  #scale_fill_viridis(discrete=TRUE) +
  scale_color_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  scale_fill_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
 scale_x_continuous(limits = c(0, 10), expand = c(0,0)) +  #Bacillus has one point by 20!! removing for graph but note in capt
  labs(title = expression('%'~N~Assimilated), y = '', x = '%', caption = "") + #y axis "Genus" removed for ggarrange
  theme_clean() +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 16), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.5, 'cm'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12 ),
        axis.text.x = element_text(size = 16 ) )  
genus.ridgeplot

ggsave("PubGraphs/NAi_ridgeline_alt_upd.png") #purple and green saved as alt
ggsave("PubGraphs/NAi_ridgeline_alt_upd.svg")
#add sig levels later in illustrator

#Now EAF

EAF.t20.LabIncs <- bf.bi.long.filt %>%
  filter(Incubation == "Lab") %>%
  group_by(taxon) %>%
  summarise(median = median(EAF, na.rm =  FALSE)) %>%  #note that na.rm = TRUE here completely changes things. Keeping False because all data with NAs has at least 3 reps.
  slice_max(n = 20, order_by = median)

EAF.t20.FieldIncs <- bf.bi.long.filt %>%
  filter(Incubation == "Field") %>%
  group_by(taxon) %>%
  summarise(median = median(EAF, na.rm = FALSE)) %>%
              slice_max(n = 20, order_by = median)



combined.EAF <- full_join(EAF.t20.FieldIncs, EAF.t20.LabIncs, by = "taxon")


EAF.names <- combined.EAF[, taxon]


genus.EAF.t20 <- merge(bf.bi.long, EAF.names, by.y = taxon)
length(unique(genus.EAF.t20$Genus)) #37 /50

#top EAF assim. genus with highest median NAi for discussion/results
winnersEAF <- genus.EAF.t20 %>%
  group_by(taxon, Genus) %>%
  summarise(median.EAF = median.na_rm(EAF), median.NAi = median.na_rm(NAi))
# is Nitrososphaeraceae Unknown 0.276558796 Nai, 0.2231532 EAF 
#others with NAi over 0.2 = Terrimonas EAF 0.149784 NAi 0.221999107
#Luteolibacter EAF 0.2338113 NAi = 0.294410017


#5/14 switch to one-way.... since sites aren't comparable .  one site at a time
aov_ValleyEAF <- genus.EAF.t20 %>%
  filter(Management == "Valley") %>%
  group_by(Genus) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(EAF ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay


#copy and pasted results to .doc titled "EAF_genus_ridgeline_stats"

#NEWDATA <- data.frame(matrix(unlist(DATA), nrow=length(DATA), byrow=T))

unlist.EAF.Valley <- data.frame(matrix(unlist(aov_ValleyEAF$aov_results), nrow = length(aov_ValleyEAF$aov_results), byrow = T))
write.csv(unlist.EAF.Valley, "EAF.valley.1wayanova.csv")


aov_RidgeEAF <- genus.EAF.t20 %>%
  filter(Management == "Ridge") %>%
  group_by(Genus) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(EAF ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay

aov_RidgeEAF$aov_results[[1]]#change ref for each taxa


#NEWDATA <- data.frame(matrix(unlist(DATA), nrow=length(DATA), byrow=T))

unlist.EAF.Ridge <- data.frame(matrix(unlist(aov_RidgeEAF$aov_results), nrow = length(aov_RidgeEAF$aov_results), byrow = T))
write.csv(unlist.EAF.Ridge, "EAF.ridge.1wayanova.csv")

#2 way anova for discussion of some

aov_2_EAF <- genus.EAF.t20 %>%
  group_by(Genus) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(EAF ~ Incubation * Management, data = .), na.rm = FALSE)))


#see google sheet for 1 way aov results


genus.EAF.ridgeplot <- genus.EAF.t20 %>%
  mutate(Genus = fct_reorder(Genus, EAF, .fun = 'median', .na_rm = TRUE)) %>%  #na_rm important here for ordering.. This median is across BOTH farms.
  ggplot(aes(x = EAF, y = Genus, fill = Management_Method, color = Management_Method)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.0001, alpha = 0.4) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
 # scale_color_viridis(discrete=TRUE) +
  #scale_fill_viridis(discrete=TRUE) +
  scale_color_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  scale_fill_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  labs(title = expression(Relative~N~Assimilation~Rate), x = 'per 5 days', y = 'Genus', caption = "") +
  theme_clean() +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold", size = 12),
        plot.caption = element_text(hjust = 0),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 16), 
        legend.title = element_blank(), 
        legend.key.size = unit(0.5, 'cm'),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 12 ),
        axis.text.x = element_text(size = 16 )) 
genus.EAF.ridgeplot

ggsave("PubGraphs/20EAF_ridgeline_alt.png")
ggsave("PubGraphs/20EAF_ridgeline_alt.svg")


ggarrange( genus.EAF.ridgeplot, genus.ridgeplot, labels = c("A", "B"), ncol = 2, nrow =1, common.legend = TRUE, legend = "bottom")

ggsave("PubGraphs/Ridgeline_Median_2.png")
ggsave("PubGraphs/Ridgeline_Median_2.svg")


#### Phylum Level Supp Fig 3

#sum NAi and median EAF
#very important 

bf.bi.long_Phylum <- bf.bi.long.filt %>%
  group_by(tube, Phylum, Management_Method, Management, Incubation) %>%
  summarise(summed.NAi = sum.na_rm(NAi), median.EAF = median.na_rm(EAF))
#this func is na.rm = TRUE

bf.bi.long_Phylum$Phylum <- as_factor(bf.bi.long_Phylum$Phylum)

#top median for field vs top median for lab?

#NAi start

NAi.t20.LabIncs_Phylum <- bf.bi.long_Phylum %>%
  filter(Incubation == "Lab") %>%
  group_by(Phylum) %>%
  summarise(median = median(summed.NAi, na.rm = TRUE)) %>%
  slice_max(n = 20, order_by = median)

NAi.t20.FieldIncs_Phylum <- bf.bi.long_Phylum %>%
  filter(Incubation == "Field") %>%
  group_by(Phylum) %>%
  summarise(median = median(summed.NAi, na.rm = TRUE)) %>%
  slice_max(n = 20, order_by = median)



combined.NAi.Phylum <- full_join(NAi.t20.FieldIncs_Phylum, NAi.t20.LabIncs_Phylum, by = "Phylum")


#only need Phylum col remove vals

Phylum <- c("Phylum")

NAi.names_Phylum<- combined.NAi.Phylum[, Phylum]


NAi.t20_Phylum <- merge(bf.bi.long_Phylum, NAi.names_Phylum, by.y = Phylum)

length(unique(NAi.t20_Phylum$Phylum)) #20



aov_ValleyNAi_Phy <- NAi.t20_Phylum %>%
  filter(Management == "Valley") %>%
  group_by(Phylum) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(summed.NAi ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay


#NEWDATA <- data.frame(matrix(unlist(DATA), nrow=length(DATA), byrow=T))

unlist.NAi.Valley.Phy <- data.frame(matrix(unlist(aov_ValleyNAi_Phy$aov_results), nrow = length(aov_ValleyNAi_Phy$aov_results), byrow = T))
write.csv(unlist.NAi.Valley.Phy, "NAi.Valley.Phy.1wayanova.csv")


aov_RidgeNAi_Phy <- NAi.t20_Phylum %>%
  filter(Management == "Ridge") %>%
  group_by(Phylum) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(summed.NAi ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay



unlist.NAi.Ridge.Phy <- data.frame(matrix(unlist(aov_RidgeNAi_Phy$aov_results), nrow = length(aov_RidgeNAi_Phy$aov_results), byrow = T))
write.csv(unlist.NAi.Ridge.Phy, "NAi.Ridge.Phy.1wayanova.csv")
#none sig
#see google sheet for 1 way aov results

#Now EAF
EAF.t20.LabIncs_Phylum <- bf.bi.long_Phylum %>%
  filter(Incubation == "Lab") %>%
  group_by(Phylum) %>%
  summarise(median = median(median.EAF, na.rm = TRUE)) %>%
  slice_max(n = 20, order_by = median)

EAF.t20.FieldIncs_Phylum <- bf.bi.long_Phylum %>%
  filter(Incubation == "Field") %>%
  group_by(Phylum) %>%
  summarise(median = median(median.EAF, na.rm = TRUE)) %>%
  slice_max(n = 20, order_by = median)



combined.EAF.Phylum <- full_join(EAF.t20.FieldIncs_Phylum, EAF.t20.LabIncs_Phylum, by = "Phylum")


#only need Phylum col remove vals

Phylum <- c("Phylum")

EAF.names_Phylum<- combined.EAF.Phylum[, Phylum]


EAF.t20_Phylum <- merge(bf.bi.long_Phylum, EAF.names_Phylum, by.y = Phylum)

length(unique(EAF.t20_Phylum$Phylum)) #24



aov_ValleyEAF_Phy <- EAF.t20_Phylum %>%
  filter(Management == "Valley") %>%
  group_by(Phylum) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(median.EAF ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay


#NEWDATA <- data.frame(matrix(unlist(DATA), nrow=length(DATA), byrow=T))

unlist.EAF.Valley.Phy <- data.frame(matrix(unlist(aov_ValleyEAF_Phy$aov_results), nrow = length(aov_ValleyEAF_Phy$aov_results), byrow = T))
write.csv(unlist.EAF.Valley.Phy, "EAF.Valley.Phy.1wayanova.csv")


aov_RidgeEAF_Phy <- EAF.t20_Phylum %>%
  filter(Management == "Ridge") %>%
  group_by(Phylum) %>%
  nest() %>%
  mutate(.data = .,
         aov_results = data %>% purrr::map(.x = ., .f = ~ summary(aov(median.EAF ~ Incubation, data = .), na.rm = FALSE)))

#hopefully NA removal = FALSE here is okay


#NEWDATA <- data.frame(matrix(unlist(DATA), nrow=length(DATA), byrow=T))

unlist.EAF.Ridge.Phy <- data.frame(matrix(unlist(aov_RidgeEAF_Phy$aov_results), nrow = length(aov_RidgeEAF_Phy$aov_results), byrow = T))
write.csv(unlist.EAF.Ridge.Phy, "EAF.Ridge.Phy.1wayanova.csv")
#none sig
#see google sheet for 1 way aov results

#end

#over both incs
winners_Phy <- bf.bi.long_Phylum %>%
  group_by(Phylum) %>%
  summarise(median.EAF = median.na_rm(median.EAF), median.NAi = median.na_rm(summed.NAi), IQR.NAi = IQR(summed.NAi, na.rm = TRUE), IQR.EAF = IQR(median.EAF, na.rm = TRUE))
#lab
winners_Phy_Lab <- bf.bi.long_Phylum %>%
  filter(Incubation == "Lab") %>%
  group_by(Phylum) %>%
  summarise(med.EAF = median.na_rm(median.EAF), median.NAi = median.na_rm(summed.NAi), IQR.NAi = IQR(summed.NAi, na.rm = TRUE), IQR.EAF = IQR(median.EAF, na.rm = TRUE))
#field
winners_Phy_Field <- bf.bi.long_Phylum %>%
  filter(Incubation == "Field") %>%
  group_by(Phylum) %>%
  summarise(med.EAF = median.na_rm(median.EAF), median.NAi = median.na_rm(summed.NAi), IQR.NAi = IQR(summed.NAi, na.rm = TRUE), IQR.EAF = IQR(median.EAF, na.rm = TRUE))



#library(ggridges)
Phylum.ridgeplot <- NAi.t20_Phylum %>%
  mutate(Phylum = fct_reorder(Phylum, summed.NAi, .fun = 'median', .na_rm = FALSE)) %>%
  ggplot(aes(x = summed.NAi, y = Phylum, fill = Management_Method, color = Management_Method)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.0001, alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  #scale_color_viridis(discrete=TRUE) +
  #scale_fill_viridis(discrete=TRUE) +
  scale_color_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  scale_fill_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  #scale_x_continuous(limits = c(0, 10), expand = c(0,0)) + #Bacillaceae has one Ridge Field rep around 20, Chitinophagaceae around 10
  labs(title = expression('%'~N~Assimilated), y = 'Phylum', x = '%', caption = "") + #removed "Phylum" for ggarrange
  theme_clean() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16),
        plot.caption = element_text(hjust = 0),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 12), 
        legend.title = element_blank(),
        axis.text.y = element_text(size = 12 ),
        axis.text.x = element_text(size = 16 ),
        plot.margin = margin(r = 20, unit = "pt"))

Phylum.ridgeplot

ggsave("PubGraphs/NAi_ridgeline_Phylum.png")
ggsave("PubGraphs/NAi_ridgeline_Phylum.svg")
#add sig levels later.. see word doc



Phylum.EAF.ridgeplot <- EAF.t20_Phylum %>%
  mutate(Phylum = fct_reorder(Phylum, median.EAF, .fun = 'median', .na_rm = FALSE)) %>%
  ggplot(aes(x = median.EAF, y = Phylum, fill = Management_Method, color = Management_Method)) +
  geom_density_ridges(scale = 1, rel_min_height = 0.0001, alpha = 0.4) +
  scale_x_continuous(limits = c(0, 0.8), expand = c(0,0)) +
  #scale_color_viridis(discrete=TRUE) +
  scale_color_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  scale_fill_manual(
    values = c(`Ridge Field` = "#92CA96",
               `Ridge Lab` = "#D1A1F5",
               `Valley Field` = "#085A12",
               `Valley Lab` = "#7F3FD0")) +
  # scale_fill_viridis(discrete=TRUE) +
  labs(title = expression(Relative~N~Assimilation~Rate), x = 'per 5 days', y = 'Phylum', caption = "") +
  theme_clean() +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", size = 16),
        plot.caption = element_text(hjust = 0),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        legend.text = element_text(size = 12), 
        legend.title = element_blank(),
        axis.text.y = element_text(size = 12 ),
        axis.text.x = element_text(size = 16 ),
        plot.margin = margin(r = 20, unit = "pt"))

Phylum.EAF.ridgeplot

ggsave("PubGraphs/EAF_ridgeline_Phylum.png") #purple and green saved as alt
ggsave("PubGraphs/EAF_ridgeline_Phylum.svg")


ggarrange(Phylum.EAF.ridgeplot, Phylum.ridgeplot, labels = c("A", "B"), ncol = 1, common.legend = TRUE, legend = "bottom")

ggsave("PubGraphs/20Ridgeline_Median_Phylum.png")
ggsave("PubGraphs/20Ridgeline_Median_Phylum.svg")




#sum microbes with less than XX ? abundance and see what percentage N they assimilated
#

#####
#NMDS

#upd Mar 13 2024
#USE CORRECTED FILES WITHOUT DUPLICATE
Lv6.filt.tube <- read.csv("Lv6-filt-tube_corrected.csv", header = TRUE)

tail(Lv6.filt.tube[,c(702:707)])

#will need to change "Farm" to "Site"
#make a vector of metadata names *change these to match yours*
data_meta_col_names <- c("index", "Site", "Method", "time_hrs", "Site_Method")
metadata_tube <- Lv6.filt.tube[, names(Lv6.filt.tube) %in% data_meta_col_names]
#remove all except index
matrix <- Lv6.filt.tube[,-c( 703:707)]

dim(matrix)
tail(matrix[,701:702]) #should end in bacillus
tail(matrix[,1:3])
#Make row names (what about "Control"? or time 0? Don't need to add?)
rownames(matrix) <- matrix$index
matrix <- matrix[,-1]


metadata_tube$incubation <- as.factor(metadata_tube$Method)
metadata_tube$management <- as.factor(metadata_tube$mFarm)

#already fixed below in "corrected Lv6-filt-tube"
#levels(metadata_tube$management) <- c('Valley', 'Ridge')
#levels(metadata_tube$management)
#If want to change from "None"
#levels(metadata_tube$incubation) <- c('Field', 'Lab', 'Initial')

#metadata_tube <- metadata_tube %>%
  #rename("Farm" = "management") %>%
 # rename("Method" = "incubation")

#metadata_tube$Farm_Incubation <- paste(metadata_tube$Farm,metadata_tube$Method)
metadata_tube$Farm_Incubation <- as.factor(metadata_tube$Farm_Method)

library(vegan)

NMDS <- metaMDS(matrix)



plot_df <- scores(NMDS, display = "sites") %>% 
  as.data.frame() %>% 
  rownames_to_column("site") %>% 
  full_join(metadata_tube, by = c("site" = "index"))



plot_nmds <- ggplot(plot_df, aes(x = NMDS1, y = NMDS2, color = Method, shape = Farm)) +
  geom_point(size = 6, alpha = 0.8) +
  # stat_ellipse(linetype = 2, size = 1) +
  scale_color_manual(values = c("Lab" = "#7F3FD0", "Field" ='#085A12', "Pre-Inc" = "darkgrey")) +
 # scale_color_manual(
   # values = c(`Ridge Field` = "#92CA96",
   #            `Ridge Lab` = "#D1A1F5",
    #           `Valley Field` = "#085A12",
     #          `Valley Lab` = "#7F3FD0")) +
  theme_clean() +
  labs(title = "NMDS") +
  theme(legend.position = "right",
        axis.title = element_text(face = "bold"),
        plot.caption = element_text(hjust = 0),
        axis.ticks = element_blank(),
        plot.background = element_blank(),
        panel.background = element_blank() ,
        panel.border = element_blank(),
        panel.grid.major.y = element_blank()) 
plot_nmds



ggsave("PubGraphs/NMDS_newcolors.png")
ggsave("PubGraphs/NMDS_newcolors.svg")



#Table 1 - 1-way anova for site properties
SiteProps <- read.csv("FxP.csv")

library(esquisse)
esquisser(SiteProps)
#rstatix

#add unique ID col 
SiteProps <- SiteProps %>%
  mutate(id = row_number())

pivot_wider(SiteProps, id_cols = id, names_from = colnames(SiteProps[2]), values_from = )

molten <- molten %>%
  mutate(id = row_number())

t <- pivot_wider(molten, id_cols = variable, names_from = Farm, values_from = value)

#data.table
molten <- melt(SiteProps)

SiteProps$Farm <- as.factor(SiteProps$Farm)

plantmass <- aov(data = SiteProps, Farm ~ plant_massg)


t.test(SiteProps$plant_massg ~ SiteProps$Farm)
t.test(SiteProps$plant_ht.m ~ SiteProps$Farm)
t.test(SiteProps$perc_15NPlantUptake ~ SiteProps$Farm)
t.test(SiteProps$C.NRhizosphere ~ SiteProps$Farm)
t.test(SiteProps$PercN.Rhizo ~ SiteProps$Farm)
t.test(SiteProps$PercC.Rhizo ~ SiteProps$Farm)
t.test(SiteProps$PercS.Rhizo ~ SiteProps$Farm)
t.test(SiteProps$Perc.Moisture ~ SiteProps$Farm)
t.test(SiteProps$SOM ~ SiteProps$Farm)
t.test(SiteProps$pH ~ SiteProps$Farm)
t.test(SiteProps$WHC ~ SiteProps$Farm)

#% higher among all Phyla
allphy <- bf.bi.long_Phylum %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(summed.NAi), meanNAi = mean(summed.NAi, na.rm = TRUE), medianEAF = median.na_rm(median.EAF), meanEAF = mean(median.EAF, na.rm = TRUE))
#mean EAF 15% higher lab
#median NAi 8% higher in lab

#Ridge
allphyridge <- bf.bi.long_Phylum %>%
  filter(Management == "Ridge") %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(summed.NAi), meanNAi = mean(summed.NAi, na.rm = TRUE), medianEAF = median.na_rm(median.EAF), meanEAF = mean(median.EAF, na.rm = TRUE))

#meanEAF 24% higher lab 
#median NAi 17% higher lab

#Valley
allphyvalley <- bf.bi.long_Phylum%>%
  filter(Management == "Valley") %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(summed.NAi), meanNAi = mean(summed.NAi, na.rm = TRUE), medianEAF = median.na_rm(median.EAF), meanEAF = mean(median.EAF, na.rm = TRUE))

#meanEAF 4% higher lab 
#median NAi 6% higher field

#% higher among top Phylum only 

esquisser(NAi.t20_Phylum)
total <- NAi.t20_Phylum %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(summed.NAi), meanNAi = mean(summed.NAi, na.rm = TRUE), medianEAF = median.na_rm(median.EAF), meanEAF = mean(median.EAF, na.rm = TRUE))

esquisser(total)
#EAF mean ~14% increase in the lab
#EAF median ~4% incease in lab
#median NAi ~2% increase in the field
#mean NAi 0.4% increase in lab

#across all Phylum
total <- NAi.t20_Phylum %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(summed.NAi), meanNAi = mean(summed.NAi, na.rm = TRUE), medianEAF = median.na_rm(median.EAF), meanEAF = mean(median.EAF, na.rm = TRUE))


#% change among top genera only for NAi


total2 <- genus.NAi.t20 %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))

#EAF mean ~18% increase in lab
#EAF median ~8% increase in lab
#NAi mean ~3% increase in field
#NAi median ~0.8% increase in field


#% change among top genera for each site

#Ridge
total4 <- genus.NAi.t20 %>%
  filter(Management == "Ridge") %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))

#EAF mean ~28% increase in lab
#EAF median ~21% increase in lab
#NAi mean ~2% increase in field
#NAi median ~6% increase in lab


#Valley

total5 <- genus.NAi.t20 %>%
  filter(Management == "Valley") %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))

#EAF mean ~4% increase in lab
#EAF median ~2% increase in field
#NAi mean ~5% increase in field
#NAi median ~4% increase in field


#% change among all genera

total3 <- bf.bi.long.filt %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))

#EAF mean 16% increase in lab
#EAF median ~13% increase in lab
#NAi mean ~1.3% higher in field
#NAi median ~3% increase in lab

#Valley
total6 <- bf.bi.long.filt %>%
  filter(Management == "Valley") %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))

#EAF mean ~7% increase in lab
#EAF median ~4% increase in lab
#NAi mean ~0.7% increase in lab
#NAi median ~5% increase in lab ***NAi median and mean VERY different. This is true for all of these

#Ridge
total7 <- bf.bi.long.filt %>%
  filter(Management == "Ridge") %>%
  group_by(Incubation) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))

#EAF mean ~23% increase in lab
#EAF median ~26% increase in lab
#NAi mean ~3% increase in field
#NAi median ~0.07% increase in lab



#genus % change 

EAF.t20.sum <- genus.EAF.t20 %>%
  group_by(Genus, Incubation, Management) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))

NAi.t20.sum <- genus.NAi.t20 %>%
  group_by(Genus, Incubation, Management) %>%
  summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))


#Phylum % change
NAi.t20.sum.Phy <- NAi.t20_Phylum %>%
  group_by(Phylum, Incubation, Management) %>%
  summarise(medianNAi = median.na_rm(summed.NAi), meanNAi = mean(summed.NAi, na.rm = TRUE), medianEAF = median.na_rm(median.EAF), meanEAF = mean(median.EAF, na.rm = TRUE))


EAF.t20.sum.Phy <- EAF.t20_Phylum %>%
  group_by(Phylum, Incubation, Management) %>%
  summarise(medianNAi = median.na_rm(summed.NAi), meanNAi = mean(summed.NAi, na.rm = TRUE), medianEAF = median.na_rm(median.EAF), meanEAF = mean(median.EAF, na.rm = TRUE))



write.csv(EAF.t20.sum.Phy, "EAF.t20.sum.Phy.csv")
write.csv(NAi.t20.sum.Phy, "NAi.t20.sum.Phy.csv")

write.csv(EAF.t20.sum, "EAF.t20.sum.genus.csv")
write.csv(NAi.t20.sum, "NAi.t20.sum.genus.csv")


#plot high assimilating taxa similarly to fig 2
#Chanso assign



wide.t20.NAi <- pivot_wider(genus.NAi.t20, names_from = "Incubation", values_from = c("NAi", "EAF"), id_cols = c("Management", "tube", "tube", "Genus", "Phylum", "Class", "Order", "Family", "taxon"))
esquisser(wide.t20.NAi)


wide.t20.EAF <- pivot_wider(genus.EAF.t20, names_from = "Incubation", values_from = c("NAi", "EAF"), id_cols = c("Management", "tube", "tube", "Genus", "Phylum", "Class", "Order", "Family", "taxon"))

ggplot(wide.t20.EAF) +
 aes(x = EAF_Lab, y = EAF_Field, colour = Management) +
 geom_point(shape = "circle", 
 size = 1.5) +
#scale_x_continuous(limits = c(0, 8), expand = c(0,0)) +
 scale_color_hue(direction = 1) +
 labs(title = "Highest EAF only Field v Lab") +
 theme_minimal() +
 theme(legend.position = "bottom") +
 facet_wrap(vars(Management))

ggplot(wide.t20.NAi) +
 aes(x = NAi_Lab, y = NAi_Field, colour = Management) +
 geom_point(shape = "circle", 
 size = 1.5) +
 scale_y_continuous(limits = c(0, 6), expand = c(0,0)) + #removes outlier bacillus
 scale_color_hue(direction = 1) +
 labs(title = "NAi Top Assim Field v Lab") +
 theme_minimal() +
 theme(legend.position = "bottom") +
 facet_wrap(vars(Management))


#looks the same as fig 2

#taxa with less than 0.XX % rel abundance, assimilated X% of the N

ncopies.tube.taxa <- read.csv("ncopies.tube.taxa.csv")
ncopies.tube.taxa <- subset(ncopies.tube.taxa, rel.ncopies > 0.0001) #was 0.0005. changed 29 Mar 2023 for graphing purposes
length(unique(ncopies.tube.taxa$taxon)) #692 -  We initially started with 703 taxa


#change headers to match bf.bi.long

ncopies.tube.taxa <- ncopies.tube.taxa %>%
  rename("Management" = "management") %>%
  rename("Incubation" = "incubation") 

ncopies.tube.taxa$Management <- as.factor(ncopies.tube.taxa$Management)
ncopies.tube.taxa$Incubation <- as.factor(ncopies.tube.taxa$Incubation)
library(forcats)

ncopies.tube.taxa.new <- ncopies.tube.taxa %>%
  filter(Incubation == c("Field", "Lab")) %>%
  mutate(Management = recode(Management, Conventional = "Valley",
                    Organic = "Ridge" ))

#esquisser(ncopies.tube.taxa.new) # heavily skewed like NAi


#summarise each genus by site and method (no reps)
#both data frames (ncopies one will remove time 0s)
#so match

sum.ncopies <- ncopies.tube.taxa.new %>%
  group_by(Management, Incubation, taxon) %>%
 summarise(median.rel.ncopies = median(rel.ncopies))
  
bf.bi.long.filt.sum <-bf.bi.long.filt %>%
    group_by(Management, Incubation, taxon, Management_Method, Genus, Phylum) %>%
summarise(medianNAi = median.na_rm(NAi), meanNAi = mean(NAi, na.rm = TRUE), medianEAF = median.na_rm(EAF), meanEAF = mean(EAF, na.rm = TRUE))


bf.bi.long.filt_wcopies <- merge(bf.bi.long.filt.sum,
                                 sum.ncopies, by = c("taxon", "Incubation", "Management"))


sum.relabund <- bf.bi.long.filt_wcopies %>%
  group_by(Management_Method) %>%
  get_summary_stats(median.rel.ncopies)

#median - IQR/2 = Q1 (0 in summary stat not helpful...)

#Do for each site and each inc. 

org.str <- bf.bi.long.filt_wcopies %>%
  count(Management_Method)
## of obs for each
#Ridge Field 398
#Ridge Lab 396
#Valley Field 412
#Valley Lab 411

pt5percabund <- bf.bi.long.filt_wcopies %>%
  filter(median.rel.ncopies < 0.005) %>% #those with less than 0.5% median rel. abund
  count(Management_Method) 

#Ridge Field 344
#Ridge Lab 344
#Valley Field 359
#Valley Lab 358
  
#344/398 = 0.86
#344/396 = 0.87
#359/412 = 0.87
#358/411 = 0.87

onepercabund <- bf.bi.long.filt_wcopies %>%
  filter(median.rel.ncopies < 0.01) %>% #those with less than 1% median rel. abund
  count(Management_Method)

#Ridge Field 379
#Ridge Lab 376
#Valley Field 395
#Valley Lab 390

#379/398 = .95
#376/396 = .95
#395/412 = .96
#390/411 = .95


t0.5 <- bf.bi.long.filt_wcopies %>%
  filter(median.rel.ncopies < 0.005) %>% 
  group_by(Management_Method) %>%
  summarise(summed.NAi = sum(medianNAi))
#33-38% of N assimilated by those with rel. abund under .5%

#33.5/90.1 = .37
#34.7/91.5 = .38
#37.5/94.5 = 0.40
#38.6/93.1 = .41

t1 <- bf.bi.long.filt_wcopies %>%
  filter(median.rel.ncopies < 0.01) %>% #only have those with rel abund < 1% to start)
  group_by(Management_Method) %>%
  summarise(summed.NAi = sum(medianNAi))
#between 55-63% by those under 1%

over1 <- bf.bi.long.filt_wcopies %>%
  filter(median.rel.ncopies > 0.01) %>% 
  group_by(Management_Method) %>%
  summarise(summed.NAi = sum(medianNAi))
#31-36% over 1% ... missing ~20% of isotope in summaries


bf.bi.long.filt_wcopies %>%
  group_by(Management_Method) %>%
  summarise(summed.NAi = sum(medianNAi))
#no filtering = 
#90.1% Ridge field
#91.5% Ridge Lab
#94.5% Valley Field
#93.1% Valley Lab



#Bacillus lab vs field rel abundance
bacillus <- ncopies.tube.taxa %>%
  filter(taxon == "d__Bacteria.p__Firmicutes.c__Bacilli.o__Bacillales.f__Bacillaceae.g__Bacillus") 


#F v L table with gene copies

FvL <- read.csv("FvL.csv")

summaryFvL <- FvL %>%
  group_by(Farm, Incubation) %>%
  get_summary_stats()

copiesaov <- aov(data = FvL, X16Scopies_gsoil ~ Farm * Incubation)
summary(copiesaov)

#               Df    Sum Sq   Mean Sq F value  Pr(>F)   
#Farm             1 3.999e+17 3.999e+17   3.904 0.05977 . 
#Incubation       1 5.578e+17 5.578e+17   5.445 0.02833 * 
#  Farm:Incubation  1 1.300e+18 1.300e+18  12.695 0.00158 **
 # Residuals       24 2.458e+18 1.024e+17                   
---
  
TukeyHSD(copiesaov, "Farm:Incubation")

#try sep for each farm 1 way anova

Valleyaov <- FvL %>%
  filter(Farm == "Valley") 
 
summary(aov(data = Valleyaov,X16Scopies_gsoil ~ Incubation))

Ridgeaov <- FvL %>%
  filter(Farm == "Ridge") 

summary(aov(data = Ridgeaov,X16Scopies_gsoil ~ Incubation))

#still not sig

esquisser(FvL)


#Valley lab has outliers??

ggplot(FvL) +
 aes(x = X16Scopies_gsoil, fill = Incubation) +
 geom_histogram(bins = 62L) +
 scale_fill_hue(direction = 1) +
 theme_minimal()


ggplot(FvL) +
 aes(x = X16Scopies_gsoil, y = Incubation, fill = Incubation) +
 geom_boxplot() +
 scale_fill_hue(direction = 1) +
 theme_minimal() +
 facet_wrap(vars(Farm))



copiesaov <- aov(data = FvL, X16Scopies_gsoil ~ Farm * Incubation)
summary(aov(data = FvL, CO2Flux ~ Farm * Incubation))

TukeyHSD(aov(data = FvL, CO2Flux ~ Farm * Incubation), "Farm:Incubation")


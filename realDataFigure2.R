options(timeout = 600)

library(mutoss) ## function BH

## library(hierarchicalFDR) 
## from Schildknecht, Tabelow and Dickhaus (2016) 
## we only use the part in the original two-step procedure 
## "hierasymptkappa()" that calculates the family p-values 
## this part is extracted and put into the function "calculate_family_p_value()"

calculate_family_p_value <- function(pValues, kappa = 0.001) {
  fpValues <- sort(pValues)
  m <- length(pValues)
  u <- as.integer(kappa*m) + 1
  conjpValues <- numeric(m-u+1)
  for (i in 1:(m-u+1)) {
    conjpValues[i] <- fpValues[[(u-1+i)]]*(m-u+1)/i
  }
  min(conjpValues) ## family p-value
}


sink(file = "Aparc-combination-results-Variante2.txt")

data <- read.csv("https://archive.wias-berlin.de/servlets/MCRFileNodeServlet/wias_derivate_00003786/pvalueFromRealData.csv")[, c(-6, -8)] 

## including p-values ## without test statistics
data$Participant <- factor(data$Participant)
data$Aparc <- factor(data$Aparc)


kappa <- 0.001 
## default value 
## tuning parameter introduced in the function "hierasymptkappa()" 
## influences two things 
## 1) influences the threshold of active voxels needed for a family 
##    hypotheses to be not filtered out due to containing only "noise" 
## 2) influences the significance level needed for rejecting family p-values


Aparc_labels <- levels(data$Aparc)
m <- length(Aparc_labels)
participants <- levels(data$Participant)
n <- length(participants) ## 11 participants in PC2
family_p_values.sorted_by_participants <- list() 
## will contain a vector of all Aparc p-values for each participant
family_p_values.sorted_by_Aparc_labels <- list() 
## will contain a vector of Aparc p-values over all participants for each Aparc label


t <- Sys.time()
lower_bound <- 1.0e-15;
for(i in 1:n) {
  data_i <- data[data$Participant == participants[i], -1]
  data_i$PC2_Maske_p[data_i$PC2_Maske_p==0] <- lower_bound; 
  ## p-values can be zero due to rounding otherwise 
  ## few p-values with value zero can dominate the combined test result 
  ## => family p-value = 0 and Fisher X2 = Inf
  data_i$PC2_Compr_p[data_i$PC2_Compr_p==0] <- lower_bound;
  
  
  #######################################################
  #### Step 1: (Bottom-Up) Comprehension versus Rest ####
  #######################################################
  
  ## in this step, we construct a filter by searching for significant voxels 
  ## when we compare comprehention with rest 
  alpha_step1 <- 0.05
  p_values_step1 <- data_i$PC2_Maske_p ## two-sided p-values
  BH_step1 <- BH(p_values_step1, alpha_step1, silent = TRUE)
  adjusted_p_values_step1 <- BH_step1$adjPValues ## adjust the p-values using Benjamini-Hochberg
  
  active_p_values_step1 <- BH_step1$rejected
  
  
  #########################################################
  #### Step 2: Bottom-Up (Comprehension) versus Syntax ####
  #########################################################
  
  ## in this step, we apply the filter from step 1 ##
  alpha_step2 <- 0.05 
  ## changed from 0.01 to 0.05 in the second study PC2 since there are less participants
  p_values_step2 <- data_i$PC2_Compr_p[active_p_values_step1] 
  ## filtered by significant p-values from step1 
  ## two-sided p-values
  Aparc_step2 <- droplevels(data_i$Aparc[active_p_values_step1]) 
  ## Aparc labels filtered by significant p-values from step1 
  ## Aparc labels are dropped when they contain only non-significant voxel from step1;
  
  
  ###############################################################################
  #### Step 3: Data Analysis using Schildknecht's function hierasymptkappa() ####
  ###############################################################################
  
  ## in this step, we calculate the family p-values (i.e. the Aparc p-values) using 
  ## the filtered p-values from step 2
  family_p_values <- numeric(m) + NA_real_
  for(j in 1:m) {
    index_j <- which(Aparc_step2 == Aparc_labels[j])
    if(length(index_j) > 0)
      family_p_values[j] <- calculate_family_p_value(p_values_step2[index_j], kappa)
  }
  family_p_values.sorted_by_participants <- c(family_p_values.sorted_by_participants, list(family_p_values = family_p_values))
}
print(Sys.time() - t)


#####################################################
#### Step 4: Combine p-Values by Fisher's Method ####
#####################################################

## in this step, we combine the Aparc p-values over all participants 
## and use them to test the Aparc labels 
Fisher_X2 <- numeric(m) + NA_real_
crit_values <- numeric(m) + NA_real_
rejected <- numeric(m) + NA_real_
for(j in 1:m) {
  p_values_ <- numeric(n)
  for(i in 1:n) {
    p_values_[i] <- family_p_values.sorted_by_participants[[i]][j]
  }
  family_p_values.sorted_by_Aparc_labels <- c(family_p_values.sorted_by_Aparc_labels, list(family_p_values = p_values_))
  p_values_ <- na.omit(p_values_)
  if(length(p_values_) > 0) {
    Fisher_X2[j] <- -2*sum(log(p_values_)) ## Fisher's test statistic X^2
    crit_values[j] <- qchisq(1-alpha_step2*kappa, df = 2*length(p_values_)) ##
    rejected[j] <- Fisher_X2[j] > crit_values[j]
  }
}


###########################
#### Print the results ####
###########################

print(list(
  Aparc_labels = Aparc_labels,
  Fisher_X2 = Fisher_X2,
  crit_values = crit_values,
  rejected = rejected,
  Aparc_labels_rejected = na.omit(Aparc_labels[as.logical(rejected)])
))


closeAllConnections()


###########################
# Create ggplot image     #
###########################
Aparc_labels_rejected = na.omit(Aparc_labels[as.logical(rejected)]);

library(ggplot2)

family_p_values_rejected <- as.data.frame(family_p_values.sorted_by_participants)[!is.na(rejected) & rejected == 1, ]
names(family_p_values_rejected) <- participants

ggplot_data_frame  <- data.frame(Aparc_label = Aparc_labels_rejected, family_p_values = family_p_values_rejected[, 1])
for (i in 2:n) {ggplot_data_frame <- rbind(ggplot_data_frame, cbind(data.frame(Aparc_label = Aparc_labels_rejected, family_p_values = family_p_values_rejected[, i])))}

gg <- ggplot(ggplot_data_frame, aes(x = family_p_values,y = Aparc_label)) + geom_point()
ggsave("./2021-01-06-rejected-family-p-values-Plot.png", width = 5, height = 4)
print(gg);


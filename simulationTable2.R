library(mutoss) ## function BH
library(neuRosim)
library(fmri)
library(oro.nifti)

##################################################
# GenerateFMRI data                              #
##################################################
generate_fMRI_data <-function(nsubjects, my_SNR=7.0)
{
 # for each subject create a datatset with a blob.

 # spatial dimension
 ddim <- c(20, 20, 20)


 # temporal design
 nscan <- 195
 TR <- 2
 on1 <- c(16, 76, 136) * TR
 on2 <- c(46, 106, 166) * TR
 on3 <- c(16, 46, 76, 106, 136, 166) * TR
 onsets <- list(list(on1), list(on2), list(on3))
 dur <- list(list(15*TR), list(15*TR), list(15*TR))

 hrf1 <- fmri.stimulus(nscan, c(16, 76, 136), 15, 2, type = "gamma")
 hrf2 <- fmri.stimulus(nscan, c(46, 106, 166), 15, 2, type = "gamma")
 xdesign <- fmri.design(cbind(hrf1, hrf2), order = 2)

 design <- simprepTemporal(totaltime = nscan*TR, 
                           regions = 3, 
                           onsets = onsets, 
                           durations = dur, 
                           TR = TR, 
                           effectsize = list(40, 20, 40),
                           hrf = "double-gamma")

 ttt <- matrix(0, nsubjects * prod(ddim), 9)

 i <- 0
 for (subj in 1: nsubjects) {

   jitter1 <- as.integer(4*runif(3)-1)
   rc1 <- c(5, 5, 5) + jitter1
   rc2 <- c(5, 5, 5) + jitter1
   jitter2 <- as.integer(4*runif(3)-1)
   rc3 <- c(15, 15, 15) + jitter2
  
   regions <- simprepSpatial(regions = 3, 
                             coord = list(rc1, rc2, rc3), 
                             radius = c(2, 2, 4), 
                             form = "sphere")
  
   ds.sim <- simVOLfmri(design = design, 
                        image = regions, 
                        base = 100, 
                        dim = ddim,
                        noise = "mixture", 
                        type = "rician", 
                        spat = "gaussRF",
                        weights = c(0.1, 0.3, 0.2, 0.1, 0.1, 0.2), 
                        rho.temp = 0.3, 
                        SNR = my_SNR)
  
   # plot(ds.sim[1, 1, 1,], type = "l",xlab="Frame",ylab="BOLD signal",main="Background")
   # plot(ds.sim[5, 5, 5,], type = "l",xlab="Frame",ylab="BOLD signal",main="Active region")
   # plot(ds.sim[15, 15, 15,], type = "l",xlab="Frame",ylab="BOLD signal",main="Active region")
  
   nii <- as.nifti(ds.sim)
   alpha <- 0.05  

   # analyze each of the dataset and create p-values

   ds <- oro2fmri(nii, setmask = FALSE)
   spm <- fmri.lm(ds, xdesign, contrast = c(1, 1))
   cluster_p_value1 <- fmri.cluster(spm, alpha, ncmin=5, ncmax = 15, verbose = FALSE); #ncmin was 5 before
   
  
   theta <- spm$cbeta/sqrt(spm$var)
   pvalue <- pt(-theta, spm$df)
   #pvalue[pvalue >= alpha] <- NA
   dim(pvalue) <- spm$dim[1:3]
  
   spm2 <- fmri.lm(ds, xdesign, contrast = c(1, -1))
   cluster_p_value2 <- fmri.cluster(spm2, alpha, ncmin=5, ncmax = 15, verbose = FALSE); #ncmin was 5 before
   theta2 <- spm2$cbeta/sqrt(spm2$var)
   pvalue2 <- pt(-theta2, spm2$df)
   #pvalue2[pvalue2 >= alpha] <- NA
   dim(pvalue2) <- spm2$dim[1:3]
  
   # build table for hierarchical analysis
   for (xd in 1:ddim[1]) {
     for (yd in 1:ddim[2]) {
       for (zd in 1:ddim[3]) {
         i <- i +1
         if (zd > 10) {
           if (yd > 10) {
             if (xd > 10) {
               aparc <- 1L
             } else {
               aparc <- 2L
             }
           } else {
             if (xd > 10) {
               aparc <- 3L
             } else {
               aparc <- 4L
             }
           }
         } else {
           if (yd > 10) {
             if (xd > 10) {
               aparc <- 5L
             } else {
               aparc <- 6L
             }
           } else {
             if (xd > 10) {
               aparc <- 7L
             } else {
               aparc <- 8L
             }
           }
         }
	    cluster_p_value1$pvalue[is.na(cluster_p_value1$pvalue)] <- pvalue[is.na(cluster_p_value1$pvalue)];
         cluster_p_value2$pvalue[is.na(cluster_p_value2$pvalue)] <- pvalue2[is.na(cluster_p_value2$pvalue)];
         ttt[i, ] <- c(as.integer(subj), #subject
                       as.integer(xd), # x-coord
                       as.integer(yd), # y-coord
                       as.integer(zd), # z-coord
                       aparc, # APARC
                       pvalue[xd, yd, zd], # pvalue first contrast
                       pvalue2[xd, yd, zd], # pvalue second contrast
                       cluster_p_value1$pvalue[xd, yd, zd], # cluster pvalue first contrast
                       cluster_p_value2$pvalue[xd, yd, zd]) # cluster pvalue second contrast
       }
     }
   } 
 }
 return(ttt);
}


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



MC <- 500; #Number of Monte Carlo runs
Number_of_Aparc_Labels <- 8;
my_n <- 11; #number of subjects
signal_to_noise <- 1.5;

kappa <- 0.01; #10 out of 1000
## tuning parameter introduced in the function "hierasymptkappa()" 
## influences two things 
## 1) influences the threshold of active voxels needed for a family 
##    hypotheses to be not filtered out due to containing only "noise" 
## 2) influences the significance level needed for rejecting family p-values

sink("out.txt")

for (signal_to_noise in c(0.75, 1.0, 1.25, 1.5, 1.75)) {
#for (signal_to_noise in c(1.5)) {
    
REGION_rejected_over_MC <-rep(0L, Number_of_Aparc_Labels); 
CLUSTER_rejected_over_MC <-rep(0L, Number_of_Aparc_Labels);
VOXEL_rejected_over_MC <-rep(0L, Number_of_Aparc_Labels);

#############################################
#     Loop over simulation runs             #
#############################################
for(run in (1:MC))
{
 print(run);
 data_matrix <- generate_fMRI_data(nsubjects=my_n, my_SNR=signal_to_noise);
 data <- data.frame(data_matrix);
 ## including p-values ## without test statistics
 data$Participant <- factor(data$X1);
 data$X <- data$X2;
 data$Y <- data$X3;
 data$Z <- data$X4;
 data$Aparc <- factor(data$X5);
 data$p_values.contrast1 <- data$X6;
 data$p_values.contrast2 <- data$X7;

 Aparc_labels <- levels(data$Aparc);
 m <- length(Aparc_labels);
 participants <- levels(data$Participant);
 n <- length(participants); ## 11 participants in PC2

 family_p_values.sorted_by_participants <- list() 
 ## will contain a vector of all Aparc p-values for each participant
 family_p_values.sorted_by_Aparc_labels <- list() 
 ## will contain a vector of Aparc p-values over all participants for each Aparc label

 #t <- Sys.time()
 lower_bound <- 1.0e-15;
 for(i in 1:n) {
   data_i <- data[data$Participant == participants[i], -1];
   data_i$p_values.contrast1[data_i$p_values.contrast1==0] <- lower_bound; 
   ## p-values can be zero due to rounding otherwise 
   ## few p-values with value zero can dominate the combined test result 
   ## => family p-value = 0 and Fisher X2 = Inf
   data_i$p_values.contrast2[data_i$p_values.contrast2==0] <- lower_bound;


   #######################################################
   #### Step 1: (Bottom-Up) Comprehension versus Rest ####
   #######################################################

   ## in this step, we construct a filter by searching for significant voxels 
   ## when we compare comprehention with rest 
   alpha_step1 <- 0.05
   p_values_step1 <- data_i$p_values.contrast1 ## two-sided p-values
   BH_step1 <- BH(p_values_step1, alpha_step1, silent = TRUE)
   adjusted_p_values_step1 <- BH_step1$adjPValues ## adjust the p-values using Benjamini-Hochberg

   active_p_values_step1 <- BH_step1$rejected


   #########################################################
   #### Step 2: Bottom-Up (Comprehension) versus Syntax ####
   #########################################################

   ## in this step, we apply the filter from step 1 ##
   alpha_step2 <- 0.05 
   ## changed from 0.01 to 0.05 in the second study PC2 since there are less participants
   p_values_step2 <- data_i$p_values.contrast2[active_p_values_step1]; 
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
 #print(Sys.time() - t)


 #####################################################
 #### Step 4: Combine p-Values by Fisher's Method ####
 #####################################################

 ## in this step, we combine the Aparc p-values over all participants 
 ## and use them to test the Aparc labels 
 Fisher_X2 <- numeric(m) + NA_real_
 crit_values <- numeric(m) + NA_real_
 #rejected <- numeric(m) + NA_real_
 rejected <- rep(0L, m);
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


 ##########################################
 #### Store the results (Aparc method) ####
 ##########################################
 for(j in (1:m))
 {
  if(rejected[j])
  {
   REGION_rejected_over_MC[j] <- REGION_rejected_over_MC[j] + 1;
  }
 }

 ############################
 #       Mask method        #
 ############################
 Z_scores_contrast1 <- qnorm(data_matrix[, 6]);
 dim(Z_scores_contrast1) <- c(8000, my_n);
 average_Z_scores_contrast1 <- rowMeans(Z_scores_contrast1);
 aggregated_p_values_contrast1 <- pnorm(average_Z_scores_contrast1, sd=1/sqrt(n)); 
 aggregated_BH_contrast1 <- BH(aggregated_p_values_contrast1, alpha_step1, silent = TRUE);
 active_aggregated_p_values_contrast1 <- aggregated_BH_contrast1$rejected;
 
 Z_scores_constrast2 <- qnorm(data_matrix[, 7]);
 dim(Z_scores_constrast2) <- c(8000, 11);
 average_Z_scores_constrast2 <- rowMeans(Z_scores_constrast2);
 aggregated_p_values_contrast2 <- pnorm(average_Z_scores_constrast2, sd=1/sqrt(n)); 
 masked_aggregated_p_values_contrast2 <- aggregated_p_values_contrast2[active_aggregated_p_values_contrast1];
 Aparc_labels_subject1 <- data_matrix[1:8000, 5];
 masked_Aparc_labels <- Aparc_labels_subject1[active_aggregated_p_values_contrast1];

 alpha_aparc <- 0.05;
 length_mask <- length(masked_Aparc_labels);
 for(alabel in (1:m))
 {
  masked_p_values_label <- masked_aggregated_p_values_contrast2[masked_Aparc_labels == alabel];
  if(sum(masked_p_values_label < alpha_aparc / length_mask) > 1000 * kappa)
  {
   VOXEL_rejected_over_MC[alabel] <- VOXEL_rejected_over_MC[alabel] + 1;
  } 
 }

 ############################
 #       Cluster method     #
 ############################
 Cluster_Z_scores_constrast1 <- qnorm(data_matrix[, 8]);
 dim(Cluster_Z_scores_constrast1) <- c(8000, 11);
 average_Cluster_Z_scores_contrast1 <- rowMeans(Cluster_Z_scores_constrast1);
 aggregated_Cluster_p_values_contrast1 <- pnorm(average_Cluster_Z_scores_contrast1, sd=1/sqrt(n)); 
 aggregated_Cluster_BH_contrast1 <- BH(aggregated_Cluster_p_values_contrast1, alpha_step1, silent = TRUE);
 active_aggregated_Cluster_p_values_contrast1 <- aggregated_Cluster_BH_contrast1$rejected;
 
 Cluster_Z_scores_constrast2 <- qnorm(data_matrix[, 9]);
 dim(Cluster_Z_scores_constrast2) <- c(8000, 11);
 average_Cluster_Z_scores_constrast2 <- rowMeans(Cluster_Z_scores_constrast2);
 aggregated_Cluster_p_values_contrast2 <- pnorm(average_Cluster_Z_scores_constrast2, sd=1/sqrt(n)); 
 masked_aggregated_Cluster_p_values_contrast2 <- aggregated_Cluster_p_values_contrast2[active_aggregated_Cluster_p_values_contrast1];
 Aparc_labels_subject1 <- data_matrix[1:8000, 5];
 masked_Cluster_Aparc_labels <- Aparc_labels_subject1[active_aggregated_Cluster_p_values_contrast1];

 alpha_aparc <- 0.05;
 length_Cluster_mask <- length(masked_Cluster_Aparc_labels);
 for(alabel in (1:m))
 {
  masked_Cluster_p_values_label <- masked_aggregated_Cluster_p_values_contrast2[masked_Cluster_Aparc_labels == alabel];
  if(sum(masked_Cluster_p_values_label < alpha_aparc / length_Cluster_mask) > 1000 * kappa)
  {
   CLUSTER_rejected_over_MC[alabel] <- CLUSTER_rejected_over_MC[alabel] + 1;
  } 
 }

} # end of main simulation loop

print(paste("SNR", signal_to_noise))
print(REGION_rejected_over_MC)
print(CLUSTER_rejected_over_MC)
print(VOXEL_rejected_over_MC)

}

closeAllConnections()
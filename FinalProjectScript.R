#This is the code to accompany Jeffrey Saxon's CYO Capstone Project for HarvardX PH125.9x


#This initial code will download the relevant files from Github and set up the large date file and pedigree file used in later analysis

if(!require(tidyverse)){
  install.packages("tidyverse")
}

if(!require(caret)){
  install.packages("caret")
}

if(!require(e1071)){
  install.packages("e1071")
}

library(tidyverse)
library(caret)
library(e1071)

setwd("/home/jeffrey/R/CapstonePaternity2")

pedigree_filepath <- 'https://github.com/JeffSaxon/PH125.9x-Capstone-Saxon/raw/master/Pedigree.csv'
saved_results_filepath <- 'https://github.com/JeffSaxon/PH125.9x-Capstone-Saxon/raw/master/saved_results'
saved_data_filepath <- 'https://github.com/JeffSaxon/PH125.9x-Capstone-Saxon/raw/master/saved_data.csv'

download.file(pedigree_filepath, 'Saxon_Capstone_pedigree')
download.file(saved_results_filepath, 'Saxon_Capstone_saved_results')
download.file(saved_data_filepath, 'Saxon_Capstone_saved_data')

pedigree<-read_tsv("Saxon_Capstone_pedigree")
individuals<-c(pedigree$person_1[1:15], pedigree$person_2[1:10])
parents<-pedigree[1:10,1:2]

large_dat<-read_csv('Saxon_Capstone_saved_data')
allRsids<-unique(large_dat$rsid)
rsidSample<- allRsids

#This code creates an unpopulated dataframe to hold, for each of 300 relationships among 25 people, the number of identical genotype values, the number of wholly-different genotype values, and whether or not the relationship is one of a parent and child 
relationships<- t(combn(individuals,2))
dat<-data.frame("person_1"=relationships[,1], "person_2"=relationships[,2])
complete_data<- dat %>% mutate(diff_0=0, diff_1=0, diff_2=0, isParent=FALSE)

#This helper function receives an individual's ID and returns the ID number of its family
get_family <- function (x) {
  result<-0
  for (i in 1:10){
    if (x==as.character(parents[i,1]) | x==as.character(parents[i,2])) {
      result <- i
    }
  }
  result
}

#This helper function receives the IDs of two individuals and uses the parents file to return a boolean value as to whether their relationship is of parent and child
is_parent <- function (x,y) {
  result<-FALSE
  for (i in 1:nrow(parents)){
    if ((x==as.character(parents[i,1]) & y==as.character(parents[i,2])) | (y==as.character(parents[i,1]) & x==as.character(parents[i,2]))) {
      result<-TRUE
    }
  }
  result
}

#This helper function receives two genotype values and returns the difference between them as an integer 
gt_diff <- function (x,y) {
  alleles<-c("A", "C", "T", "G")
  sum(abs(str_count(x,alleles)-str_count(y,alleles)))/2
}

#This helper function receives IDs for two individuals and a vector of RSID values and returns a vector of the differences between their genotype values
diff_vector <- function (x, y, z) {
  a<-large_dat %>% filter(as.character(person)==x) %>% filter(rsid %in% z)
  b<-large_dat %>% filter(as.character(person)==y) %>% filter(rsid %in% z)
  c<-a$genotype
  d<-b$genotype
  if (identical(a$rsid, b$rsid)==FALSE) {
    print("error")
  }
  e<-1:length(c)
  for (i in 1:length(c)) {
    e[i]<-gt_diff(as.character(c[i]),as.character(d[i]))
  }
  e
}

#This helper function receives a data frame and vector of RSID values and updates the data frame for differences among the genotype of individuals in the data frame.
update_data <- function (x, y){
  for (i in 1:nrow(x)) {
    person1<-as.character(x[i,1])
    person2<-as.character(x[i,2])
    diffVector<-diff_vector(person1, person2, y)
    x[i,3]<-sum(diffVector==0)
    x[i,4]<-sum(diffVector==1)
    x[i,5]<-sum(diffVector==2)
    x[i,6]<-is_parent(person1, person2)
  }
  x
}


#This code identifies the distribution of genotypes in the large data file
large_dat %>% group_by(genotype) %>% summarize(n=n())

#This code counts the number of RSIDs that are homozygous identical among all of the individuals in the large data file.
a<-large_dat %>% group_by(rsid) %>% summarize(g_count=sum(str_count(genotype, "G")), t_count=sum(str_count(genotype, "T")), c_count=sum(str_count(genotype, "C")), a_count=sum(str_count(genotype, "A")))
sum(a$g_count==52) + sum(a$a_count==52) + sum(a$t_count==52) + sum(a$c_count==52)




#This code uses the helper functions to prepare a data frame of differences in genotype for sample sizes of 5, 40 and 500.
rsidSample<- sample(allRsids, 5, replace= FALSE)
second_test<-update_data(complete_data, rsidSample)
test_data_5<-second_test %>% mutate(NumSnps=5)
test_data_5<-test_data_5 %>% mutate(y=as.integer(isParent))
test_data_5$y<-as.factor(test_data_5$y)

rsidSample<- sample(allRsids, 40, replace= FALSE)
second_test<-update_data(complete_data, rsidSample)
test_data_40<-second_test %>% mutate(NumSnps=40)
test_data_40<-test_data_40 %>% mutate(y=as.integer(isParent))
test_data_40$y<-as.factor(test_data_40$y)

rsidSample<- sample(allRsids, 500, replace= FALSE)
second_test<-update_data(complete_data, rsidSample)
test_data_500<-second_test %>% mutate(NumSnps=500)
test_data_500<-test_data_500 %>% mutate(y=as.integer(isParent))
test_data_500$y<-as.factor(test_data_500$y)

test_data_all <- rbind(test_data_5, test_data_40, test_data_500) %>% mutate(Proportion_Diff_0=diff_0/NumSnps, Proportion_Diff_2=diff_2/NumSnps)

#This code plots the densities of identical SNPs for the three sample sizes
test_data_all %>% mutate(Relationship_Type = ifelse(y == 1, "Parental", "Non_Parental")) %>% ggplot(aes(Proportion_Diff_0, fill=Relationship_Type, position="identity", color=Relationship_Type)) + geom_density(alpha=0.2, position="identity", adjust=2.0) + facet_grid(NumSnps ~ .) + ggtitle("Proportion of Identical Snps for Different Sample Sizes")

#This code re-plots the densities of identical SNPs for the three sample sizes to overcome intermittent bug
test_data_all %>% mutate(Relationship_Type = ifelse(y == 1, "Parental", "Non_Parental")) %>% ggplot(aes(Proportion_Diff_0, fill=Relationship_Type, position="identity", color=Relationship_Type)) + geom_density(alpha=0.2, position="identity", adjust=2.0) + facet_grid(NumSnps ~ .) + ggtitle("Proportion of Identical Snps for Different Sample Sizes")


#This code plots the densities of wholly different SNPs for the three sample sizes
test_data_all %>% mutate(Relationship_Type = ifelse(y == 1, "Parental", "Non_Parental")) %>% ggplot(aes(diff_2, fill=Relationship_Type, position="identity", color=Relationship_Type)) + geom_density(alpha=0.2, position="identity", adjust=2.5) + facet_grid(NumSnps ~ .) + ggtitle("Number of Mendelian Errors for Different Sample Sizes")



#Only the 10 parent child pairs are included in the heatpmap
distance_people<-c(pedigree$person_1[1:10], pedigree$person_2[1:10])
distance_data<- large_dat %>% filter(person %in% distance_people) %>% mutate(family=0, g_count=0) %>% filter(str_count(genotype, "C")==0) %>% filter(str_count(genotype, "T")==0)


for (i in 1:nrow(distance_data)){
  distance_data[i,7]=get_family(distance_data[i,2])
  distance_data[i,8]=str_count(distance_data[i,6], "G") #Distance is be based on the count of "G" values
}

distance_wide <- distance_data %>% select(person, family, rsid, g_count) %>% spread(rsid, g_count)
a<-(colSums(is.na(distance_wide))==0) #Excludes SNP values that are not available for everyone, likely because they include "AC" or similar for at least one person
b<-distance_wide[,a]
c<-dist(b)
d<-as.matrix(c)
heatmap(d, Rowv=NA, Colv=NA, keep.dendro=FALSE, labRow=b$person, labCol=b$person, main="Distance by Guanine Count")

#This code determines and plots principal components using the same distance measure as above
pca<-prcomp(as.matrix(b[,2:2158]))
e<-data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], family = as.vector(distance_wide$family), label=factor(as.vector(distance_wide$family)))
e %>% ggplot(aes(PC1, PC2, fill=label)) + geom_point(cex=3, pch=21) + ggtitle("Principal Components Labelled by Family")



NumSnps<-seq(5, 200, by=5)
results<-data.frame("glm_sens"=1:length(NumSnps), "glm_spec"=1:length(NumSnps), "knn_sens"=1:length(NumSnps), "knn_spec"=1:length(NumSnps), "rf_sens"=1:length(NumSnps), "rf_spec"=1:length(NumSnps), "rpart_sens"=1:length(NumSnps), "rpart_spec"=1:length(NumSnps))

#The commented code below will run the main simulation of the project, but it takes up to 24 hours, so the alternative code below will download the results from Github
#for (i in 1:length(NumSnps)) {
#  temp_glm_sens <- 1:200
#  temp_glm_spec <- 1:200
#  temp_knn_sens <- 1:200
#  temp_knn_spec <- 1:200
#  temp_rf_sens <- 1:200
#  temp_rf_spec <- 1:200
#  temp_rpart_sens <- 1:200
#  temp_rpart_spec <- 1:200
#  print("i=")
#  print(i)
#  for (j in 1:200) {
#    print("  j= ")
#    print(j)
#    rsidSample<- sample(allRsids, NumSnps[i], replace= FALSE)
#    first_test<-update_data(complete_data, rsidSample)
#    pre_split_data<-first_test %>% mutate(y=as.integer(isParent))
#    a<-pre_split_data %>% filter(isParent==TRUE)
#    b<-pre_split_data %>% filter(isParent==FALSE)
#    test_index_pos<-sample(1:10, 2, replace= FALSE)
#    test_index_neg<-sample(1:290, 58, replace= FALSE)
#    train_index_pos <- setdiff(1:10, test_index_pos)
#    train_index_neg <- setdiff(1:290, test_index_neg)
#    test_data<-rbind(a[test_index_pos,], b[test_index_neg,])
#    train_data<-rbind(a[train_index_pos,], b[train_index_neg,])
#    test_data$y<-as.factor(test_data$y)
#    train_data$y<-as.factor(train_data$y)
#    train_glm<-train(y~diff_0+diff_2, method="glm", data=train_data)
#    y_hat_glm <- predict(train_glm, test_data, type="raw")
#    a<- confusionMatrix(y_hat_glm, test_data$y, positive="1")
#    train_knn<-train(y~diff_0+diff_2, method="knn", data=train_data)
#    y_hat_knn <- predict(train_knn, test_data)
#    b<- confusionMatrix(y_hat_knn, test_data$y, positive="1")
#    train_rf<-train(y~diff_0+diff_2, method="rf", data=train_data)
#    y_hat_rf <- predict(train_rf, test_data)
#    c<- confusionMatrix(y_hat_rf, test_data$y, positive="1")
#    train_rpart<-train(y~diff_0+diff_2, method="rpart", data=train_data)
#    y_hat_rpart <- predict(train_rpart, test_data)
#    d<- confusionMatrix(y_hat_rpart, test_data$y, positive="1")
#    temp_glm_sens[j] <- a$byClass[1]
#    temp_glm_spec[j] <- a$byClass[2]
#    temp_knn_sens[j] <- b$byClass[1]
#    temp_knn_spec[j] <- b$byClass[2]
#    temp_rf_sens[j] <- c$byClass[1]
#    temp_rf_spec[j] <- c$byClass[2]
#    temp_rpart_sens[j] <- d$byClass[1]
#    temp_rpart_spec[j] <- d$byClass[2]
#  }
#  results[i,1]<-mean(temp_glm_sens)
#  results[i,2]<-mean(temp_glm_spec)
#  results[i,3]<-mean(temp_knn_sens)
#  results[i,4]<-mean(temp_knn_spec)
#  results[i,5]<-mean(temp_rf_sens)
#  results[i,6]<-mean(temp_rf_spec)
#  results[i,7]<-mean(temp_rpart_sens)
#  results[i,8]<-mean(temp_rpart_spec)
#  print(results[i,])
#}

#write.csv(results, 'Saxon_Capstone_saved_results')

results<-read.csv('Saxon_Capstone_saved_results')
NumSnps<-seq(5, 200, by=5)
results2<-results %>% mutate(NumSnps=NumSnps, glm_balacc=(glm_sens+glm_spec)/2, knn_balacc=(knn_sens+knn_spec)/2, rf_balacc=(rf_sens+rf_spec)/2, rpart_balacc=(rpart_sens+rpart_spec)/2)


# Plots the sensitivity and specificity obtained for the glm model for each sample size
final_results <- results2 %>% gather(MeasType, Meas, glm_sens, glm_spec) %>% select(NumSnps, MeasType, Meas)
final_results %>% ggplot(aes(NumSnps, Meas, col=MeasType)) + geom_line() + ggtitle("Sensitivity and Specificity for GLM Model")

# Plots the balanced accuracy of each of the four models
final_results <- results2 %>% gather(MeasType, Meas, glm_balacc, knn_balacc, rf_balacc, rpart_balacc) %>% select(NumSnps, MeasType, Meas)
final_results %>% ggplot(aes(NumSnps, Meas, col=MeasType)) + geom_line() + ggtitle("Balanced Accuracy for Four Models")


# Re-runs the knn model over a broader tuning grid, with often inconsistent results
rsidSample<- sample(allRsids, 50, replace= FALSE)
first_test<-update_data(complete_data, rsidSample)
pre_split_data<-first_test %>% mutate(y=as.integer(isParent))
a<-pre_split_data %>% filter(isParent==TRUE)
b<-pre_split_data %>% filter(isParent==FALSE)
test_index_pos<-sample(1:10, 2, replace= FALSE)
test_index_neg<-sample(1:290, 58, replace= FALSE)
train_index_pos <- setdiff(1:10, test_index_pos)
train_index_neg <- setdiff(1:290, test_index_neg)
test_data<-rbind(a[test_index_pos,], b[test_index_neg,])
train_data<-rbind(a[train_index_pos,], b[train_index_neg,])
test_data$y<-as.factor(test_data$y)
train_data$y<-as.factor(train_data$y)
train_knn<-train(y~diff_0+diff_2, method="knn", data=train_data, tuneGrid = data.frame(k = seq(1, 10, 1)))
y_hat_knn <- predict(train_knn, test_data)
b<- confusionMatrix(y_hat_knn, test_data$y, positive="1")

ggplot(train_knn) + ggtitle("Knn Tuning Results for 50 RSID Sample Size")



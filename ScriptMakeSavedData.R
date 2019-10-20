library(tidyverse)
setwd("/home/jeffrey/R/CapstonePaternity2")
pedigree<-read_tsv("Pedigree.csv")
permitted_gts<- c("AA", "AC", "AG", "AT", "CC", "CG", "CT", "GG", "TT")
individuals<-c(pedigree$person_1[1:15], pedigree$person_2[1:10], pedigree$person_2[15])


filename<-paste("./raw23/", individuals[1], ".txt", sep="")
a<-read_tsv(filename)
b<-a %>% filter(chromosome>2 & chromosome < 20)
c<-b %>% filter(genotype %in% permitted_gts)
ranSamp<-sample(c$rsid, 10000, replace= FALSE)
d<- c %>% filter(rsid %in% ranSamp)
e<-d %>% mutate(person=individuals[1])
dat<-select(e, person, rsid, chromosome, position, genotype)

build_data<- function(x){
  result<-x
  for (i in 2:length(individuals)) {
    filename<-paste("./raw23/", individuals[i], ".txt", sep="")
    a<-read_tsv(filename)
    b<-a %>% filter(chromosome>2 & chromosome < 20)
    c<-b %>% filter(genotype %in% permitted_gts)
    d<- c %>% filter(rsid %in% ranSamp)
    e<-d %>% mutate(person=individuals[i])
    dat2<-select(e, person, rsid, chromosome, position, genotype)
    result<-rbind(result,dat2)
  }
  result
}

find_common_snps <- function (x) {
  a <- unique(x$person)
  b <- x %>% filter(person==a[1])
  result <- b$rsid
  for (i in 2:length(a)) {
    c <- x %>% filter(person==a[i])
    result <- intersect(result, c$rsid)
  }
  result
}

CommonSnps <- find_common_snps(large_dat)
full_snps<-intersect(CommonSnps, ranSamp)
full_data<-large_dat %>% filter(rsid %in% full_snps)
final_list_snps<- sample(unique(full_data$rsid), 5000, replace= FALSE)
saved_data <- full_data %>% filter(rsid %in% final_list_snps)
write.csv(saved_data, 'saved_data.csv')

setwd("/Users/yang/Desktop/Disertation_2/data")

data1 = read.table("CryptoWakeUpSampleSheetPlusTags.txt" , header = TRUE)

data2 = read.table("CW-kallisto-abundance-foldchange-long-bygene.txt" , header = F)

data3 = read.table("H99_allorfs_promoter500nt_5mercounts.txt",header = T)



library('glmnet')
library('nnet')
library('permutations')
library('ggplot2')
library('tidyverse')
library('factoextra')

############################TEST SECTION############################
data2_sub = data2
data2_sub$replicates = substr(data2$V1,4,4)

p<-ggplot(data=data2_sub, aes(x=V1,y=log(V5)))+geom_boxplot(aes(fill=replicates))
p +theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) + 
  xlab("Samples") + ylab("Log(fold change)") 
#ggtitle("Comparison of the log(fold change) of two replicates") 

#set the threshold, we say the estimated gene count is 0 or 1 is useless in our analysis

data2_threshold_A = data2_sub[which(data2_sub[,3]>=1 & data2_sub$replicates == 'A'),] #discard those genes whose estimated counts <= 1
data2_threshold_B = data2_sub[which(data2_sub[,3]>=1 & data2_sub$replicates == 'B'),]


#clustering the rep A
q1_A = data2_threshold_A[,c(1,2,5)]
names(q1_A) [c(1,2,3)]=c("code","gene","fold_change")

q1_A$fold_change = log(q1_A$fold_change)#take log

q1_tidy_A= pivot_wider(q1_A, id_cols = gene, names_from = code, values_from = fold_change)#transformation
#After tidying up the data, rows are genes and columns are each samples from various conditions.
#The numbers are log(fold change)

q1_cluster_A = q1_tidy_A[,2:19]
q1_cluster_A = na.omit(q1_cluster_A) #remove NAs


#kmeans
set.seed(1234)
fviz_nbclust(q1_cluster_A, kmeans, method = "silhouette") 

q1_kmeans<-kmeans(q1_cluster_A,2)
fviz_cluster(q1_kmeans, q1_cluster_A) 

#clustering the rep B
q1_B = data2_threshold_B[,c(1,2,5)]
names(q1_B) [c(1,2,3)]=c("code","gene","fold_change")

q1_B$fold_change = log(q1_B$fold_change)#take log

q1_tidy_B= pivot_wider(q1_B, id_cols = gene, names_from = code, values_from = fold_change)#transformation
#After tidying up the data, rows are genes and columns are each samples from various conditions.
#The numbers are log(fold change)

q1_cluster_B = q1_tidy_B[,2:19]
q1_cluster_B = na.omit(q1_cluster_B) #remove NAs

#kmeans
set.seed(1234)
fviz_nbclust(q1_cluster_B, kmeans, method = "silhouette") 

q1_kmeans<-kmeans(q1_cluster_B,2)
fviz_cluster(q1_kmeans, q1_cluster_B)

############################TEST DONE############################




###################start###################
################preprossess################

#combine the two replicates
data2_sub = data2
data2_sub$replicates = substr(data2$V1,4,4)

data2_sub_A = data2_sub[which(data2_sub$replicates == 'A'),]
data2_sub_B = data2_sub[which(data2_sub$replicates == 'B'),]

#to check if the subset A and B are identical in order
check = numeric(146736)
for(i in 1:146736){
  if(data2_sub_A$V2[i] == data2_sub_B$V2[i]){
    check[i] = T
  }
}
sum(check) #sum is 146736, so the orders are identical, so we can simply average then

data2_sub_A$V5[is.na(data2_sub_A$V5)]=0
data2_sub_B$V5[is.na(data2_sub_B$V5)]=0


gene_count = round((data2_sub_A$V3+data2_sub_B$V3)/2) # use this to remove the row with relatively small counts
fold_change = as.numeric((data2_sub_A$V5+data2_sub_B$V5)/2)#the combined fold change
code =  substr(data2_sub_A$V1,1,3)#sample code without A and B replicates 
gene = as.numeric(data2_sub_A$V2)#gene code with specific order

#obtain a new dataframe with two replicates combined
kali_foldchange = as.data.frame(code)
kali_foldchange$gene = data2_sub_A$V2
kali_foldchange$log_fold_change = log(fold_change)

kali_foldchange = kali_foldchange[which(gene_count>1),]

#tidy the data
kali_foldchange_tidy= pivot_wider(kali_foldchange, id_cols = gene, names_from = code, values_from = log_fold_change)#transformation
#After tidying up the data, rows are genes and columns are each samples from various conditions.
#The numbers are log(fold change)
kali_foldchange_tidy = na.omit(kali_foldchange_tidy)  #remove NAs

#Normality test for the observations
#######normality test####
p_value = numeric(nrow(kali_foldchange_tidy))

for (i in 1:nrow(kali_foldchange_tidy)) {
  test = as.vector(unlist(kali_foldchange_tidy[i,-1]))
  p_value[i] = shapiro.test(test)$p.value
}

is_significant = rep(0,nrow(kali_foldchange_tidy))
is_significant[which(p_value>0.05)] = 1
sum(is_significant == 1)/nrow(kali_foldchange_tidy)

norm_test = as.data.frame(cbind(c(1:length(p_value)),p_value,is_significant))

ggplot(norm_test, aes(x = V1,y=p_value)) +
  geom_point(size=1,aes(color=is_significant)) + 
  xlab("Index") + ylab("P-value for the Shapiro-Wilk test") +theme(legend.title = element_blank()) 


############################build a new data3 from the tidy data#############################
#now we only have a subset from all the genes,do the same to data3
col_in_row = data3$Gene %in% kali_foldchange_tidy$gene
data3_sub = data.frame(data3,col_in_row)
data3_sub = data3_sub[col_in_row,-1026]

#we find the rows of data3_sub is less than the tidy data,so there are some genes included in tidy data
#but not in the data3_sub
#Again

col_in_row_2 = kali_foldchange_tidy$gene %in% data3_sub$Gene
kali_foldchange_tidy = data.frame(kali_foldchange_tidy,col_in_row_2)
kali_foldchange_tidy = kali_foldchange_tidy[col_in_row_2,-20]

kali_foldchange_cluster = kali_foldchange_tidy[,2:19]

#check the optimal cluster num
set.seed(1234)
fviz_nbclust(kali_foldchange_cluster, kmeans, method = "silhouette")  # number of clusters 
fviz_nbclust(kali_foldchange_cluster, kmeans, method = "wss") 

kali_foldchange_kmeans<-kmeans(kali_foldchange_cluster,2)
fviz_cluster(kali_foldchange_kmeans, kali_foldchange_cluster)

#however, lets see if we can imporve the performance

#do the PCA 
PC<-princomp(kali_foldchange_cluster,cor=TRUE) 
summary(PC,loadings=TRUE)

screeplot(PC,type = 'lines',npcs = 18)#from the fourth one, the variance change is slight
#pick first ten


pc = PC$loadings[,c(1:6)]

kali_foldchange_cluster_pca = as.matrix(kali_foldchange_cluster) %*% as.matrix(pc)#make the new data

set.seed(1234)

fviz_nbclust(kali_foldchange_cluster_pca, kmeans, method = "silhouette")  # number of clusters 
fviz_nbclust(kali_foldchange_cluster_pca, kmeans, method = "wss")

kali_foldchange_kmeans_pca<-kmeans(kali_foldchange_cluster_pca,2)
fviz_cluster(kali_foldchange_kmeans_pca, kali_foldchange_cluster_pca)

#calculate the mean of log(fold change) in each row for two clusters

kali_foldchange_tidy$cluster = kali_foldchange_kmeans_pca$cluster
mean_c_1 = mean_c_2 = numeric(18)
for (i in 1:18) {
  mean_c_1[i] = mean(kali_foldchange_tidy[,i+1][which(kali_foldchange_tidy$cluster == 1)])
  mean_c_2[i] = mean(kali_foldchange_tidy[,i+1][which(kali_foldchange_tidy$cluster == 2)])
}

mean(mean_c_1)
mean(mean_c_2)
#clsuter one is higher in gene expression amount




######################3deal with the 5-mers data######################

#function to transform ATGC
trans = function(a){
  if(a == 'A'){
    return('T')
  }else if(a == 'T'){
    return('A')
  }else if(a == 'C'){
    return('G')
  }else{return('C')}
}

gene_cluster = as.data.frame(kali_foldchange_tidy[,1])
names(gene_cluster)[1] = 'Gene'

gene_cluster$cluster = kali_foldchange_kmeans_pca$cluster
gene_cluster = cbind(gene_cluster,class.ind(gene_cluster$cluster))

#This is the function to combine the motifs according to the k 
combine_the_pair = function(data,k){
  motifs_num = 4^k
  data_combined = data.frame(matrix(0,nrow(gene_cluster),motifs_num/2+1))
  data_combined[,1] = data[gene_cluster$Gene,1]
  names(data_combined)[1] = 'Gene'
  name = colnames(data)
  #loop to combine the pair
  for (i in 2:(motifs_num/2+1)) {
    str = strsplit(name[i],split = "")
    temp = character(k)
    for (j in k:1) {
      temp[j]=trans(unlist(str)[abs(j-k-1)])
    }
    temp = paste(temp,collapse = '')
    num_order = which(name == temp)
    data_combined[,i] = (data[gene_cluster$Gene,i]+data[gene_cluster$Gene,num_order])/2
    names(data_combined)[i] = paste(name[i],temp,collapse = ',')
  }
  return(data_combined)
}

data3_combined = combine_the_pair(data3,5)

#q3
set.seed(4321)
x = as.matrix(data3_combined[,c(2:513)])


#also we apply the pca to the tidied data
kali_foldchange_tidy_q3 <- kali_foldchange_tidy[!is.infinite(rowSums(kali_foldchange_tidy[,2:19])),2:19]

PC1<-princomp(kali_foldchange_tidy_q3,cor=TRUE) 
summary(PC1,loadings=TRUE)

screeplot(PC1,type = 'lines',npcs = 18)#from the fourth one, the variance change is slight

y_pca = as.matrix(kali_foldchange_tidy[,2:19]) %*% PC1$loadings[,c(1:6)]#pick the first 6 components


#pre set the lambda value and test the amount significant features for each model
lambda_test = seq(0,0.1,0.002)
num_significant_motifs = numeric(51)

for (i in 1:51) {
  lasso_test = glmnet(x = x,y = y_pca,
                      lambda = lambda_test[i], family = 'mgaussian')
  fit_coefficient_test = coef(lasso_test)
  num_significant_motifs[i] = length(colnames(x)[which(unlist(fit_coefficient_test[[1]]!=0))])
}
#combine the output to draw a plot
num = as.data.frame(cbind(lambda_test,num_significant_motifs))

#draw the picture to show the relation between num of significant features and lambda
ggplot(num, aes(x = lambda_test,y=num_significant_motifs)) +
  geom_point(size=1) + 
  xlab("Lambda") + ylab("Number of significant motifs") +geom_hline(yintercept = c(10))

#fit the model with 10 significant features
lasso_fit = glmnet(x = x,y = y_pca,
                   lambda = 0.16, family = 'mgaussian')



fit_coefficient = coef(lasso_fit)
colnames(x)[which(unlist(fit_coefficient[[1]]!=0))]
fit_coefficient[[1]][which(unlist(fit_coefficient[[1]]!=0))]

#significant features and its corresponding coef
feature_coef = as.data.frame(cbind(colnames(x)[which(unlist(fit_coefficient[[1]]!=0))],
                                   fit_coefficient[[1]][which(unlist(fit_coefficient[[1]]!=0))],
                                   abs(fit_coefficient[[1]][which(unlist(fit_coefficient[[1]]!=0))])))


#import the 4_mers and 6_mers data
six_mers = read.table("H99_allorfs_promoter500nt_6mercounts.txt",header = T)
four_mers = read.table("H99_allorfs_promoter500nt_4mercounts.txt",header = T)


################################ 6_mer ################################
six_mers_combined = combine_the_pair(six_mers,6)
x_6_mers = as.matrix(six_mers_combined[,c(2:length(six_mers_combined))])

lambda_test_six = seq(0,0.30,0.01)
#same operation to the  6-mers

num_significant_motifs_six = numeric(31)
for (i in 1:31) {
  lasso_test_six = glmnet(x = x_6_mers,y = y_pca,
                      lambda = lambda_test_six[i], family = 'mgaussian')
  fit_coefficient_test_six = coef(lasso_test_six)
  num_significant_motifs_six[i] = length(colnames(x_6_mers)[which(unlist(fit_coefficient_test_six[[1]]!=0))])
}
#combine the output to draw a plot
num_six = as.data.frame(cbind(lambda_test_six,num_significant_motifs_six))

#draw the picture to show the relation between num of significant features and lambda
ggplot(num_six, aes(x = lambda_test_six,y=num_significant_motifs_six)) +
  geom_point(size=1) + 
  xlab("Lambda") + ylab("Number of significant motifs") +geom_hline(yintercept = c(10))

lasso_fit_6 = glmnet(x = x_6_mers,y = y_pca,
                     lambda = 0.18, 
                     family = 'mgaussian')

#significant features and its corresponding coef
fit_coefficient_6 = coef(lasso_fit_6)
feature_coef_6 = as.data.frame(cbind(colnames(x_6_mers)[which(unlist(fit_coefficient_6[[1]]!=0))],
                                     fit_coefficient_6[[1]][which(unlist(fit_coefficient_6[[1]]!=0))],
                               abs(fit_coefficient_6[[1]][which(unlist(fit_coefficient_6[[1]]!=0))])))

################################ 4_mer ################################
four_mers_combined = combine_the_pair(four_mers,4)
x_4_mers = as.matrix(four_mers_combined[,c(2:length(four_mers_combined))])

sum(is.na(x_4_mers))

#same operation to the  4-mers
lambda_test_four = seq(0,0.3,0.01)
num_significant_motifs_four = numeric(31)
for (i in 1:31) {
  lasso_test_four = glmnet(x = x_4_mers,y = y_pca,
                          lambda = lambda_test_four[i], family = 'mgaussian')
  fit_coefficient_test_four = coef(lasso_test_four)
  num_significant_motifs_four[i] = length(colnames(x_4_mers)[which(unlist(fit_coefficient_test_four[[1]]!=0))])
}
#combine the output to draw a plot
num_four = as.data.frame(cbind(lambda_test_four,num_significant_motifs_four))

#draw the picture to show the relation between num of significant features and lambda
ggplot(num_four, aes(x = lambda_test_four,y=num_significant_motifs_four)) +
  geom_point(size=1) + 
  xlab("Lambda") + ylab("Number of significant motifs") +geom_hline(yintercept = c(10))

lasso_fit_4 = glmnet(x = x_4_mers,y = y_pca,
                     lambda = 0.14, family = 'mgaussian')
#significant features and its corresponding coef
fit_coefficient_4 = coef(lasso_fit_4)
feature_coef_4 = as.data.frame(cbind(colnames(x_4_mers)[which(unlist(fit_coefficient_4[[1]]!=0))],
                                   fit_coefficient_4[[1]][which(unlist(fit_coefficient_4[[1]]!=0))],
                                   abs(fit_coefficient_4[[1]][which(unlist(fit_coefficient_4[[1]]!=0))])))




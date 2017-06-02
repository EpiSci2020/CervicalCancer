
#This is an analysis of data from the paper:
#RNA discovery and discrete statistical biomarker analysis 
#in a collection of cervical tumours and matched controls
#BMC Biology 2010, 8:58

#I have no domain expertise in gene expression technology.
#This is my attempt at a general exploratory analysis of 
#the miRNA data used in the Witten et al. study.

#### Preliminaries

library(caret)
library(e1071)
library(tidyverse)
library(MASS)
library(stringr)

# Set working directory
setwd("C:/Users/Dave/Desktop/Kaggle/CervicalCancer")
getwd()
set.seed(123)

#*** miRNA dictionary  ***********
#https://docs.google.com/spreadsheets/d/1eKXHuNRyOunGhbV3JVG3KVbKmMIYe6x5SBsBwPbyZNk/edit?usp=sharing

cervical <- read.csv("cervical.csv")
str(cervical, max.level = 0)
data <- cervical # shorter name to work with later

col_names <- cervical[1]
str(col_names) # I like to do a lot of structure checks
col_names1 <- as.matrix(col_names)
#Substitute letters for the *'s and-'s
col_names2 <- as.matrix(str_replace_all(col_names1, pattern = fixed("*"), "x"))
col_names3 <- as.matrix(str_replace_all(col_names2, pattern = fixed("-"), "z"))
# Change the matrix to a data frame for later use
nams <- data.frame(matrix(nrow = 1, ncol = 714))
colnames(nams) <-col_names3 #Installs the column names


# Read in the data file
#data <- read_csv('cervical_ok.csv')
#data <- data[-1,]    # remove dx (diagnosis) row

# Remove the numerics from the index and dependent variable (dx)
data_nums <- data[-1]
# Put the numerics in a data matrix (not a matrix, need
#to preserve columnar structure)
data_nums_mat <- data.matrix(data_nums)

#Transpose the numerics
data_nums_mat_tr <- t(data_nums_mat)

# Convert the data matrix to a data frame to hold non-numeric data
data_nums_mat_tr_df <- data.frame(data_nums_mat_tr)
colnames(data_nums_mat_tr_df) <- col_names3

# Center and scale the numerics  - gives a matrix
data_nums_mat_tr_df_sc <- scale(data_nums_mat_tr_df) 

#Find low variance columns to remove
nzv <- nearZeroVar(data_nums_mat_tr_df_sc, saveMetrics = TRUE)
print(paste("Range: ", range(nzv$percentUnique)))

# Get the indices of the low-variance columns
inNums <- as.matrix(which(nzv$nzv, arr.ind = FALSE, useNames = TRUE))

# Keep the non-low variance columns
nums <- data_nums_mat_tr_df_sc[,-inNums]

# Re-convert numerics to a data frame
nums_df <- data.frame(nums)
dim(nums_df)



# Create a data frame of diagnosis labels for each subject
dx <- data.frame(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))
colnames(dx) <- "dx"
dim(dx)

# *** We now have a scaled data frame of numbers - nums_df and
#                 a data frame of diagnosis codes - dx

#  **** combine the dx and nums   ******
cancer <- cbind(nums_df, dx)
dim(cancer)



#Change the data from wide to long form for later use
cancer1 <- cbind(dx,data_nums_mat_tr_df)
dim(cancer1)
cancer_melt <- gather(cancer1,miRNA,counts,2:715)
dim(cancer_melt)
cancer_melt$dx <- factor(cancer_melt$dx, labels = c("Normal", "Tumor")) # factorize dx
sum(is.na(cancer_melt$counts)) # check for NA's


# ***** Data Exploration  *******************************
#Similar variance between diagnosis groups?
std_dev <- cancer_melt %>%  group_by(dx) %>%
      summarise(std_dev = sd(counts))
round(std_dev$std_dev,0)

#Histogram of log of counts in each diagnosis group
cancer_p <- ggplot(cancer_melt, aes(x = log(counts))) +
   geom_histogram(binwidth = .5) +
   facet_wrap(~dx)
cancer_p

# The samples are independent, of equal number,
# and not normally distributed
#Mann-Whitney U-test for independent, non-normal geoups
U_test <- wilcox.test(cancer_melt$counts ~ cancer_melt$dx)
U_test$p.value
#Indicates a significant differnce between count groups

   
# ***  Heatmap of miRNA data  *******
library(heatmaply)
# Selecting the highest variance miRNA's
y <- heatmaply::normalize(data_nums_mat_tr_df) %>% select_if(function(col) var(col) >.06)
heatmap(data.matrix(y), Rowv=NA, Colv=NA, col = heat.colors(256), scale="column", margins=c(5,10))


# For the miRNA data, the residuals of a linear model are not normally distributed.
# The count means and variances are greatly dissimilar.
# Therefore, neither a linear OLS model or a poisson model is appropriate, and a 
# negative binomial model is chosen.

# ** Negative Binomial Model  *******

# The cube root of the counts was used in the model to reduce the skewness of the data
ca_glmnb <- glm.nb(counts^(1/3) ~ dx, data = cancer_melt)
summary(ca_glmnb)

exp(coef(ca_glmnb))


# *** t-Tests for individual miRNA's by dx *******


a = 1
b = 28
c = 29
d = 58
df = cancer_melt

take <- function(a,b,c,d, df) {
   df = cancer_melt
   taken1 <- df[a:b,3]
   taken2 <- df[c:d,3]
   tt <- t.test(taken1,taken2)
   mirna <- ifelse(tt$p.value < 0.05,cancer_melt[[a,2]], '')
   print(tt$p.value)
   print(mirna)
}
for (i in 1:714) {
    take(a,b,c,d,df)
    print(i)
    a = a + 58
    b = b + 58
    c = c + 58
    d = d + 58
}



#***** Prediction  *************************************

#***********  PCA Model **********************************************************


#The principal components are the linear combinations of the
#miRNAs that have the largest variance and provide alternative 
#independent variables for the predictive model

# Create the PCA model (caret)
pca <- preProcess(x = nums, method = 'pca', pcaComp = 3)
nums_pca = data.frame(predict(pca, nums))
# the data frame of principal components

#Scree Plot 
vpc1 <- var(nums_pca$PC1)
vpc2 <- var(nums_pca$PC2)   
vpc3 <- var(nums_pca$PC3)
vars <- c(vpc1,vpc2,vpc3)
index <- c('PC1','PC2','PC3')
vars_df <- data.frame(index,vars)

scree <- ggplot(vars_df, aes(x = index, y = vars)) + geom_bar(stat='identity')
dev.off()
scree


# ******** PREDICTION - PCA model  ******************

cancer <- cbind(nums_pca,dx)
colnames(cancer) <- c("PC1", "PC2", "PC3", "dx")

cancer$dx <- factor(cancer$dx) 

# Split the data for prediction (from caret package)
inTrain <- createDataPartition(cancer$dx, p = 0.75, list = FALSE)
train <- cancer[inTrain,]
test <- cancer[-inTrain,]

# Using a support vector machine algorithm
model.svm = svm(formula = dx ~ .,
                 data = cancer,
                 type = 'C-classification',
                 kernel = 'linear')

# Predicting the Test set results
y_pred = predict(model.svm, newdata = test)

# Making the Confusion Matrix
cm <- confusionMatrix(y_pred, test$dx)
cm



# ****** XGBoost *******************
library(xgboost)
# data_nums_mat_tr_df
# xgboost does not require scaled data
#colnames(ca_xgb[1]) <- as.matrix('dx')
dx_num <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
class(dx_num)
dx_chr <- ifelse(dx_num == 0,"normal", "tumor")
dx_fac <- factor(dx_chr,levels = c("normal", "tumor"))
str(dx_fac)
str(data_nums_mat_tr_df, max.level = 0)
data_nums_mat_tr_df$dx <- dx_fac  
#data_nums_mat_tr_df$dx <- factor(as.character(data_nums_mat_tr_df$dx), levels = c("normal", "tumor"))
inTrain <- createDataPartition(data_nums_mat_tr_df$dx, p = 0.75, list = FALSE)
train <- data_nums_mat_tr_df[inTrain,]
#train$dx <- dx  # add dx column to train
test <- data_nums_mat_tr_df[-inTrain,] 
#train <- as.matrix(train)
dim(train)

xgbGrid <- expand.grid(nrounds = c(1, 10),
                       max_depth = c(1, 4),
                       eta = c(.1, .4),
                       gamma = 0,
                       colsample_bytree = .7,
                       min_child_weight = 1,
                       subsample = c(.8, 1))

ctrl1 <- trainControl(method = "LOOCV",
                      classProbs = TRUE,
                      summaryFunction = twoClassSummary)
m_xgb <- train(dx ~ ., data = train, 
                     method = "xgbTree", 
                     trControl = ctrl1,
                     metric = "ROC", 
                     #preProc = c("center", "scale")
                     tuneGrid = xgbGrid)

y_pred = data.frame(predict(m_xgb, newdata = data.matrix(test)))
dx <- data.frame(unlist(test$X1))
# Making the Confusion Matrix
test_dx_df <- data_frame(test$dx)
cm <- confusionMatrix(y_pred, test_dx_df)
cm
#*************************************************


eval <- function(x) {
  a <- is.matrix(x)
  b <- is.atomic(x)
  c <- is.array(x)
  d <- is.character(x)
  e <- is.data.frame(x)
  f <- is.vector(x)
  g <- is.factor(x)
  print(c("matrix    ",a));#  print(a)
  print(c("atomic    ",b));#  print(b)
  print(c("array     ",c)); # print(c)
  print(c("character ",d)); # print(d)
  print(c("data frame",e)); # print(e)
  print(c("vector    ",f)); # print(f)
  print(c("factor    ",g)); # print(g)
}



















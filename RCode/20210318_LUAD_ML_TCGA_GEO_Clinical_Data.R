


#############################################################################
## 
## Feature selection with SVM-RFE
## 04/14/2019 - 01/17/2020
## Jeff Du
##  
#############################################################################
## 
## use caret package to build classfication model
## use the DESeq2 generated significant list as pre-selected genes
## use pROC to plot ROC for training set
## 
#############################################################################
## 
#############################################################################
# install necessar packages
# library(devtools)
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("shiny")
# install_github("StatsWithR/statsr")
# BiocManager::install("sigFeature", version = "3.8")
# install.packages("mlbench")
# install.packages("caret")
# install.packages("ellipse")   ## for feature plot; 
#############################################################################

### load packages
### library(mlbench)
library(caret)
library(ggplot2)
# library(ellipse) ## one of the dependencies for pROC package
library(pROC)
library(doParallel)


set.seed(12345)

####################################################
## Data manipulation and pre-processing
## #1 deal with missing/NA data
## #2 remove constant rows
## #3 optimize training set and testing set
####################################################


########################################################################################################
## data of RNAseq for each tumor type from TCGA
## count <- read.table("D:/WorkRecord/Companies/SH_Hospital/201904/SVM_RFE_Prj/svm_rfe_prj/THCA_RNA_Set.RnaSeq_Transcript.Genes.txt",
##                    header = TRUE,
##                    row.names = 1,
##                    sep = "\t"
##                    )
########################################################################################################
## design <- read.table("D:/WorkRecord/Companies/SH_Hospital/201904/SVM_RFE_Prj/svm_rfe_prj/THCA_RNA_Set.RnaSeq_Transcript.Genes_Design.txt",
##                     header = TRUE,
##                     row.names = 1,
##                     sep = "\t"
##                     )
########################################################################################################


########################################################################################
## 
## mRNA data 
## 
## #1 input a 49 * 43 matrix, which is the count table of a mRNA dataset
## originally, there are 784 genes, after deseq2 One-Way-Test, only 49 of them showed
## slighly diffrencial expression;
## 
## #2 a design talbe, within the table there are several columns for survival plot
## and the classification ml model works on progression column;
## 
## #3 the goal is to trains a model to predict progrsssion Yes/No
## 
########################################################################################



########################################################
## Section 0: save ROC plot
## 
## function plot ROC cure and save as jpg file
## 
## pass an roc_object together with routine and file name
## Plot ROC curve, add AUC, 95% CI, save to target folder
## Also, pass the ML algorithm to the legend;

PlotROCSaveJPG <- function(roc.obj, path, file.name, ml.algo){
  
  jpeg(paste0(path, file.name, "_ROC_with_AUC.jpg"), 
       width = 8, height = 6, units = "in", res= 400)
  
  ## We can calculate the area under the curve...
  
  plot( 
    roc.obj,
    legacy.axes=TRUE, percent=TRUE, 
    main= ml.algo,
    xlab="False Positive Percentage", 
    ylab="True Postive Percentage", 
    col="darkblue", lwd=4, 
    print.auc=TRUE
  )  # end plotting
  
  # the 95% CI could be printed directly by print.auc=T; 
  #ci.95 <- paste0("95% CI = ", round(ci.auc(rocobj)[1], digits = 2), " - ", round(ci.auc(rocobj)[3], digits = 2) ) 
  
  legend("bottomright", 
         legend=c( paste("model:", ml.algo) ), 
         col=c("darkblue", "darkblue"), 
         lwd=4)
  
  dev.off()
  
} ## End PlotROCSaveJPG() function  
#####################################################################################################################



########################################################################################
## Section I 
## read in the count data and the clinical information table; 
## 
########################################################################################


#########

### get tissue expression data
setwd("C:/Users/Jeff/SH_Hospital/SVM_RFE_Prj/Jeff_Shuai/20200115_input_files/")

getwd()

luad.count <- read.table("./Counts/LUAD_plus_LUSCNormal.MirnaSeq_Count.txt", row.names = 1, header = T, sep = "\t")
luad.design <- read.table("./Counts/LUAD_plus_LUSCNormal.MirnaSeq_Count_Design.txt", row.names = 1, header = T, sep = "\t")
luad.topgenes <- readLines("./Counts/TopList_TCGALUAD_DESeq2.txt")


luad.topgenes[1:10] 
length(luad.topgenes) 

# head(luad.count)
# colnames(luad.count)
luad.count[1:5, 1:5]
luad.design[1:5, ]
class( luad.count$TCGA.05.4244.01A ) 



luad.design$CD274.and.CD8A.HH.LL <- NULL
luad.design[1:5, ]
luad.count[1:5, 1:5]


### fix row names containing '-'
library(stringr)
rownames( luad.count ) <- str_replace_all( rownames(luad.count), c("-" = ".") )
rownames( luad.design ) <- str_replace_all( rownames(luad.design), c("-" = "."))
luad.count[1:5, 1:5]
luad.design[1:5, ] 

#########
### get GSE33858 miRNA expression data; 
## setwd("C:/Users/Jeff/SH_Hospital/SVM_RFE_Prj/0604")

gse.count <- read.table("./Counts/GSE33858Count.Transcript_Count.txt", row.names = 1, header = T, sep = "\t")
gse.design <- read.table("./Counts/GSE33858Count.Transcript_Count_Design.txt", row.names = 1, header = T, sep = "\t")
gse.topgenes <- readLines("./Counts/TopList_GSE33858_DESeq2.txt") 

## fix row names containg '-'

rownames( gse.count ) <- str_replace_all( rownames(gse.count), c("-" = ".") )

gse.count[1:5, 1:5]
gse.design[1:5, ]
# class(gse.count$GSM838063_159T)

gse.topgenes[1:10] 
length(gse.topgenes) 

#########
## clinical data
clinic.count <- read.table("./Counts/ClinicalPatients_0110_miRNA.Transcript_Count.txt", row.names = 1, header = T, sep = "\t")
clinic.design <- read.table("./Counts/ClinicalPatients_0110_miRNA.Transcript_Count_Design.txt", row.names = 1, header = T, sep = "\t")

### fix row names containing '-'
rownames( clinic.count) <- str_replace_all( rownames(clinic.count), c("-" = "."))

clinic.count[1:5, 1:5]
clinic.design[1:5, ]

###############
## in the clinical data, X704 is an extra item; 

dim(clinic.count)
head(clinic.count)


###############
## get 

top.genes <- intersect( luad.topgenes, gse.topgenes) 
length(top.genes)
top.genes

### replace all '-' in the gene names with '.'
top.genes <- str_replace_all( top.genes, c("-" = ".") )
top.genes


##################
## after running the code till Line 500
## we got the top 50 genes with most important contributions to the svmLinear model; 
top.genes <- top.genes.50

### since the top 21 genes have already been saved to a local txt file
### for downstream analysis, we load it to R env directory

top.genes.21 <- readLines("./Counts/top.genes.21.txt") 
top.genes.21

## further iteration confirmed top 21 genes are good 
top.genes <- top.genes.21



## it's interesting there's no overlap between the two top-gene lists; 

# rownames.ori <- rownames(clinical.data, do.NULL = T, prefix = "X")
# 
# rownames.new <- paste0('X', rownames.ori)
# row.names(clinical.data) <- rownames.new

dim(luad.count)
dim(gse.count)
dim(clinic.count)

class( luad.count$TCGA.05.4244.01A)
class( gse.count$GSM838063_159T)
class( clinic.count$X1232746_CA_GTCCGCAT_S44)

############################################################
## remove rows with constant values (zeros)
## for mRNA data, remove rows with varience less than 0.0001
## 

datamatrix.luad <- luad.count
datamatrix.gse <- gse.count
datamatrix.clinic <- clinic.count

class(datamatrix.luad$TCGA.05.4244.01A)


dim(datamatrix.luad)
dim(datamatrix.gse)
dim(datamatrix.clinic)


## 
# constantRow.luad <- row.names(as.matrix(which(apply(datamatrix.luad, MARGIN = 1, function(x) var(x) < 0.001) ) ) ) 
# constantRow.gse <- row.names(as.matrix(which(apply(datamatrix.gse, MARGIN = 1, function(x) var(x) < 0.001) ) ) ) 
# constantRow.clinic <- row.names(as.matrix(which(apply(datamatrix.clinic, MARGIN = 1, function(x) var(x) < 0.001) ) ) ) 

#exclude rows with constant value: 
#luad.count <- datamatrix.luad[!row.names(datamatrix.luad) %in% constantRow.luad, ] 
luad.count <- datamatrix.luad[ row.names(datamatrix.luad) %in% top.genes, ]
class( luad.count$TCGA.05.4249.01A )
dim(luad.count)

#exclude rows with constant value: 
# gse.count <- datamatrix.gse[!row.names(datamatrix.gse) %in% constantRow.gse, ] 
gse.count <- datamatrix.gse[ row.names(datamatrix.gse) %in% top.genes, ]
dim(gse.count)


# clinic.count <- datamatrix.clinic[!row.names(datamatrix.clinic) %in% constantRow.clinic, ] 
clinic.count <- datamatrix.clinic[ row.names(datamatrix.clinic) %in% top.genes, ]
dim(clinic.count)


print("When some variables get constant values, thus got droped during PCA scale:") 
# constantRow.luad
# constantRow.gse
# constantRow.clinic


############################################################
## within the input count table, each row represents a gene
## in the model, it will consider each column as one feature
## thus we need to transfer the per-row-a-gene into per-col-a gene
## use t() fucntion to directly transport the input table
dim(luad.count)
luad_data <- as.data.frame( t(luad.count) )

gse_data <- as.data.frame( t(gse.count))
clinic_data <- as.data.frame( t(clinic.count))

class( luad_data$hsa.miR.200b.3p )

sapply(luad_data, class)

dim(luad_data)
dim(gse_data)
dim(clinic_data)

class(luad_data$hsa.miR.224.5p)
class(gse_data$hsa.miR.200b.3p)
class(clinic_data$hsa.miR.200b.3p)

### convert factor to numeric
number.genes <- dim(luad_data)[2]
number.genes

luad_data2 <- apply(luad_data[, c(1:number.genes)], 2, function(x) as.numeric(as.factor(x)))
luad_data2 <- as.data.frame(luad_data2)

luad_data2[1:5, 1:5]
rownames(luad_data2) <- row.names(luad_data)
class(luad_data2$hsa.miR.200a.5p)

luad_data <- luad_data2

########
### only need to transfer luad_data columns to numeric
any( is.na(luad_data))



luad_data[1:5, 1:5]
gse_data[1:5, 1:5]
clinic_data[1:5, 1:5]

###############################
## merge count data and design data into on dataframe


luad_data$class <- luad.design$SampleType
summary(luad_data$class)

luad_data$class <- ifelse(luad_data$class=="Solid Tissue Normal", "Normal", "Tumor")
summary(luad_data$class)

luad_data$class <- as.factor(luad_data$class)
summary(luad_data$class)

luad_data[1:5, 1:5]

## in case there are multi factors in the class columns, remove them;
## subset the data, only keep primiary tumor and solid normal
#my_data <- subset( my_data, class!='Metastatic') # my_data$class != "Metastatic", ]
#
## check factor levels left
#summary(my_data$class) 

## drop the empty Metastatic level
#my_data$class <- factor(my_data$class)



##############
## for gse and clinic dataset, the factors have pre-defined as Normal vs Tumor
## 
gse_data$class <- gse.design$TumorNormal
summary(gse_data$class)

clinic_data$class <- clinic.design$TumorNormal
summary(clinic_data$class)

####################################################################################
## Section II 
## Create the training and test datasets
####################################################################################


###################################################################
## Section II, Step 1: Get row numbers for the training data
## The limitation in this case is: dataset too small
## So, we split the original data into 75% training + 25% testing
###################################################################



set.seed(123)
trainRowNumbers <- createDataPartition(luad_data$class, p=0.70, list=FALSE)


###################################################################
## Section II, Step 2: Create the training  dataset
trainData <- luad_data[trainRowNumbers,]

###################################################################
## Section II, Step 3: Create the test dataset
testData <- luad_data[-trainRowNumbers,] 

validData <- gse_data
applicationData <- clinic_data

## Section II, Step 3: Create the test dataset by another random sampling; 
## testData <- my_data[testRowNumbers,]


dim(trainData)
dim(testData)

num.col <- dim(trainData)[2]
num.col

# Store X and Y for later use.
x = trainData[, 1:num.col-1]
data.matrix(x)
y = trainData$class 
x[1:5, 1:5]

summary(y)

x_test = testData[ , 1:num.col-1]
data.matrix(x_test)
y_test = testData$class

x_valid = validData[ , 1:num.col-1]
y_valid <- gse_data$class

x_apply = applicationData[ , 1:num.col-1]
y_apply = clinic_data$class

dim(x)
dim(x_test)
dim( x_valid)
dim( x_apply)

x$class
# x$class should be NULL


###############################################################################################
## Section III, data, mamipulation, normalization and scalling
###############################################################################################


########################################################################
## 
## convert all the numeric variables to range between 0 and 1, 
## by setting method=range in preProcess(). 
## preProcess_range_model <- preProcess(trainData, method='range')
## trainData <- predict(preProcess_range_model, newdata = trainData)
## 
## Append the Y variable
## trainData$class <- y
########################################################################


########################################################################
## Section III, step 1 
## remove constant columns
## for top 500/1000 significant genes, this step could be ignored; 

#x <- x[,apply(x, 2, var, na.rm=TRUE) > 0.0001]
#dim(x)


########################################################################
## Section III, step 2
## convert all the numeric variables to range between 0 and 1, by setting method=range in preProcess(). 
x <- data.matrix(x)
x[1:5, 1:5]

trainData_range_model <- preProcess( x, method = c("center", "scale") )

trainData <- predict( trainData_range_model, newdata = x) 
trainData <- as.data.frame(trainData)

trainData$class <- y
# y
head(trainData)
trainData[1:5, 1:5]

########################################################################
## Section III, step 3
## convert all test data 

dim(x_test)

testData_range_model <- preProcess( x_test, method = c("center", "scale") )
testData <- predict( testData_range_model, newdata = x_test) 
testData[1:5, 1:5]

testData$class <- y_test
# check head of the test data
head(testData)

dim(x_valid)

validData_range_model <- preProcess( x_valid, method = c("center", "scale"))
validData <- predict( validData_range_model, newdata = x_valid)
validData[1:5, 1:5]
validData$class <- y_valid


dim(x_apply)
applicationData_range_model <- preProcess( x_apply, method = c("center", "scale"))
applicationData <- predict( applicationData_range_model, newdata = x_apply)
applicationData[1:5, 1:5]
applicationData$class <- y_apply

## 
## check dim for both train and test data
dim( trainData )
dim( testData )
dim( validData )
dim( applicationData )

colnames(trainData)

########################################################################
## Section IV briefly check the feature weight and density
## call featurePlot() function 
######################################################################## 
##
## 
## box plot, to see significantly differential expression genes
trainData$Row.names <- NULL
testData$Row.names <- NULL

featurePlot(x = trainData[, 1:15], 
            y = trainData$class, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free"))
)


## density plot, to visulize more important variables 
featurePlot(x = trainData[, 1:15], 
            y = trainData$class, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))



### 
## Plot correlation of clinical variables
## we have to make sure there's no bias variables introduced into the model. 

library(corrplot)


##########################################################################################################  
##### 
##### Section V: RFE 
##### The most important part of this code
##### 
##########################################################################################################  


##########################################################################################################  
##  RFE works in 3 broad steps:
##    
##  Step 1: Build a ML model on a training dataset and estimate the feature importances on the test dataset.
##          Here in this case, the features are genes; 
## 
##  Step 2: Keeping priority to the most important variables, 
##          iterate through by building models of given subset sizes, 
##          that is, subgroups of most important predictors determined from step 1. 
##          Ranking of the predictors is recalculated in each iteration.
##  
##  Step 3: The model performances are compared across different subset sizes 
##          to arrive at the optimal number and list of final predictors.
##  
##  Stop 4: It can be implemented using the rfe() function and you have the flexibility 
##          to control what algorithm rfe uses and how it cross validates by defining the rfeControl().
##########################################################################################################  




##########################################################################################################  
#### Section V, Step 3  
#### Define the training control

fitControl <- trainControl(
  method = 'repeatedcv',           # k-fold cross validation
  number = 25,                     # number of folds
  repeats = 10,                    # number of repeats
  savePredictions = T,             # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  #,  # results summary function
  # savePredictions = 'final'
) 


#################################################################################################
#### Section V, step 4
#### Tune hyper parameters by setting tuneLength



#############################
#### Step 4.1 choose svm_linear method
set.seed(1234)
model_svmLinear = train(class ~ ., 
                        data=trainData, 
                        method='svmLinear', 
                        tuneLength = 20, 
                        metric='ROC', 
                        trControl = fitControl
      )





#############################
## briefly check the svmLinear results
model_svmLinear 

# model_svmLinear$pred$yes
# model_svmLinear$pred$no


varimp_svmLinear <- varImp(model_svmLinear)
plot(varimp_svmLinear, main="Variable Importance of Merged Data with svmLinear")

varImp(model_svmLinear)
str( varimp_svmLinear )
varimp_svmLinear[1]


  ## create a new dataframe with variable importance evaluated by the svm_Linear model
  varimp.df <- varimp_svmLinear[1]$importance
  
  str(varimp.df)
  
  ## sort the variable importance dataframe by Tumor column in ascending orders
  ## in this way, we could pick top 50% genes based on the ML model
  ## then loop back to line 220 to replace the top genes. 
  varimp.sor <- varimp.df[ with( varimp.df, order(-Tumor)), ]
  varimp.sor
  top.genes.50 <- row.names(varimp.sor)[1:50]
  top.genes.50
  
  ## save top 50 ranked miRNA genes
  sink("C:/Users/Jeff/SH_Hospital/SVM_RFE_Prj/Jeff_Shuai/20200115_input_files/Counts/top.genes.50.txt")
  print(top.genes.50)
  sink()
  
####### 
  ## after getting the top 50 ranked miRNAs, go back to code line 218, run the whole workflow again; 
    
  top.genes.21 <- row.names(varimp.sor)[1:21]
  
  sink("C:/Users/Jeff/SH_Hospital/SVM_RFE_Prj/Jeff_Shuai/20200115_input_files/Counts/top.genes.21.txt")
  print(top.genes.21)
  sink()
  
########################################################################################
## step 4.2
## plot roc for the training data
## We can calculate the area under the curve...
## Select a parameter setting
## selectedIndices <- model_mars2$pred

rocobj_svmlinear <- roc(model_svmLinear$pred$obs, model_svmLinear$pred$Tumor, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="Combined Model svmLinear",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=5, 
                        print.auc=TRUE,
                        print.auc.y = 25)





############################

###############################################################
## Step 4.3  Predict on testData and Compute the confusion matrix
## 

predicted2 <- predict(model_svmLinear, testData)
predicted.valid <- predict(model_svmLinear, validData) 
predicted.apply <- predict( model_svmLinear, applicationData)

################################################################
## print out confusion matrix 

confusionMatrix(reference = testData$class, data = predicted2, mode='everything', positive='Tumor')

confusionMatrix(reference = validData$class, data = predicted.valid, mode='everything', positive='Tumor')

confusionMatrix(reference = applicationData$class, data = predicted.apply, mode='everything', positive='Tumor')

head( validData )

###################
## plot ROC-AUC for both training model and prediction result 

getwd()

tiff("Jan162020_Plot_svmLinear_Combined_All_4.tif", width = 8, height = 8, units = 'in', res = 300)

rocobj_svmlinear <- roc(model_svmLinear$pred$obs, model_svmLinear$pred$Tumor, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svm Linear ROC-AUC Plot",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=5, 
                        print.auc=TRUE,
                        print.auc.y = 25)

predicted_test <- predict(model_svmLinear, testData, type = "prob")

str( predicted_test )

rocobj_svmlinear.test <- roc( testData$class, predicted_test[[2]] , ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE,
                               #main="svmLinear",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="red", lwd=5,
                               
                               print.auc=TRUE,
                               print.auc.y = 30,
                               add=TRUE)

predicted_valid <- predict( model_svmLinear, validData, type = "prob")
rocobj_svmlinear.vliad <- roc( validData$class, predicted_valid[[2]] , ci=TRUE,
                              plot=TRUE, 
                              legacy.axes=TRUE, percent=TRUE,
                              main="svmLinear",
                              xlab="False Positive Percentage", 
                              ylab="True Postive Percentage", 
                              col="darkgreen", lwd=5,
                              
                              print.auc=TRUE,
                              print.auc.y = 35,
                              add=TRUE)


predicted_apply <- predict( model_svmLinear, applicationData, type = "prob")

rocobj_svmlinear.apply <- roc( applicationData$class, predicted_apply[[2]] , ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE,
                               main="svmLinear",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="brown", lwd=5,
                               
                               print.auc=TRUE,
                               print.auc.y = 40,
                               add=TRUE)


legend("bottomright", 
       legend=c( "svmLinear Training on TCGA", "svmLinear Testing on TCGA", "svmLinear Validation on GSE Data", 
                 "svmLinear Application on Clinical Data"), 
       col=c( "darkblue", "red", "darkgreen", "brown"), 
       lwd=4
)


dev.off()

predicted_prob



##################################################################################
####  Section VII
####  ensemble predictions from multiple models using caretEnsemble
#### 

library(caretEnsemble)

# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainControl(method="repeatedcv", 
                             number=15, 
                             repeats=5,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

fitControl <- trainControl(
  method = 'repeatedcv',           # k-fold cross validation
  number = 15,                     # number of folds
  repeats = 5,                    # number of repeats
  savePredictions = T,             # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  #,  # results summary function
  # savePredictions = 'final'
) 

algorithmList <- c('rf', 'knn', 'earth','xgbDART', 'svmRadial', 'svmLinear')

## MARS does not work well, so we replaced it with svmPoly

algorithmList <- c('rf', 'earth','xgbDART', 'svmPoly','svmRadial', 'svmLinear')

# a quick test for algorithm list 
# algorithmList <- c('rf', 'knn', 'svmRadial', 'svmLinear', 'earth')

################################################################################
## Run all algorithms in the list: 
set.seed(123)

merged.models <- caretList(class ~ ., 
                           data=trainData, 
                           metric='ROC',
                           trControl=fitControl, 
                           methodList=algorithmList
) 






################################################################# 
#  alternatively, we can use traniControl; 
# trainControl 
# set.see(123)
# 
# merged.models <- caretList(class ~ ., 
#                           data=trainData, 
#                           metric='ROC',
#                           trControl=trainControl, 
#                           methodList=algorithmList
# ) 



## sumarize resample results from 
summary( resamples(merged.models) )

################################################################################
## check resample() results
## 
results <- resamples(merged.models)
summary(results)


varImp(results)

#> summary( resamples(merged.models) )
#
#Call:
#  summary.resamples(object = resamples(merged.models))
#
#Models: rf, knn, svmRadial, svmLinear 
#Number of resamples: 75 
#
#ROC 
#Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#rf           0       1      1 0.8769231       1    1   10
#knn          0       1      1 0.8500000       1    1   10
#svmRadial    0       1      1 0.9000000       1    1   10
#svmLinear    0       1      1 0.9692308       1    1   10
#
#Sens 
#Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#rf           0     1.0      1 0.9066667       1    1    0
#knn          0     1.0      1 0.9666667       1    1    0
#svmRadial    0     0.5      1 0.7533333       1    1    0
#svmLinear    0     1.0      1 0.9000000       1    1    0
#
#Spec 
#Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#rf           0       0      1 0.5384615       1    1   10
#knn          0       0      0 0.2615385       1    1   10
#svmRadial    0       0      1 0.7384615       1    1   10
#svmLinear    0       0      1 0.7230769       1    1   10


################################################################################
## Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

## Save as multi_algo_Accuracy_Kappa_boxplot



################################################################################
## Plot multi ROCs in one plot
rocobj_models <- roc(merged.models$rf$pred$obs, 
                     merged.models$rf$pred$Tumor, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main = "Combined Data Multi Model ROCs",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(merged.models$svmRadial$pred$obs, 
                     merged.models$svmRadial$pred$Tumor, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="green", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 44,
                     add = TRUE
)


rocobj_models <- roc(merged.models$svmLinear$pred$obs, 
                     merged.models$svmLinear$pred$Tumor, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 48,
                     add = TRUE
)

rocobj_models <- roc(merged.models$xgbDART$pred$obs, 
                     merged.models$xgbDART$pred$Tumor, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="black", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 52,
                     add = TRUE
)

rocobj_models <- roc(merged.models$svmPoly$pred$obs, 
                     merged.models$svmPoly$pred$Tumor, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="yellow", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 56,
                     add = TRUE
)


rocobj_models <- roc(merged.models$knn$pred$obs, 
                     merged.models$knn$pred$Tumor, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="pink", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 60,
                     add = TRUE
)

# rocobj_models <- roc(merged.models$ada$pred$obs, 
#                     merged.models$ada$pred$yes, 
#                     ci=TRUE,
#                     plot=TRUE, 
#                     legacy.axes=TRUE, percent=TRUE, 
#                     xlab="False Positive Percentage", 
#                     ylab="True Postive Percentage", 
#                     col="purple", lwd=4, 
#                     print.auc=TRUE,
#                     print.auc.y = 64,
#                     add = TRUE
#                     )


legend("bottomright", 
       legend=c( "RandomForest", "svmRadial", "svmLinear", "xgbDART", "svmPoly", "KNN" ), 
       col=c( "darkblue", "green", "red", "black", "yellow", "pink"), 
       lwd=4
)


###################################################
## save 




merged.models$rf

merged.models$knn

merged.models$svmPoly

merged.models$svmLinear

merged.models$svmRadial 

merged.models$xgbDART

t( coords( smooth( rocobj_models_svmL ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))





### svm linear 
################################################################################
## Plot FOUR ROCs (training, test, valid, application) in one plot
## Plot multi ROCs in one plot

rocobj_models_svmL <- roc(merged.models$svmLinear$pred$obs, 
                     merged.models$svmLinear$pred$Tumor, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     main = "ROC-AUC of svmLinear",
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 25
)

# rocobj_models.rf.valid <- roc(merged.models.valid$rf$pred$obs, 
#                              merged.models.valid$rf$pred$yes, 
#                              ci=TRUE,
#                              plot=TRUE, 
#                              legacy.axes=TRUE, percent=TRUE, 
#                              xlab="False Positive Percentage", 
#                              ylab="True Postive Percentage", 
#                              col="red", lwd=4, 
#                              print.auc=TRUE,
#                              print.auc.y = 30,
#                              add = TRUE
#                             )


merged.models$svmLinear

## on test data
predicted_test <- predict(merged.models$svmLinear, testData, type = "prob")

str( predicted_test )

rocobj_svmlinear.test <- roc( testData$class, predicted_test[[2]] , ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE,
                               main="svmLinear",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="red", lwd=5,
                               
                               print.auc=TRUE,
                               print.auc.y = 30,
                               add=TRUE)

## on GSE33858 validation data

predicted_valid <- predict(merged.models$svmLinear, validData, type = "prob")

str( predicted_valid )

rocobj_svmlinear.valid <- roc( validData$class, predicted_valid[[2]] , ci=TRUE,
                              plot=TRUE, 
                              legacy.axes=TRUE, percent=TRUE,
                              main="svmLinear",
                              xlab="False Positive Percentage", 
                              ylab="True Postive Percentage", 
                              col="green", lwd=5,
                              
                              print.auc=TRUE,
                              print.auc.y = 35,
                              add=TRUE)


predicted_apply <- predict(merged.models$svmLinear, applicationData, type = "prob")

str( predicted_apply )

rocobj_svmlinear.apply <- roc( applicationData$class, predicted_apply[[2]] , ci=TRUE,
                              plot=TRUE, 
                              legacy.axes=TRUE, percent=TRUE,
                              main="svmLinear",
                              xlab="False Positive Percentage", 
                              ylab="True Postive Percentage", 
                              col="brown", lwd=5,
                              
                              print.auc=TRUE,
                              print.auc.y = 40,
                              add=TRUE)



legend("bottomright", 
       legend=c( "SVM-Linear Training TCGA","SVM-Linear Testing TCGA", "SVM-Linear Validation GSE", "SVM-Linear Application Clinic" ), 
       col=c( "darkblue", "red", "green", "brown"), 
       lwd=4
)


predicted_test <- predict(merged.models$svmLinear, testData)


### 
predicted_valid <- predict(merged.models$svmLinear, validData)
predicted_apply <- predict(merged.models$svmLinear, applicationData)

################################################################
## print out confusion matrix 

confusionMatrix(reference = testData$class,        data = predicted_test, mode='everything', positive='Tumor')
confusionMatrix(reference = validData$class,       data = predicted_valid, mode='everything', positive='Tumor')
confusionMatrix(reference = applicationData$class, data = predicted_apply, mode='everything', positive='Tumor')

############################
## END svmLinear
###########################



### svm Radial 
################################################################################
## Plot FOUR ROCs (training, test, valid, application) in one plot
## Plot multi ROCs in one plot

rocobj_models_svmR <- roc(merged.models$svmRadial$pred$obs, 
                          merged.models$svmRadial$pred$Tumor, 
                          ci=TRUE,
                          plot=TRUE, 
                          legacy.axes=TRUE, percent=TRUE, 
                          main = "ROC-AUC of svmRadial",
                          xlab="False Positive Percentage", 
                          ylab="True Postive Percentage", 
                          col="darkblue", lwd=4, 
                          print.auc=TRUE,
                          print.auc.y = 25
)



merged.models$svmRadial

## on test data
predicted_test <- predict(merged.models$svmRadial, testData, type = "prob")

str( predicted_test )

rocobj_svmRadial.test <- roc( testData$class, predicted_test[[2]] , ci=TRUE,
                              plot=TRUE, 
                              legacy.axes=TRUE, percent=TRUE,
                              main="svmRadial",
                              xlab="False Positive Percentage", 
                              ylab="True Postive Percentage", 
                              col="red", lwd=5,
                              
                              print.auc=TRUE,
                              print.auc.y = 30,
                              add=TRUE)

## on GSE33858 validation data

predicted_valid <- predict(merged.models$svmRadial, validData, type = "prob")

str( predicted_valid )

rocobj_svmRadial.valid <- roc( validData$class, predicted_valid[[2]] , ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE,
                               main="svmLinear",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="green", lwd=5,
                               
                               print.auc=TRUE,
                               print.auc.y = 35,
                               add=TRUE)


predicted_apply <- predict(merged.models$svmRadial, applicationData, type = "prob")

str( predicted_apply )

rocobj_svmRadial.apply <- roc( applicationData$class, predicted_apply[[2]] , ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE,
                               main="svmLinear",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="brown", lwd=5,
                               
                               print.auc=TRUE,
                               print.auc.y = 40,
                               add=TRUE)



legend("bottomright", 
       legend=c( "SVM-Radial Training TCGA","SVM-Radial Testing TCGA", "SVM-Radial Validation GSE", "SVM-Radial Application Clinic" ), 
       col=c( "darkblue", "red", "green", "brown"), 
       lwd=4
)


predicted_test <- predict(merged.models$svmRadial, testData)


### 
predicted_valid <- predict(merged.models$svmRadial, validData)
predicted_apply <- predict(merged.models$svmRadial, applicationData)

################################################################
## print out confusion matrix 

confusionMatrix(reference = testData$class,        data = predicted_test,  mode='everything', positive='Tumor')
confusionMatrix(reference = validData$class,       data = predicted_valid, mode='everything', positive='Tumor')
confusionMatrix(reference = applicationData$class, data = predicted_apply, mode='everything', positive='Tumor')

#############
## END svm Radial
#############





#### random forest 
################################################################################
## Plot FOUR ROCs (training, test, valid, application) in one plot
## Plot multi ROCs in one plot

rocobj_models_rf <- roc(merged.models$rf$pred$obs, 
                          merged.models$rf$pred$Tumor, 
                          ci=TRUE,
                          plot=TRUE, 
                          legacy.axes=TRUE, percent=TRUE, 
                          main = "ROC-AUC of Random Forest",
                          xlab="False Positive Percentage", 
                          ylab="True Postive Percentage", 
                          col="darkblue", lwd=4, 
                          print.auc=TRUE,
                          print.auc.y = 25
)



merged.models$rf

## on test data
predicted_test <- predict(merged.models$rf, testData, type = "prob")

str( predicted_test )

rocobj_rf.test <- roc( testData$class, predicted_test[[2]] , ci=TRUE,
                              plot=TRUE, 
                              legacy.axes=TRUE, percent=TRUE,
                              main="svmRadial",
                              xlab="False Positive Percentage", 
                              ylab="True Postive Percentage", 
                              col="red", lwd=5,
                              
                              print.auc=TRUE,
                              print.auc.y = 30,
                              add=TRUE)

## on GSE33858 validation data

predicted_valid <- predict(merged.models$rf, validData, type = "prob")

str( predicted_valid )

rocobj_rf.valid <- roc( validData$class, predicted_valid[[2]] , ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE,
                               main="svmLinear",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="green", lwd=5,
                               
                               print.auc=TRUE,
                               print.auc.y = 35,
                               add=TRUE)


predicted_apply <- predict(merged.models$rf, applicationData, type = "prob")

str( predicted_apply )

rocobj_rf.apply <- roc( applicationData$class, predicted_apply[[2]] , ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE,
                               main="svmLinear",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="brown", lwd=5,
                               
                               print.auc=TRUE,
                               print.auc.y = 40,
                               add=TRUE)



legend("bottomright", 
       legend=c( "RF Training TCGA","RF Testing TCGA", "RF Validation GSE", "RF Application Clinic" ), 
       col=c( "darkblue", "red", "green", "brown"), 
       lwd=4
)


predicted_test <- predict(merged.models$rf, testData)


### 
predicted_valid <- predict(merged.models$rf, validData)
predicted_apply <- predict(merged.models$rf, applicationData)

################################################################
## print out confusion matrix 

confusionMatrix(reference = testData$class,        data = predicted_test,  mode='everything', positive='Tumor')
confusionMatrix(reference = validData$class,       data = predicted_valid, mode='everything', positive='Tumor')
confusionMatrix(reference = applicationData$class, data = predicted_apply, mode='everything', positive='Tumor')

#############
## END Random Forest
#############






### KNN is the worst
################################################################################
## Plot FOUR ROCs (training, test, valid, application) in one plot
## Plot multi ROCs in one plot

rocobj_models_knn <- roc(merged.models$knn$pred$obs, 
                        merged.models$knn$pred$Tumor, 
                        ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE, 
                        main = "ROC-AUC of KNN",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE,
                        print.auc.y = 25
)



merged.models$knn

## on test data
predicted_test <- predict(merged.models$knn, testData, type = "prob")

str( predicted_test )

rocobj_knn.test <- roc( testData$class, predicted_test[[2]] , ci=TRUE,
                       plot=TRUE, 
                       legacy.axes=TRUE, percent=TRUE,
                       main="svmRadial",
                       xlab="False Positive Percentage", 
                       ylab="True Postive Percentage", 
                       col="red", lwd=5,
                       
                       print.auc=TRUE,
                       print.auc.y = 30,
                       add=TRUE)

## on GSE33858 validation data

predicted_valid <- predict(merged.models$knn, validData, type = "prob")

str( predicted_valid )

rocobj_knn.valid <- roc( validData$class, predicted_valid[[2]] , ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="green", lwd=5,
                        
                        print.auc=TRUE,
                        print.auc.y = 35,
                        add=TRUE)


predicted_apply <- predict(merged.models$knn, applicationData, type = "prob")

str( predicted_apply )

rocobj_rf.apply <- roc( applicationData$class, predicted_apply[[2]] , ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="brown", lwd=5,
                        
                        print.auc=TRUE,
                        print.auc.y = 40,
                        add=TRUE)



legend("bottomright", 
       legend=c( "KNN Training TCGA","KNN Testing TCGA", "KNN Validation GSE", "KNN Application Clinic" ), 
       col=c( "darkblue", "red", "green", "brown"), 
       lwd=4
)


predicted_test <- predict(merged.models$knn, testData)


### 
predicted_valid <- predict(merged.models$knn, validData)
predicted_apply <- predict(merged.models$knn, applicationData)

################################################################
## print out confusion matrix 

confusionMatrix(reference = testData$class,        data = predicted_test,  mode='everything', positive='Tumor')
confusionMatrix(reference = validData$class,       data = predicted_valid, mode='everything', positive='Tumor')
confusionMatrix(reference = applicationData$class, data = predicted_apply, mode='everything', positive='Tumor')

#############
## END KNN
#############



### svmPoly
################################################################################
## Plot FOUR ROCs (training, test, valid, application) in one plot
## Plot multi ROCs in one plot

rocobj_models_svmP <- roc(merged.models$svmPoly$pred$obs, 
                         merged.models$svmPoly$pred$Tumor, 
                         ci=TRUE,
                         plot=TRUE, 
                         legacy.axes=TRUE, percent=TRUE, 
                         main = "ROC-AUC of svmPoly",
                         xlab="False Positive Percentage", 
                         ylab="True Postive Percentage", 
                         col="darkblue", lwd=4, 
                         print.auc=TRUE,
                         print.auc.y = 25
)



merged.models$svmPoly

## on test data
predicted_test <- predict(merged.models$svmPoly, testData, type = "prob")

str( predicted_test )

rocobj_mars.test <- roc( testData$class, predicted_test[[2]] , ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmRadial",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="red", lwd=5,
                        
                        print.auc=TRUE,
                        print.auc.y = 30,
                        add=TRUE)

## on GSE33858 validation data

predicted_valid <- predict(merged.models$svmPoly, validData, type = "prob")

str( predicted_valid )

rocobj_mars.valid <- roc( validData$class, predicted_valid[[2]] , ci=TRUE,
                         plot=TRUE, 
                         legacy.axes=TRUE, percent=TRUE,
                         main="svmLinear",
                         xlab="False Positive Percentage", 
                         ylab="True Postive Percentage", 
                         col="green", lwd=5,
                         
                         print.auc=TRUE,
                         print.auc.y = 35,
                         add=TRUE)


predicted_apply <- predict(merged.models$svmPoly, applicationData, type = "prob")

str( predicted_apply )

rocobj_marth.apply <- roc( applicationData$class, predicted_apply[[2]] , ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="brown", lwd=5,
                        
                        print.auc=TRUE,
                        print.auc.y = 40,
                        add=TRUE)



legend("bottomright", 
       legend=c( "svm-Poly Training TCGA","svm-Poly Testing TCGA", "svm-Poly Validation GSE", "svm-Poly Application Clinic" ), 
       col=c( "darkblue", "red", "green", "brown"), 
       lwd=4
)


predicted_test <- predict(merged.models$svmPoly, testData)


### 
predicted_valid <- predict(merged.models$svmPoly, validData)
predicted_apply <- predict(merged.models$svmPoly, applicationData)

################################################################
## print out confusion matrix 


merged.models$svmPoly

confusionMatrix(reference = testData$class,        data = predicted_test,  mode='everything', positive='Tumor')
confusionMatrix(reference = validData$class,       data = predicted_valid, mode='everything', positive='Tumor')
confusionMatrix(reference = applicationData$class, data = predicted_apply, mode='everything', positive='Tumor')

#############
## END MARS
#############



### xgbDART
################################################################################
## Plot FOUR ROCs (training, test, valid, application) in one plot
## Plot multi ROCs in one plot

rocobj_models_xgb <- roc(merged.models$xgbDART$pred$obs, 
                          merged.models$xgbDART$pred$Tumor, 
                          ci=TRUE,
                          plot=TRUE, 
                          legacy.axes=TRUE, percent=TRUE, 
                          main = "ROC-AUC of xgbDART",
                          xlab="False Positive Percentage", 
                          ylab="True Postive Percentage", 
                          col="darkblue", lwd=4, 
                          print.auc=TRUE,
                          print.auc.y = 25
)



merged.models$xgbDART

## on test data
predicted_test <- predict(merged.models$xgbDART, testData, type = "prob")

str( predicted_test )

rocobj_xgb.test <- roc( testData$class, predicted_test[[2]] , ci=TRUE,
                         plot=TRUE, 
                         legacy.axes=TRUE, percent=TRUE,
                         main="svmRadial",
                         xlab="False Positive Percentage", 
                         ylab="True Postive Percentage", 
                         col="red", lwd=5,
                         
                         print.auc=TRUE,
                         print.auc.y = 30,
                         add=TRUE)

## on GSE33858 validation data

predicted_valid <- predict(merged.models$xgbDART, validData, type = "prob")

str( predicted_valid )

rocobj_xgb.valid <- roc( validData$class, predicted_valid[[2]] , ci=TRUE,
                          plot=TRUE, 
                          legacy.axes=TRUE, percent=TRUE,
                          main="svmLinear",
                          xlab="False Positive Percentage", 
                          ylab="True Postive Percentage", 
                          col="green", lwd=5,
                          
                          print.auc=TRUE,
                          print.auc.y = 35,
                          add=TRUE)


predicted_apply <- predict(merged.models$xgbDART, applicationData, type = "prob")

str( predicted_apply )

rocobj_xgb.apply <- roc( applicationData$class, predicted_apply[[2]] , ci=TRUE,
                           plot=TRUE, 
                           legacy.axes=TRUE, percent=TRUE,
                           main="svmLinear",
                           xlab="False Positive Percentage", 
                           ylab="True Postive Percentage", 
                           col="brown", lwd=5,
                           
                           print.auc=TRUE,
                           print.auc.y = 40,
                           add=TRUE)



legend("bottomright", 
       legend=c( "xgbDART Training TCGA","xgbDART Testing TCGA", "xgbDART Validation GSE", "xgbDART Application Clinic" ), 
       col=c( "darkblue", "red", "green", "brown"), 
       lwd=4
)


predicted_test <- predict(merged.models$xgbDART, testData)


### 
predicted_valid <- predict(merged.models$xgbDART, validData)
predicted_apply <- predict(merged.models$xgbDART, applicationData)

################################################################
## print out confusion matrix 


merged.models$xgbDART

confusionMatrix(reference = testData$class,        data = predicted_test,  mode='everything', positive='Tumor')
confusionMatrix(reference = validData$class,       data = predicted_valid, mode='everything', positive='Tumor')
confusionMatrix(reference = applicationData$class, data = predicted_apply, mode='everything', positive='Tumor')

#############
## END xgbDART
#############



























################################################################################
## Plot multi ROCs in one plot
  rocobj_models.svmRadial <- roc(merged.models$svmRadial$pred$obs, 
                                 merged.models$svmRadial$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 main="svmRadial Combined ROC",
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="darkblue", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 40
  )
  
  rocobj_models.svmRadial.valid <- roc(merged.models.valid$svmRadial$pred$obs, 
                                       merged.models.valid$svmRadial$pred$yes, 
                                       ci=TRUE,
                                       plot=TRUE, 
                                       legacy.axes=TRUE, percent=TRUE, 
                                       xlab="False Positive Percentage", 
                                       ylab="True Postive Percentage", 
                                       col="red", lwd=4, 
                                       print.auc=TRUE,
                                       print.auc.y = 30,
                                       add = TRUE
  )
  
  
  
  
  legend("bottomright", 
         legend=c( "svmRadial Training Set", "svmRadial Validation" ), 
         col=c( "darkblue", "red" ), 
         lwd=4
  )
  
  
  #######################
  
  
  rocobj_models.svmLinear <- roc(merged.models$svmLinear$pred$obs, 
                                 merged.models$svmLinear$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 main="svmLinear Combined ROC",
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 48
  )
  
  
  ## Plot multi ROCs in one plot
  rocobj_models.svmLinear <- roc(merged.models$svmLinear$pred$obs, 
                                 merged.models$svmLinear$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 main="svmLinear Combined ROC",
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="darkblue", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 40
  )
  
  rocobj_models.svmLinear.valid <- roc(merged.models.valid$svmLinear$pred$obs, 
                                       merged.models.valid$svmLinear$pred$yes, 
                                       ci=TRUE,
                                       plot=TRUE, 
                                       legacy.axes=TRUE, percent=TRUE, 
                                       xlab="False Positive Percentage", 
                                       ylab="True Postive Percentage", 
                                       col="red", lwd=4, 
                                       print.auc=TRUE,
                                       print.auc.y = 30,
                                       add = TRUE
  )
  
  
  
  
  legend("bottomright", 
         legend=c( "svmLinear Training Set", "svmLinear Validation" ), 
         col=c( "darkblue", "red" ), 
         lwd=4
  )
  
  #######################################################
  
  
  rocobj_models.xgbDART <- roc(merged.models$xgbDART$pred$obs, 
                               merged.models$xgbDART$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               main="xgbDART Clinical ROC",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="black", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 52,
                               add = TRUE
  )
  
  
  ## Plot multi ROCs in one plot
  rocobj_models.xgbDART <- roc(merged.models$xgbDART$pred$obs, 
                               merged.models$xgbDART$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               main="xgbDART Combined ROC",
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 40
  )
  
  rocobj_models.xgbDART.valid <- roc(merged.models.valid$xgbDART$pred$obs, 
                                     merged.models.valid$xgbDART$pred$yes, 
                                     ci=TRUE,
                                     plot=TRUE, 
                                     legacy.axes=TRUE, percent=TRUE, 
                                     xlab="False Positive Percentage", 
                                     ylab="True Postive Percentage", 
                                     col="red", lwd=4, 
                                     print.auc=TRUE,
                                     print.auc.y = 30,
                                     add = TRUE
  )
  
  
  
  
  legend("bottomright", 
         legend=c( "xgbDART Training Set", "xgbDART Validation Set" ), 
         col=c( "darkblue", "red" ), 
         lwd=4
  )
  
  
  
  
  ####################################################
  rocobj_models <- roc(clin.models$earth$pred$obs, 
                       clin.models$earth$pred$yes, 
                       ci=TRUE,
                       plot=TRUE, 
                       legacy.axes=TRUE, percent=TRUE, 
                       xlab="False Positive Percentage", 
                       ylab="True Postive Percentage", 
                       col="yellow", lwd=4, 
                       print.auc=TRUE,
                       print.auc.y = 56,
                       add = TRUE
  )
  
  ## Plot multi ROCs in one plot
  rocobj_models.mars <- roc(merged.models$earth$pred$obs, 
                            merged.models$earth$pred$yes, 
                            ci=TRUE,
                            plot=TRUE, 
                            legacy.axes=TRUE, percent=TRUE, 
                            main="MARS Combined ROC",
                            xlab="False Positive Percentage", 
                            ylab="True Postive Percentage", 
                            col="darkblue", lwd=4, 
                            print.auc=TRUE,
                            print.auc.y = 40
  )
  
  rocobj_models.mars.valid <- roc(merged.models.valid$earth$pred$obs, 
                                  merged.models.valid$earth$pred$yes, 
                                  ci=TRUE,
                                  plot=TRUE, 
                                  legacy.axes=TRUE, percent=TRUE, 
                                  xlab="False Positive Percentage", 
                                  ylab="True Postive Percentage", 
                                  col="red", lwd=4, 
                                  print.auc=TRUE,
                                  print.auc.y = 30,
                                  add = TRUE
  )
  
  
  
  
  legend("bottomright", 
         legend=c( "MARS Training Set", "MARS Validation" ), 
         col=c( "darkblue", "red" ), 
         lwd=4
  )
  
  
  
  ##############################################################
  rocobj_models.knn <- roc(clin.models$knn$pred$obs, 
                           clin.models$knn$pred$yes, 
                           ci=TRUE,
                           plot=TRUE, 
                           legacy.axes=TRUE, percent=TRUE, 
                           xlab="False Positive Percentage", 
                           ylab="True Postive Percentage", 
                           col="pink", lwd=4, 
                           print.auc=TRUE,
                           print.auc.y = 35.1,
                           add = TRUE
  )
  
  ###
  ## Plot multi ROCs in one plot
  rocobj_models.knn <- roc(merged.models$knn$pred$obs, 
                           merged.models$knn$pred$yes, 
                           ci=TRUE,
                           plot=TRUE, 
                           legacy.axes=TRUE, percent=TRUE, 
                           main="KNN Combined ROC",
                           xlab="False Positive Percentage", 
                           ylab="True Postive Percentage", 
                           col="darkblue", lwd=4, 
                           print.auc=TRUE,
                           print.auc.y = 40
  )
  
  rocobj_models.knn.valid <- roc(merged.models.valid$knn$pred$obs, 
                                 merged.models.valid$knn$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 30,
                                 add = TRUE
  )
  
  
  
  
  legend("bottomright", 
         legend=c( "KNN Training Set", "KNN Validation Set" ), 
         col=c( "darkblue", "red" ), 
         lwd=4
  )
  
  
  ## END UPDATE 0525 2019 
  ## THE amazing client wanted ROC-AUCs for both training set and validation set, 
  ## instead of confusion matrix;
  ############################################################################################
  
  
  ## PRINT OUT model spec/sens/roc/ppn/cpn table:
  
  t( coords( smooth( rocobj_models.knn ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  t( coords( smooth( rocobj_models.knn.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  
  t( coords( smooth( rocobj_models.mars ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  t( coords( smooth( rocobj_models.mars.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  
  t( coords( smooth( rocobj_models.svmLinear ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  t( coords( smooth( rocobj_models.svmLinear.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  
  t( coords( smooth( rocobj_models.svmRadial ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  t( coords( smooth( rocobj_models.svmRadial.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  
  t( coords( smooth( rocobj_models.rf ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  t( coords( smooth( rocobj_models.rf.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  
  t( coords( smooth( rocobj_models.xgbDART ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  t( coords( smooth( rocobj_models.xgbDART.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  
  t( coords( smooth( rocobj_models.knn ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  t( coords( smooth( rocobj_models.knn.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  
  
  # 
  # > t( coords( smooth( rocobj_models.knn ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity     ppv   npv
  # [1,]    99.99285         0.9 98.7924 60.82
  # > t( coords( smooth( rocobj_models.knn.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]    99.99794         0.9 99.65745 60.21629
  # > t( coords( smooth( rocobj_models.mars ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity ppv     npv
  # [1,]         100         0.9 100 60.8217
  # > t( coords( smooth( rocobj_models.mars.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity     ppv      npv
  # [1,]     99.2355         0.9 43.9723 60.03279
  # > t( coords( smooth( rocobj_models.svmLinear ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]    99.99764         0.9 99.59895 60.82114
  # > t( coords( smooth( rocobj_models.svmLinear.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]    99.82948         0.9 77.86917 60.17589
  # > t( coords( smooth( rocobj_models.svmRadial ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]    99.98568         0.9 97.61021 60.81829
  # > t( coords( smooth( rocobj_models.svmRadial.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]    99.89311         0.9 84.87855 60.19116
  # > t( coords( smooth( rocobj_models.rf ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity     ppv     npv
  # [1,]         100         0.9 99.9998 60.8217
  # > t( coords( smooth( rocobj_models.rf.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]    99.98276         0.9 97.20752 60.21265
  # > t( coords( smooth( rocobj_models.xgbDART ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]     99.9999         0.9 99.98294 60.82168
  # > t( coords( smooth( rocobj_models.xgbDART.valid ), x = 0.9, input = "sensitivity", ret = c("specificity", "sensitivity", "ppv", "npv") ))
  # specificity sensitivity      ppv      npv
  # [1,]    99.99739         0.9 99.56637 60.21615
  # > 
  
  getwd()
  save.image(file = "C:/Users/Jeff/SH_Hospital/SVM_RFE_Prj/Jeff_Shuai/20200115_input_files/workplace_7models.RData")


  
  
  
  
  
  
  
################################################################################
####  Section VIII 
####  combine the predictions of multiple models to form a final prediction


# Create the trainControl
set.seed(100)

stackControl <- trainControl(method="repeatedcv", 
                             number=20, 
                             repeats=5,
                             savePredictions=TRUE, 
                             classProbs=TRUE
)

# Ensemble the predictions of `models` to form a new combined prediction based on glm
stack.glm <- caretStack(merged.models, method="glm", metric="Accuracy", trControl=stackControl)

stack.glm
summary(stack.glm)
varImp(stack.glm)


### the current release of cretStack() would reverse the prediction
### SO,,, we have to re-assign Tumor vs. Normal to the group. 

testData.stack <- testData
validData.stack <- validData
applyData.stack <- applicationData

#luad_data$class <- ifelse(luad_data$class=="Solid Tissue Normal", "Normal", "Tumor")
testData.stack$class  <- ifelse(testData.stack$class  == "Tumor", "Normal", "Tumor")
validData.stack$class <- ifelse(validData.stack$class == "Tumor", "Normal", "Tumor")
applyData.stack$class <- ifelse(applyData.stack$class == "Tumor", "Normal", "Tumor")

testData.stack$class <- as.factor(testData.stack$class)
summary(testData.stack$class)
validData.stack$class <- as.factor(validData.stack$class)
applyData.stack$class <- as.factor(applyData.stack$class)

# Predict on testData 
stack_predicteds <- predict(stack.glm, newdata=testData.stack) 
head(stack_predicteds) 

confusionMatrix(reference = testData.stack$class, data = stack_predicteds, mode='everything', positive='Tumor') 


# predict on validation data GSE33858
stack_predicted.valid <- predict(stack.glm, newdata=validData.stack) 
head(stack_predicted.valid) 

confusionMatrix(reference = validData.stack$class, data = stack_predicted.valid, mode='everything', positive='Tumor') 


# predict on application clinical data
stack_predicted.apply <- predict(stack.glm, newdata=applyData.stack) 
head(stack_predicted.apply) 

confusionMatrix(reference = applyData.stack$class, data = stack_predicted.apply, mode='everything', positive='Tumor') 


varImp(merged.models)
summary(stack.glm)

plot(stack.glm)

################################################################################
## Plot TWO ROCs in one plot
## Plot multi ROCs in one plot
rocobj_models_stack <- roc(stack.glm$pred$obs, 
                          stack.glm$pred$Tumor, 
                          ci=TRUE,
                          plot=TRUE, 
                          legacy.axes=TRUE, percent=TRUE, 
                          main = "Combined Data Stacked Model ROCs",
                          xlab="False Positive Percentage", 
                          ylab="True Postive Percentage", 
                          col="darkblue", lwd=4, 
                          print.auc=TRUE,
                          print.auc.y = 25
)

stack.glm$models

models <- caretList()
ens <- caretEnsemble(merged.models)
plot(ens)


xyplot(resamples(merged.models), what = "BlandAltman")
xyplot(results, what = "BlandAltman")
splom( results )

########################################################################3
####
####         END OF THE MAIN CODE SECTION                           ####3
####
########################################################################3
































































































my_data$class <- design$Progressed
my_data$class <- ifelse(my_data$class==1, "yes", "no")

my_data[1:5, 1:5]



sort( row.names(my_data) )
sort( row.names(clinical.data) )


## in case there are multi factors in the class columns, remove them;
## subset the data, only keep primiary tumor and solid normal
my_data <- subset( my_data, class!='Metastatic') # my_data$class != "Metastatic", ]

## check factor levels left
summary(my_data$class) 

## drop the empty Metastatic level
my_data$class <- factor(my_data$class)




####################################################################################
## Section II 
## Create the training and test datasets
####################################################################################

set.seed(1234)

###################################################################
## Section II, Step 1: Get row numbers for the training data
## The limitation in this case is: dataset too small
## So, we split the original data into 75% training + 25% testing
###################################################################

trainRowNumbers <- createDataPartition(my_data$class, p=0.60, list=FALSE)

set.seed(321)
testRowNumbers <- createDataPartition(my_data$class, p=0.55, list=FALSE)
testData <- my_data[testRowNumbers, ]

###################################################################
## Section II, Step 2: Create the training  dataset
trainData <- my_data[trainRowNumbers,]

###################################################################
## Section II, Step 3: Create the test dataset
#testData <- my_data[-trainRowNumbers,]

dim(trainData)
dim(testData)

# Store X and Y for later use.
x = trainData[, 1:dim(trainData)[2]-1]
y = trainData$class

x_test = testData[ , 1:dim(testData)[2] - 1]
y_test = testData$class
dim(x)
dim(x_test)

x$class
# x$class should be NULL


###############################################################################################
## Section III, data, mamipulation, normalization and scalling
###############################################################################################


########################################################################
## 
## convert all the numeric variables to range between 0 and 1, 
## by setting method=range in preProcess(). 
## preProcess_range_model <- preProcess(trainData, method='range')
## trainData <- predict(preProcess_range_model, newdata = trainData)
## 
## Append the Y variable
## trainData$class <- y
########################################################################


########################################################################
## Section III, step 1 
## remove constant columns
## for top 500/1000 significant genes, this step could be ignored; 

x <- x[,apply(x, 2, var, na.rm=TRUE) > 0.0001]
dim(x)


########################################################################
## Section III, step 2
## convert all the numeric variables to range between 0 and 1, by setting method=range in preProcess(). 

trainData_range_model <- preProcess( x, method = c("center", "scale") )

trainData <- predict( trainData_range_model, newdata = x) 
trainData$class <- y

head(trainData)

########################################################################
## Section III, step 3
## convert all test data 

dim(x_test)

testData_range_model <- preProcess( x_test, method = c("center", "scale") )
testData <- predict( testData_range_model, newdata = x_test) 
testData$class <- y_test
# check head of the test data
head(testData)

## 
## check dim for both train and test data
dim( trainData )
dim( testData )


########################################################################
## Section IV briefly check the feature weight and density
## call featurePlot() function 
######################################################################## 
##
## 
## box plot, to see significantly differential expression genes 
featurePlot(x = trainData[, 1:15], 
            y = trainData$class, 
            plot = "box",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free"))
)


## density plot, to visulize more important variables 
featurePlot(x = trainData[, 1:15], 
            y = trainData$class, 
            plot = "density",
            strip=strip.custom(par.strip.text=list(cex=.7)),
            scales = list(x = list(relation="free"), 
                          y = list(relation="free")))



##########################################################################################################  
##### 
##### Section V: RFE 
##### The most important part of this code
##### 
##########################################################################################################  


##########################################################################################################  
##  RFE works in 3 broad steps:
##    
##  Step 1: Build a ML model on a training dataset and estimate the feature importances on the test dataset.
##          Here in this case, the features are genes; 
## 
##  Step 2: Keeping priority to the most important variables, 
##          iterate through by building models of given subset sizes, 
##          that is, subgroups of most important predictors determined from step 1. 
##          Ranking of the predictors is recalculated in each iteration.
##  
##  Step 3: The model performances are compared across different subset sizes 
##          to arrive at the optimal number and list of final predictors.
##  
##  Stop 4: It can be implemented using the rfe() function and you have the flexibility 
##          to control what algorithm rfe uses and how it cross validates by defining the rfeControl().
##########################################################################################################  


set.seed(1234)

options(warn=-1)

subsets <- c(1:45)

dim(trainData)

##########################################################################################################  
## Section V, Step 1
## build a control model on training dataset; estimate the feautre importances 

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",      ## cross-validation
                   number = 10,                ## fold of cross validation
                   repeats = 3,                ## repeats of cv
                   verbose = FALSE
)


#########################################
## Recursive Feature Elimination rfe() 

lmProfile <- rfe(x=trainData[, 1:45], 
                 y=trainData$class,
                 sizes = subsets,        ## model sizes (the number of most important genes to choose)
                 rfeControl = ctrl       ## use rfeControl output as reference: algorithm and cv to use;
)

## check lm_profile results
lmProfile 

## check the optimized variables (gens in this case)
lmProfile$optVariables

#####################################################
## tested by 4/15/2019 evening;
##################################################### 



################################################################################ 
#### Section V, step 2
#### 
#### train() the model and interpret the results
#### 
## Set the seed for reproducibility
set.seed(100)

################################################################################
## Train the model using MARS and predict on the training data itself.
## train() will cross validate the model 'earth' = Multivariate Adaptive Regression Splines (MARS); 
model_mars = train(class ~ ., 
                   data=trainData, 
                   method='earth'
                   #savePredictions = T
)

## check the results generated by train() with MARS algo:
model_mars
# model_mars$pred
# lmProfile$optsize

########################################
## fit the predicted model_mars
fitted <- predict(model_mars)

#######################################
## plot fit accuracies 
plot(model_mars, main="Model Accuracies with MARS") 

## 
## compute variable importance (genes more important to the model)

varimp_mars <- varImp(model_mars)
plot(varimp_mars, main="Variable Importance with MARS")
varimp_mars$importance 

## > varimp_mars$importance
## Overall
## PLAU                      100.00000
## CCL14                      81.73621
## IFITM1                     69.61859
## IL1RN                      42.42643
## ABCB1                       0.00000
## BST1                        0.00000
## C4B                         0.00000
## CCL17                       0.00000

##########################################################################################

##########################################################################################################  
#### Section V, Step 3  
#### Define the training control

fitControl <- trainControl(
  method = 'repeatedcv',           # k-fold cross validation
  number = 15,                     # number of folds
  repeats = 10,                    # number of repeats
  savePredictions = T,             # saves predictions for optimal tuning parameter
  classProbs = T,                  # should class probabilities be returned
  summaryFunction=twoClassSummary  #,  # results summary function
  # savePredictions = 'final'
) 


#################################################################################################
#### Section V, step 4
#### Tune hyper parameters by setting tuneLength

set.seed(12345)

#############################
#### Step 4.1 choose svm_linear method
model_svmLinear = train(class ~ ., 
                        data=trainData, 
                        method='svmLinear', 
                        tuneLength = 20, 
                        metric='ROC', 
                        trControl = fitControl
)

#############################
## briefly check the svmLinear results
model_svmLinear 
# model_svmLinear$pred$yes
# model_svmLinear$pred$no


varimp_svmLinear <- varImp(model_svmLinear)
varimp_svmLinear

varImp( model_svmLinear )$importance

plot(varimp_svmLinear, main="Variable Importance with svmLinear")

########################################################################################
## step 4.2
## plot roc for the training data
## We can calculate the area under the curve...
## Select a parameter setting
## selectedIndices <- model_mars2$pred

rocobj_svmlinear <- roc(model_svmLinear$pred$obs, model_svmLinear$pred$yes, ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE,
                        main="svmLinear",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=5, 
                        print.auc=TRUE)


############################

###############################################################
## Step 4.3  Predict on testData and Compute the confusion matrix
## 

predicted2 <- predict(model_svmLinear, testData)

################################################################
## print out confusion matrix 

confusionMatrix(reference = testData$class, data = predicted2, mode='everything', positive='yes')

#######################################
# Confusion Matrix and Statistics
#
#             Reference
# Prediction  No Yes
#         No   7   1
#         Yes  0   2
# #####################################
# Accuracy : 0.9            
# 95% CI : (0.555, 0.9975)
# No Information Rate : 0.7            
# P-Value [Acc > NIR] : 0.1493  
# 
### save the ROC with svmRadiacl plot to local drive
path <- "C:/Users/dug/SH_Hospital/SVM_RFE_Prj/0507CleanedData/Plots/"
path <- paste0(getwd(), "/60vs40Plots/" )
path

PlotROCSaveJPG(rocobj_svmlinear, path, "svmLinear_tune18_train60", "svmLinear" )


################################################################
## print out the whole prediction results
predicted2


######################################################################
#### Section V, step 5
#### Hyper Parameter Tuning using tuneGrid
#### Alternately, you can set the tuneGrid instead of tuneLength.

## Step 5.1: Define the tuneGrid
marsGrid <-  expand.grid(nprune = c(1:10), 
                         degree = c(1:3))

## Step 5.2: Tune hyper parameters by setting tuneGrid
set.seed(123)
model_marsG = train(class ~ ., 
                    data=trainData, 
                    method='earth', 
                    metric='ROC', 
                    tuneGrid = marsGrid, 
                    trControl = fitControl
)

## check tuned MARS model results
model_marsG

## Step 5.2 plot roc for tuned MARS model
rocobj_marsG<- roc(model_marsG$pred$obs, 
                   model_marsG$pred$yes, ci=TRUE,
                   plot=TRUE, 
                   legacy.axes=TRUE, percent=TRUE, 
                   main="MARS_with_Grid",
                   xlab="False Positive Percentage", 
                   ylab="True Postive Percentage", 
                   col="darkblue", lwd=4, 
                   print.auc=TRUE)

## Step 5.3: Predict on testData and Compute the confusion matrix 
predictedG <- predict(model_marsG, testData)
confusionMatrix(reference = testData$class, data = predictedG, mode='everything', positive='yes')

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_marsG, path, "MARS_Grid_failed", "MARS_Grid" )


######################################################################
#### Section V, Step 6  
#### use random forest method
####
###################################################################### 
## the gbm algo does not work well on this 'small' training set;
## gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 5), 
##                        n.trees = (1:5)*5, 
##                        shrinkage = 0.1,
##                        n.minobsinnode = 15)

##################################################
### Step 6.1 train with RF
set.seed(100)
model_rf = train(class ~ ., 
                 data=trainData, 
                 method='rf', 
                 metric='ROC', 
                 #tuneGrid = gbmGrid, 
                 tuneLength = 1, 
                 trControl = fitControl
)

### Step 6.2 plot ROC 
rocobj_rf<- roc(model_rf$pred$obs, model_rf$pred$yes, 
                ci=TRUE,
                plot=TRUE, 
                legacy.axes=TRUE, percent=TRUE,
                main="Random Forest",
                xlab="False Positive Percentage", 
                ylab="True Postive Percentage", 
                col="darkblue", lwd=4, 
                print.auc=TRUE)



### Step 6.3 Predict on testData and Compute the confusion matrix
predicted_rf <- predict(model_rf, testData)
confusionMatrix(reference = testData$class, data = predicted_rf, mode='everything', positive='yes')

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_rf, path, "RandomForest_tune1", "Random Forest" )



########################################################################################
#### Section V, step 7 Training SVM_Radial
#### based on the few training methods used, svm_Radial is the optimized one
############################################################

set.seed(100) 

### step 7.1 Train the model using SVMRadial
model_svmRadial = train(class ~ ., 
                        data=trainData, 
                        method='svmRadial', 
                        tuneLength=7,
                        trControl = fitControl
)

## check the prediction
## model_svmRadial$pred

########################################################
## step #7.2 create one roc object, with AUC and 95% CI

rocobj_svmRadial <- roc(model_svmRadial$pred$obs, 
                        model_svmRadial$pred$yes, 
                        ci=TRUE,
                        plot=TRUE, 
                        legacy.axes=TRUE, percent=TRUE, 
                        main="svmRadial",
                        xlab="False Positive Percentage", 
                        ylab="True Postive Percentage", 
                        col="darkblue", lwd=4, 
                        print.auc=TRUE)

## plot ROC 

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_svmRadial, path, "svmRadical_tune5", "svmRadical" )


###################################################################
## Step 7.3: Predict on testData and Compute the confusion matrix
predict_svmRadical <- predict(model_svmRadial, testData)
confusionMatrix(reference = testData$class, data = predict_svmRadical, mode='everything', positive='yes') 


#######################################################
#      Confusion Matrix and Statistics
#      
#                   Reference
#      Prediction   no yes
#             no    6   0
#             yes   1   3
#######################################################





#############################################
#### Section V, step 8 Training KNN

set.seed(100)

## step 8.1 Train the model using KNN
model_knn = train(class ~ ., 
                  data=trainData, 
                  method='knn', 
                  tuneLength=4, 
                  trControl = fitControl
)

## check the prediction
# model_knn$pred
model_knn$result

############################################################
## step #8.2 plot roc
rocobj_knn <- roc(model_knn$pred$obs, 
                  model_knn$pred$yes, 
                  ci=TRUE,
                  plot=TRUE, 
                  legacy.axes=TRUE, percent=TRUE, 
                  main="KNN",
                  xlab="False Positive Percentage", 
                  ylab="True Postive Percentage", 
                  col="darkblue", lwd=4, 
                  print.auc=TRUE
)


############################################################
## Step 8.3: Predict on testData and Compute the confusion matrix
predict_knn <- predict(model_knn, testData)
confusionMatrix(reference = testData$class, data = predict_knn, mode='everything', positive='yes') 

### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_knn, path, "KNN_tune4", "KNN" )




#############################################
#### Section V, step 9 Training adaBoost & xgbDART

set.seed(100)

##########################################3
## Unfortunately, adaboost failed 
## Train the model using adaboost
## 
## got an error: Error in { : task 1 failed - "object 'X.HLA.DPB1.' not found"
## hopythesis, the gene name HLA-DPB1 could not be parssed as table header, HLA_DPB1 works;

## change HLA-DPB1 column header into HLA_DPB1:
## trainData$HLA_DPB1 <- trainData$`HLA-DPB1`
## trainData$`HLA-DPB1` <- NULL

## 
model_ada = train(class ~ ., 
                  data=trainData, 
                  method='ada', 
                  tuneLength=2, 
                  trControl = fitControl
) 


############################################################
## step #9.2 plot roc for ada
## this step took about 90 mins

rocobj_ada <- roc(model_ada$pred$obs, 
                  model_ada$pred$yes, 
                  ci=TRUE,
                  plot=TRUE, 
                  legacy.axes=TRUE, percent=TRUE, 
                  main="ADA",
                  xlab="False Positive Percentage", 
                  ylab="True Postive Percentage", 
                  col="darkblue", lwd=4, 
                  print.auc=TRUE
)


############################################################
## Step #9.3: Predict on testData and Compute the confusion matrix
# testData$HLA_DPB1 <- testData$`HLA-DPB1`
# testData$`HLA-DPB1`<- NULL

predict_ada <- predict(model_ada, testData)
confusionMatrix(reference = testData$class, data = predict_ada, mode='everything', positive='yes') 


### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_ada, path, "ada_turne2", "ADA" )




#############################################
## Section V Step 10 train with sbgDART
## 
#############################################
#### Step #10. 1 Train the model using xgbDART 
#### this step might take 60+ mins for a small input data matrix
set.seed(100)

model_xgbDART = train(class ~ ., 
                      data=trainData, 
                      method='xgbDART', 
                      tuneLength=2, 
                      trControl = fitControl, 
                      verbose=F
)

model_xgbDART


############################################################
## step #9.2 plot roc
rocobj_xgbDART <- roc(model_xgbDART$pred$obs, 
                      model_xgbDART$pred$yes, 
                      ci=TRUE,
                      plot=TRUE, 
                      legacy.axes=TRUE, percent=TRUE, 
                      xlab="False Positive Percentage", 
                      ylab="True Postive Percentage", 
                      col="darkblue", lwd=4, 
                      print.auc=TRUE
)


############################################################
## Step #9.3: Predict on testData and Compute the confusion matrix
predict_xgbDART <- predict(model_xgbDART, testData)
confusionMatrix(reference = testData$class, data = predict_xgbDART, mode='everything', positive='yes') 


### save the ROC with svmRadiacl plot to local drive
PlotROCSaveJPG(rocobj_xgbDART, path, "xgbDART", "xgbDART" )






#####################################################################################
#### Section VI: Compare models
#### call resample() function



# # Compare model performances using resample()

models_compare <- resamples(
  
  list(RF=model_rf, XGBDART=model_xgbDART,   
       SVM=model_svmRadial, SVML = model_svmLinear, knn = model_knn
  )
  
)

model_ada$resample        #150
model_rf$resample         #10
model_xgbDART$resample    #10
model_mars$resample       #25
model_marsG$resample      #150
model_svmLinear$resample  #10
model_svmRadial$resample  #10 
model_knn$resample        #10


# Summary of the models performances
summary(models_compare)

# > summary(models_compare)
# 
#      Call:
#        summary.resamples(object = models_compare)
#      
#      Models: RF, XGBDART, SVM, SVML, knn 
#      Number of resamples: 10 
#      
#      ROC 
#      Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#      RF      0.00 1.00000      1 0.8833333       1    1    0
#      XGBDART 0.75 1.00000      1 0.9750000       1    1    0
#      SVM     0.50 1.00000      1 0.9000000       1    1    0
#      SVML    0.00 0.62500      1 0.8000000       1    1    0
#      knn     0.50 0.90625      1 0.9125000       1    1    0
#      
#      Sens 
#      Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#      RF       0.5   1.000   1.00 0.9000000       1    1    0
#      XGBDART  0.5   0.625   1.00 0.8500000       1    1    0
#      SVM      0.5   0.500   0.75 0.7500000       1    1    0
#      SVML     0.5   0.750   1.00 0.8666667       1    1    0
#      knn      1.0   1.000   1.00 1.0000000       1    1    0
#      
#      Spec 
#      Min. 1st Qu. Median      Mean 3rd Qu. Max. NA's
#      RF         0   0.000   1.00 0.60           1.00    1    0
#      XGBDART    0   0.625   1.00 0.75           1.00    1    0
#      SVM        0   1.000   1.00 0.85           1.00    1    0
#      SVML       0   0.000   0.75 0.55           1.00    1    0
#      knn        0   0.000   0.00 0.30           0.75    1    0
########################################################################      

## Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(models_compare, scales=scales)





##################################################################################
####  Section VII
####  ensemble predictions from multiple models using caretEnsemble
#### 

library(caretEnsemble)

# Stacking Algorithms - Run multiple algos in one call.
trainControl <- trainControl(method="repeatedcv", 
                             number=25, 
                             repeats=15,
                             savePredictions=TRUE, 
                             classProbs=TRUE)

algorithmList <- c('rf', 'knn', 'earth', 'xgbDART', 'svmRadial', 'svmLinear')

algorithmList <- c( 'rf', 'svmLinear')


################################################################################
## Run all algorithms in the list: 
set.seed(1234)

models <- caretList(class ~ ., 
                    data=trainData, 
                    trControl=trainControl, 
                    methodList=algorithmList) 


################################################################################
## check resample() results

results <- resamples(models)
summary(results)


#    > summary(results)
#    
#    Call:
#      summary.resamples(object = results)
#    
#    Models: rf, knn, earth, xgbDART, svmRadial, svmLinear 
#    Number of resamples: 30 
#    
#    Accuracy 
#    Min.   1st Qu.    Median      Mean 3rd Qu. Max. NA's
#    rf        0.3333333 0.6666667 0.7500000 0.7822222  1.0000    1    0
#    knn       0.3333333 0.6666667 0.6666667 0.7350000  0.7500    1    0
#    earth     0.0000000 0.6666667 0.6666667 0.7100000  1.0000    1    0
#    xgbDART   0.3333333 0.6666667 0.9000000 0.8211111  1.0000    1    0
#    svmRadial 0.3333333 0.6666667 0.7750000 0.8072222  1.0000    1    0
#    svmLinear 0.2500000 0.6666667 0.6666667 0.7355556  0.9375    1    0
#    
#    Kappa 
#              Min. 1st Qu.    Median      Mean 3rd Qu. Max. NA's
#    rf        -0.5     0.0 0.5000000 0.4181818   1.000    1    0
#    knn       -0.5     0.0 0.0000000 0.2681818   0.500    1    0
#    earth     -0.8     0.0 0.2000000 0.3515152   1.000    1    0
#    xgbDART   -0.5     0.0 0.7727273 0.5448485   1.000    1    0
#    svmRadial  0.0     0.4 0.5227273 0.6115152   1.000    1    0
#    svmLinear -0.5     0.0 0.2000000 0.3245455   0.875    1    0


################################################################################
## Box plots to compare models
scales <- list(x=list(relation="free"), y=list(relation="free"))
bwplot(results, scales=scales)

## Save as multi_algo_Accuracy_Kappa_boxplot



################################################################################
## Plot multi ROCs in one plot
rocobj_models <- roc(models$rf$pred$obs, 
                     models$rf$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="darkblue", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 40
)

rocobj_models <- roc(models$svmRadial$pred$obs, 
                     models$svmRadial$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="green", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 44,
                     add = TRUE
)


rocobj_models <- roc(models$svmLinear$pred$obs, 
                     models$svmLinear$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="red", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 48,
                     add = TRUE
)

rocobj_models <- roc(models$xgbDART$pred$obs, 
                     models$xgbDART$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="black", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 52,
                     add = TRUE
)

rocobj_models <- roc(models$earth$pred$obs, 
                     models$earth$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="yellow", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 56,
                     add = TRUE
)


rocobj_models <- roc(models$knn$pred$obs, 
                     models$knn$pred$yes, 
                     ci=TRUE,
                     plot=TRUE, 
                     legacy.axes=TRUE, percent=TRUE, 
                     xlab="False Positive Percentage", 
                     ylab="True Postive Percentage", 
                     col="pink", lwd=4, 
                     print.auc=TRUE,
                     print.auc.y = 35.1,
                     add = TRUE
)


legend("bottomright", 
       legend=c( "rf", "svmRadial", "svmLinear", "xgbDART", "MARS", "knn" ), 
       col=c( "darkblue", "green", "red", "black", "yellow", "pink" ), 
       lwd=4
)


########################################################################
## Run all algorithms in the list: against validation dataset

set.seed(100)

valid.models <- caretList(class ~ ., 
                          data=testData, 
                          trControl=trainControl, 
                          #trControl = fitControl,
                          methodList=algorithmList) 



###################3
## Plot two ROC curves together, 
## compare ROC train and ROC validation, although this is a pretty odd comparison
## 


rocobj_models_svmLinear <- roc(models$svmLinear$pred$obs, 
                               models$svmLinear$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               main = "svm-Linear ROC-AUC Expression",
                               xlab=" 1 - Specificity", 
                               ylab="Sensitivity", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$svmLinear$pred$obs, 
                                 valid.models$svmLinear$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 main = "svm-Linear ROC-AUC Expression",
                                 xlab="1 - Specifiity", 
                                 ylab="Sensitivity", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "svm Linear Training", "svm Linear Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)

## the variable importance would change a bit if we use different seeds; 
## for example, KLRC1 and IRF4 contricute almost equally to the model, 
## so the rank of those two genes might switch if we use another seed. 
plot( varImp(models$svmLinear) ) 


frame()

rocobj_models_svmRadial <- roc(models$svmRadial$pred$obs, 
                               models$svmRadial$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$svmRadial$pred$obs, 
                                 valid.models$svmRadial$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "svm Radial Training", "svm Radial Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)





frame()

rocobj_models_svmLinear <- roc(models$rf$pred$obs, 
                               models$rf$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$rf$pred$obs, 
                                 valid.models$rf$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "Random Forest Training", "Random Forest Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)





### KNN
frame()

rocobj_models_ <- roc(models$knn$pred$obs, 
                      models$knn$pred$yes, 
                      ci=TRUE,
                      plot=TRUE, 
                      legacy.axes=TRUE, percent=TRUE, 
                      xlab="False Positive Percentage", 
                      ylab="True Postive Percentage", 
                      col="darkblue", lwd=4, 
                      print.auc=TRUE,
                      print.auc.y = 35.1
                      #add = TRUE
)

rocobj_models_knn <- roc(valid.models$knn$pred$obs, 
                         valid.models$knn$pred$yes, 
                         ci=TRUE,
                         plot=TRUE, 
                         legacy.axes=TRUE, percent=TRUE, 
                         xlab="False Positive Percentage", 
                         ylab="True Postive Percentage", 
                         col="red", lwd=4, 
                         print.auc=TRUE,
                         print.auc.y = 25.1,
                         add = TRUE
)



legend("bottomright", 
       legend=c( "KNN Training", "KNN Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)



### xgbDART
frame()

rocobj_models_svmLinear <- roc(models$xgbDART$pred$obs, 
                               models$xgbDART$pred$yes, 
                               ci=TRUE,
                               plot=TRUE, 
                               legacy.axes=TRUE, percent=TRUE, 
                               xlab="False Positive Percentage", 
                               ylab="True Postive Percentage", 
                               col="darkblue", lwd=4, 
                               print.auc=TRUE,
                               print.auc.y = 35.1
                               #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$xgbDART$pred$obs, 
                                 valid.models$xgbDART$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "xgbDART Training", "xgbDART Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)


## MARS
frame()

rocobj_models_mars <- roc(models$earth$pred$obs, 
                          models$earth$pred$yes, 
                          ci=TRUE,
                          plot=TRUE, 
                          legacy.axes=TRUE, percent=TRUE, 
                          xlab="False Positive Percentage", 
                          ylab="True Postive Percentage", 
                          col="darkblue", lwd=4, 
                          print.auc=TRUE,
                          print.auc.y = 35.1
                          #add = TRUE
)

rocobj_models_svmLinearVa <- roc(valid.models$earth$pred$obs, 
                                 valid.models$earth$pred$yes, 
                                 ci=TRUE,
                                 plot=TRUE, 
                                 legacy.axes=TRUE, percent=TRUE, 
                                 xlab="False Positive Percentage", 
                                 ylab="True Postive Percentage", 
                                 col="red", lwd=4, 
                                 print.auc=TRUE,
                                 print.auc.y = 25.1,
                                 add = TRUE
)



legend("bottomright", 
       legend=c( "MARS Training", "MARS Validation" ), 
       col=c( "darkblue", "red" ), 
       lwd=4
)



################################################################################
####  Section VIII 
####  combine the predictions of multiple models to form a final prediction


# Create the trainControl
set.seed(100)

stackControl <- trainControl(method="repeatedcv", 
                             number=20, 
                             repeats=5,
                             savePredictions=TRUE, 
                             classProbs=TRUE
)

# Ensemble the predictions of `models` to form a new combined prediction based on glm
stack.glm <- caretStack(models, method="glm", metric="Accuracy", trControl=stackControl)

print(stack.glm)

# Predict on testData 
stack_predicteds <- predict(stack.glm, newdata=testData) 
head(stack_predicteds) 

confusionMatrix(reference = testData$class, data = stack_predicteds, mode='everything', positive='Yes') 


########################################################################3
####
####         END OF THE MAIN CODE SECTION                           ####3
####
########################################################################3





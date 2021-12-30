#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

# Neural Network Model:
# breastcancer cell severity classification with nnet function

abc <- read.csv("breast_canceer.csv")
head(abc)
tail(abc)

# column names made shorter for cabconvenint plotting of nnet model

# CPTH  clump.thick       CLSZ  cell.size 
# CLSP  cell.shape        MADH  mar.adhesion
# EPTH  epithelial        BNUC  bare.nuclei
# BLCR  bland.chromatin   NNUC  normal.nuclei 
# MTSE  mitoses           class Tumur class

colnames(abc) <- c("patient.id", "CPTH", "CLSZ", "CLSP", "MADH", "EPTH", "BNUC", 
                   "BLCR", "NNUC", "MTSE", "class")

head(abc)
tail(abc)

#----------------------------------------------

# normalization of dataset
# two methods of normalization: MinMax metod and z-score 
# Let us use MinMax method in this example

# proportion of 'benign' and 'severe' class in 'abc'
round(prop.table(table(abc$class)),2)

minval <- apply(abc[,2:11], 2, min)
maxval <- apply(abc[,2:11], 2, max)
scaledf <- scale(abc[,2:11], center=minval, scale=maxval-minval)

bcdata <- data.frame(scaledf)
head(bcdata)

#----------------------------------------------

# traindata and testdata partioning with two indices

index <- sample(2, nrow(bcdata), replace=T, prob=c(0.7,0.3))

traindata <- bcdata[index==1,]
testdata <- bcdata[index==2,]

head(traindata)

# number of instance in traindata
nrow(traindata)

# number of 'benign' and 'severe' class in traindata
table(traindata$class)

# proportion of 'benign' and 'severe' class cancer cell
round(prop.table(table(traindata$class)),2)


head(testdata)

# number of instance in testdata
nrow(testdata)

# number of 'benign' and 'severe' class in testdata
table(testdata$class)

# proportion of 'benign' and 'severe' class cancer cell
round(prop.table(table(testdata$class)),2)


#----------------------------------------------

library(nnet)
library(NeuralNetTools)

nnfit <- nnet(class~., data=traindata, size=5, decay=1e-03, linout=F, 
              skip=F, maxit=1000, Hess=T)

summary(nnfit)

plotnet(nnfit, cex=0.7)

legend("bottomright", legend=c("CPTH  clump.thick", "CLSZ  cell.size", 
                               "CLSP  cell.shape", "MADH mar.adhesion", "EPTH  epithelial",
                               "BNUC  bare.nuclei", "BLCR  bland.chromatin", "NNUC normal.nuclei", 
                               "MTSE  mitoses", "class   Tumur class"), bty="n", cex=0.6)


#----------------------------------------------

# fitted values for breastcancer with traindata

fitval <- fitted.values(nnfit)
head(fitval)

fitval <- ifelse(fitval > 0.5, 1, 0)
result1 <- data.frame(classN = bcdata$class[index==1], fitval=fitval)
head(result1)

# total number of instances in traindata
nrow(result1)

# confusion matrix for traindata
tb1 <- table(result1)
tb1

# correct classification for traindata:
sum(diag(tb1))

# misclassification for traindata:
misclass1 <- sum(tb1) - sum(diag(tb1))
misclass1

# classification accuracy for traindata
accuracy1 = sum(diag(tb1))/sum(tb1)
accuracy1

# classification error for traindata
error1 <- 1-sum(diag(tb1))/sum(tb1)
error1

#----------------------------------------------

# misclassified instances for traindata with nnet model
misfitLoc <- as.integer(rownames(result1)[result1$classN != result1$fitval])

print(misfitLoc)

# data for misclassified cancer instances & model results
misfit <- result1[(result1$classN != result1$fitval),]
cbind(abc[misfitLoc,2:11], misfit)


#----------------------------------------------

# predicted breastcancer severity for testdata

fpre <- predict(nnfit, newdata=testdata[,1:9])
predval <- fpre[,1]
head(predval)

predval <- ifelse(predval>0.5,1,0)
result2 <- data.frame(classN=bcdata$class[index==2], predval=predval)
head(result2)

# total number of instances in testdata
nrow(result2)

# confusion matrix for testdata
tb2 <- table(result2)
tb2

# correct classification for testdata:
sum(diag(tb2))

# misclassification for testdata:
misclass2 <- sum(tb2) - sum(diag(tb2))
misclass2

# classification accuracy for testdata
accuracy2 = sum(diag(tb2))/sum(tb2)
accuracy2

# classification error for testdata
error2 <- 1-sum(diag(tb2))/sum(tb2)
error2

# misclassified instances for testdata with nnet model
mispredLoc <- as.integer(rownames(result2)[result2$classN != result2$predval])

print(mispredLoc)

# data for misclassified cancer instances & model results
mispred <- result2[(result2$classN != result2$predval),]
cbind(abc[mispredLoc,2:11], mispred)


# matrix plot of the breast cancer cell data
pairs(abc[,2:11], col=abc[,11]+2, pch=(abc[,11]+16))

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------



plot(abc$CPTH, abc$CLSZ, col=abc[,11]+2, pch=(abc[,11]+16), 
     xlab = "clump thickness", ylab = "cell size") 

library(ggplot2)

#Scatter plot: Sepal Length versus Petal Width 
ggplot(abc, aes(x=CPTH, y=CLSZ)) +
  geom_point(aes(colour=factor(class), shape=factor(class)), size=2) +
  xlab("clump thickness") +
  ylab("cell size") +
  ggtitle("clth vs clsz")
## HONOURS ##
set.seed(5018090)
library(survival)
library(randomForestSRC)
library(glmnet)
library(ggplot2)
library(ggRandomForests)
library(survminer)
library(ROCR)
options(scipen = 500)


### READ IN DATA ###
rawdata <- read.csv("westpac_date.csv", header = T)
# rm(list=setdiff(ls(), "rawdata"))
data_approved <- rawdata[rawdata$Final_Decision_Summary == 'APPROVE',]   #left truncation #

#Manage rows without starting time
data_approved <- data_approved[data_approved$OPENED_DATE != '',] #removes 39367 observations without opening dates


### Process variables for modelling ###
#Change occupation codes to factors
levels(data_approved$OCCUPATION_CODE_1) = list(A = 10000:19999, B = 20000:29999, C = 30000:39999, D = 40000:49999,
                                               E = 50000:59999, F = 60000:69999, G = 70000:79999, H = 80000:89999,
                                               I = 95011, J = c(95021, 95925), K = 95031, L = c(92001, 92501),
                                               M = 98001, N = c(99999, ''))

# 9's into NAs
is.na(data_approved) <- data_approved == '99'
is.na(data_approved) <- data_approved == '999'
is.na(data_approved) <- data_approved == '9999'
is.na(data_approved) <- data_approved == '99999'
is.na(data_approved) <- data_approved == '999999'
is.na(data_approved) <- data_approved == '9999999'
is.na(data_approved) <- data_approved == '99999999'
is.na(data_approved) <- data_approved == '999999999'
is.na(data_approved) <- data_approved == '9999999999'
is.na(data_approved) <- data_approved == '99999999999'
is.na(data_approved) <- data_approved == ''


### PRE-PROCESS DATA ###
data_approved$OPENED_DATE<-as.Date(data_approved$OPENED_DATE, format = "%d-%m-%Y")              #converting factors to dates for 
data_approved$CLOSED_DATE<-as.Date(data_approved$CLOSED_DATE, format = "%d-%m-%Y")              #calculation #
data_approved$DEFAULT_MONTH<-as.Date(data_approved$DEFAULT_MONTH, format = "%d-%m-%Y")

#calculate censoring and death times
T_i <- as.numeric(difftime(data_approved$DEFAULT_MONTH, data_approved$OPENED_DATE, units = "days"))
C_i <- as.numeric(difftime(data_approved$CLOSED_DATE, data_approved$OPENED_DATE, units = "days"))
D_i <- data_approved$DEFAULT_FLAG

T_i[is.na(T_i)] = C_i[is.na(T_i)]


## REMOVING VARIABLES ##
#calculate percentage of missing data in each column
percent_missing <- function(data, n){
  percent_col_NA = c(seq(1, ncol(data)))
  for(i in 1:ncol(data)){
    percent_col_NA[i] = sum(is.na(data[, i]))/nrow(data)*100
  }
  colnames(data[,c(which(percent_col_NA >= n))])
}
percent_missing(data_approved, 25)


#remove unncessary variables
rmv1 <- c("Final_Decision_Summary", "ORIG_DECSN")
data_approved_1 <- data_approved[ , !(names(data_approved) %in% rmv1)]

#remove variables with high predictive power for default
rmv2 <- c("CLOSED_DATE", "DEFAULT_AMOUNT", "DEFAULT_MONTH", "DEFAULT_FLAG")
data_approved_1 <- data_approved_1[ , !(names(data_approved_1) %in% rmv2)]

#checking again
percent_missing(data_approved_1, 25)

sum(is.na(data_approved_1))/(nrow(data_approved_1)*ncol(data_approved_1))*100 
#percentage of NAs in entire data_approved 1 set = 9.84686% before "remove highly missing data"

#remove highly missing data
rmv3 <- percent_missing(data_approved_1, 25)
data_approved_1 <- data_approved_1[ , !(names(data_approved_1) %in% rmv3)]

#remove irrelevant variables
rmv4 <- c("APPL_MTH", "Customer_ID_anon", "ACCOUNT_NUMBER_anon", "APPL_NUM_anon")
data_approved_1 <- data_approved_1[ , !(names(data_approved_1) %in% rmv4)]

### SECOND PRE-PROCESS ###
data_approved_1$OPENED_DATE <- as.integer(format(data_approved_1$OPENED_DATE, "%Y%m%d"))

data_approved_1$APPL_DATE <- as.character(data_approved_1$APPL_DATE)
data_approved_1$APPL_DATE <- substr(data_approved_1$APPL_DATE,1,nchar(data_approved_1$APPL_DATE)-9)
data_approved_1$APPL_DATE <- as.Date.character(data_approved_1$APPL_DATE, format = "%d%b%Y")
data_approved_1$APPL_DATE <- as.integer(format(data_approved_1$APPL_DATE, "%Y%m%d"))

data_approved_1$Surv_time <- T_i
data_approved_1$Event <- D_i

data_final <- data_approved_1
rm(data_approved_1)

### IMPUTE DATA
sum(is.na(data_final))/prod(dim(data_final))*100 #0.6417696% NAs
data_imputed <- impute(Surv(Surv_time, Event)~., data = data_final)
sum(is.na(data_imputed))/prod(dim(data_imputed))*100 # equal to zero

#Example cox prop haz model
# default.coxph <- coxph(Surv(T_i, D_i)~CR_FUNDS_WITH_WBC_1 + DRIVERS_LICENCE_IND_1 + CURR_RESIDENCY_POSTCODE_1, rawdata_approved)

# percent_col_NA = c(seq(1, ncol(data_approved_1)))
# for(i in 1:ncol(data_approved_1)){
#   percent_col_NA[i] = sum(is.na(data_approved_1[, i]))/nrow(data_approved_1)*100
# }
# percent_col_NA
# colnames(data_approved_1[,which(percent_col_NA > 1)])


### IMBALANCED DATA AND COMPUTATIONAL LIMIT
#Down sampling
prop <- length(which(data_final$Event == 1))/length(data_final$Event)*100 # 5.722146% is default
data_imputed_def <- data_imputed[data_imputed$Event == 1,] #split into default
data_imputed_non_def <- data_imputed[data_imputed$Event == 0,] #split into non default

subsample_imputed_non_def <- data_imputed_non_def[sample(nrow(data_imputed_non_def), 40000, replace = FALSE),]

data_imputed_subsamp <- rbind(data_imputed_def, subsample_imputed_non_def)

#Stratified sampling
n_1 <- nrow(data_imputed_non_def)
stratsamp_imputed_non_def <- data_imputed_non_def[sample(c(1:n_1), n_1*0.32, replace = FALSE),]
n_2 <- nrow(data_imputed_def)
stratsamp_imputed_def <- data_imputed_def[sample(c(1:n_2), n_2*0.32, replace = FALSE),]

data_imputed_stratsamp <- rbind(stratsamp_imputed_def, stratsamp_imputed_non_def)

###### MODELLING AND STATISTICS #####
set.seed(5018090)
##### Random survival forests #####

##### Original whole sample
loan.rsf_1 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed, splitrule = "random",
                    ntree = 100, nsplit = 10, block.size = 5, do.trace = 5) #split rule random for fast
plot.survival(loan.rsf_1)
loan.cox_1 <- coxph(Surv(Surv_time, Event)~., data = data_imputed)
summary(loan.cox_1)

# find the most valuable variables
loan.rsf_1_maxsub <- max.subtree(loan.rsf_1)
min.depth <- loan.rsf_1_maxsub$order[, 1]
min.depth_best <- which(min.depth < loan.rsf_1_maxsub$threshold)

data_imputed_1 <- subset(data_imputed, select = c(names(min.depth_best), "Event", "Surv_time")) #removed 10 variables

#Plot top 2 and bottom 2 variables in terms of minimal depth
name_list <- colnames(data_imputed[,c(10,54,86,97)])
plot.variable(loan.rsf_1, xvar.names = name_list)
#RSF model with 10 removed from min depth criteria
loan.rsf_1_best <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed_1, splitrule = "random",
                         ntree = 100, nsplit = 10, block.size = 5, do.trace = 5)  
#####BEST 7 VARIABLES IN TERMS OF MIN DEPTH
min.depth <- loan.rsf_1_maxsub$order[, 1]
min.depth_best_7 <- which(min.depth < 6.39)
best_7 <- c("LN_AMT_1","LN_AMT_2", "NUM_BUR_ENQ_LST_7TO12_MTHS_1","OCCUPATION_CODE_1",
            "RESIDNTL_STATUS_1","TIME_AT_PREV_ADDR_1","MTH_CC_PAY_1")
data_imputed_B7 <- subset(data_imputed, select = c(best_7, "Event", "Surv_time")) 


#BEST 7
loan.rsf_2 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed_B7, splitrule = "random",
                    ntree = 100, nsplit = 10, block.size = 5, do.trace = 5)
loan.rsf_2

loan.cox_2 <- coxph(Surv(Surv_time, Event)~., data = data_imputed_B7)
summary(loan.cox_2)

##### STRATIFIED SAMPLE
loan.rsf_stratsamp_1 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed_stratsamp, splitrule = "random",
                              ntree = 100, nsplit = 10, block.size = 5, do.trace = 5)
loan.rsf_stratsamp_1

loan.cox_stratsamp_1 <- coxph(Surv(Surv_time, Event)~., data = data_imputed_stratsamp)
summary(loan.cox_stratsamp_1)

data_imputed_stratsamp_B7 <- subset(data_imputed_stratsamp, select = c(best_7, "Event", "Surv_time"))
loan.rsf_stratsamp_2 <-rfsrc(Surv(Surv_time, Event)~., data = data_imputed_stratsamp_B7, splitrule = "random",
                             ntree = 100, nsplit = 10, block.size = 5, do.trace = 5)
loan.rsf_stratsamp_2

loan.cox_stratsamp_2 <- coxph(Surv(Surv_time, Event)~., data = data_imputed_stratsamp_B7)
summary(loan.cox_stratsamp_2)


##### SUBSAMP: downsampling non-defaults
loan.rsf_subsamp_1 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed_subsamp, splitrule = "random",
                            ntree = 100, nsplit = 10, block.size = 5) #using 40000 non-defaults
loan.rsf_subsamp_1
loan.cox_subsamp_1 <- coxph(Surv(Surv_time, Event)~., data = data_imputed_subsamp)
#OOB Brier Score progressively got worst with time compared to regular rsf model

loan.rsf_subsamp_1_maxsub <- max.subtree(loan.rsf_subsamp_1)
min.depth_subsamp <- loan.rsf_subsamp_1_maxsub$order[, 1]
min.depth_subsamp_best <- which(min.depth_subsamp < loan.rsf_subsamp_1_maxsub$threshold)

data_imputed_subsamp_B7 <- subset(data_imputed_subsamp, select = c(best_7, "Event", "Surv_time")) #removed 11 variables

loan.rsf_subsamp_2 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed_subsamp_B7, splitrule = "random",
                            ntree = 100, nsplit = 10, block.size = 5, do.trace = 5)
loan.rsf_subsamp_2

loan.cox_subsamp_2 <- coxph(Surv(Surv_time, Event)~., data = data_imputed_subsamp_B7)
summary(loan.cox_subsamp_2)

##### KM Estimator
par(mfrow = c(1,2))
loan.KM_no_subsamp <- survfit(Surv(Surv_time, Event)~1, data = data_imputed)
loan.KM_w_subsamp <- survfit(Surv(Surv_time, Event)~1, data = data_imputed_subsamp)
plot(loan.KM_no_subsamp)
plot(loan.KM_w_subsamp, main = "Kaplan-Meier estimate", ylab = "Probabality to default past t",
     xlab = "Time, t", xlim=c(0,550))

##### Calculate the VIMP for 30 variables and error rate and plot
loan.rsf_stratsamp_2_30 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed_2_30, splitrule = "random",
                                 ntree = 100, nsplit = 10, block.size = 5, do.trace = 5)
vimp.loan.rsf_stratsamp_2_30 <- vimp(loan.rsf_stratsamp_2_30, do.trace = 5)

##### Logistic regression error - function
R = 10
n = nrow(data)

# empty Rx2 matrix for bootstrap results 
B = matrix(nrow = R, ncol = 2, dimnames = list(paste('Sample',1:R), c("auc_orig","auc_boot")))

set.seed(5018090)
for(i in 1:R){
  
  # draw a random sample
  obs.boot <- sample(x = 1:n, size = n, replace = T)
  data.boot <- data[obs.boot, ]
  
  # fit the model on bootstrap sample
  logit.boot <- glm(Event~. , data = data.boot, family = binomial(link = "logit"))
  
  # apply model to original data
  prob1 = predict(logit.boot, type='response', data)
  pred1 = prediction(prob1, data$Event)
  auc1 = performance(pred1,"auc")@y.values[[1]][1]
  B[i, 1] = auc1
  
  # apply model to bootstrap data
  prob2 = predict(logit.boot, type='response', data_imputed_stratsamp_B7.boot)
  pred2 = prediction(prob2, data.boot$Event)
  auc2 = performance(pred2,"auc")@y.values[[1]][1]
  B[i, 2] = auc2
}
mean(B[,1])

## HONOURS ##
set.seed(5018090)
library(survival)
library(randomForestSRC)
library(parallel)
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
percent_missing(data_approved, 40)

#remove unncessary variables
rmv1 <- c("Final_Decision_Summary", "ORIG_DECSN")
data_approved_1 <- data_approved[ , !(names(data_approved) %in% rmv1)]

#remove variables with high predictive power for default
rmv2 <- c("CLOSED_DATE", "DEFAULT_AMOUNT", "DEFAULT_MONTH", "DEFAULT_FLAG")
data_approved_1 <- data_approved_1[ , !(names(data_approved_1) %in% rmv2)]

#checking again
percent_missing(data_approved_1, 40)

#remove highly missing data
rmv3 <- percent_missing(data_approved_1, 40)
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
#Subsampling
prop <- length(which(data_final$Event == 1))/length(data_final$Event)*100 # 5.722146% is default
data_imputed_def <- data_imputed[data_imputed$Event == 1,] #split into default
data_imputed_non_def <- data_imputed[data_imputed$Event == 0,] #split into non default

subsample_imputed_non_def <- data_imputed_non_def[sample(nrow(data_imputed_non_def), 40000, replace = FALSE),]

data_imputed_subsamp <- rbind(data_imputed_def, subsample_imputed_non_def)


### MODELLING AND STATISTICS
set.seed(5018090)
## Random survival forests

#original whole sample
loan.rsf_1 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed, splitrule = "random",
                    ntree = 100, na.action = "na.omit", nsplit = 10, block.size = 100) #split rule random for fast
plot.survival(loan.rsf_1)
loan.rsf_1_maxsub <- max.subtree(loan.rsf_1)
# find the most valuable variables
min.depth <- loan.rsf_1_maxsub$order[, 1]
min.depth_best <- which(min.depth < loan.rsf_1_maxsub$threshold)

data_imputed_1 <- subset(data_imputed, select = c(names(min.depth_best), "Event", "Surv_time")) #removed 10 variables

vimp(loan.rsf_1, do.trace = )

#downsampling non-defaults
loan.rsf_subsamp_1 <- rfsrc(Surv(Surv_time, Event)~., data = data_imputed_subsamp, splitrule = "random",
                            ntree = 100, na.action = "na.omit", nsplit = 10, block.size = 100) #using 40000 non-defaults
plot.survival(loan.rsf_subsamp_1)
#OOB Brier Score progressively got worst with time compared to regular rsf model




options(scipen = 10)
loan.cox_1 <- coxph(Surv(Surv_time, Event)~., data = data_imputed_subsamp)
summary(loan.cox_1)
options(scipen = 500)

#KM Estimator
par(mfrow = c(1,2))
loan.KM_no_subsamp <- survfit(Surv(Surv_time, Event)~1, data = data_imputed)
loan.KM_w_subsamp <- survfit(Surv(Surv_time, Event)~1, data = data_imputed_subsamp)
plot(loan.KM_no_subsamp); plot(loan.KM_w_subsamp)
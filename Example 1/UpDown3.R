df <- read.table("~/Downloads/osfstorage-archive/Example 1/UpDown3.txt", header=T, stringsAsFactors=T)
str(df)
# 'data.frame':	124 obs. of  4 variables:
# $ RT   : int  421 436 390 419 359 414 311 361 379 380 ...
# $ subj : Factor w/ 31 levels "S01","S02","S03",..: 1 1 1 1 2 2 2 2 3 3 ...
# $ resp : Factor w/ 2 levels "onehand","twohands": 1 1 2 2 1 1 2 2 1 1 ...
# $ align: Factor w/ 2 levels "A","M": 1 2 1 2 1 2 1 2 1 2 ...

df$cond <- paste(df$resp, df$align, sep="") # condition
df$cond <- factor(df$cond)
mat <- reshape(df[-c(3,4)], idvar="subj", timevar="cond", direction="wide")
rownames(mat) <- mat$subj
colnames(mat)[-1] <- levels(df$cond)
mat <- mat[,c(2,4,3,5)] # wide format

round(cov(mat)) # sample covariance
#           onehandA twohandsA onehandM twohandsM
# onehandA      3937      -207     4067      -831
# twohandsA     -207      5202     -566      5922
# onehandM      4067      -566     4780      -979
# twohandsM     -831      5922     -979      7274


#############################
## Repeated-Measures ANOVA ##
#############################

summary(aov(RT ~ (resp * align) + Error(subj / (resp * align)), df))
# Error: subj
#           Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 30 270038    9001               
# 
# Error: subj:resp
#           Df Sum Sq Mean Sq F value Pr(>F)
# resp       1  21479   21479   1.854  0.183
# Residuals 30 347540   11585               
# 
# Error: subj:align
#           Df Sum Sq Mean Sq F value   Pr(>F)    
# align      1  11536   11536   28.22 9.64e-06 ***
# Residuals 30  12262     409                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: subj:resp:align
#            Df Sum Sq Mean Sq F value Pr(>F)  
# resp:align  1   1034  1033.6   5.203 0.0298 *
# Residuals  30   5960   198.7                 
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


# testing simple effects
summary(aov(RT ~ resp + Error(subj / resp), subset(df, align=="A", drop=T)))
summary(aov(RT ~ resp + Error(subj / resp), subset(df, align=="M", drop=T)))
summary(aov(RT ~ align + Error(subj / align), subset(df, resp=="onehand", drop=T)))
summary(aov(RT ~ align + Error(subj / align), subset(df, resp=="twohands", drop=T)))


################################
## Current anovaBF 0.9.12-4.7 ##
################################
library(BayesFactor)

set.seed(277)
(RIO <- anovaBF(RT ~ resp * align + subj, df, whichRandom="subj", iterations=1e5, progress=F)) # RIO
# Bayes factor analysis
# --------------
# [1] resp + subj                      : 2.075678  ±0.57%
# [2] align + subj                     : 0.6778541 ±0.85%
# [3] resp + align + subj              : 1.498859  ±0.87%
# [4] resp + align + resp:align + subj : 0.4450376 ±2.57%
# 
# Against denominator:
#   RT ~ subj 
# ---
# Bayes factor type: BFlinearModel, JZS

RIO[3] / RIO[4] # testing resp:align  : 3.367937 ±2.72% BF₀₁
RIO[2] / RIO[3] # testing resp        : 0.4522468 ±1.21% BF₀₁
RIO[1] / RIO[3] # testing align       : 1.384839 ±1.04% BF₀₁


# separate Bayesian testing
set.seed(277); anovaBF(RT ~ resp + subj, subset(df, align=="A", drop=T), 
                       whichRandom="subj", iterations=1e5, progress=F) # 1.477261 ±0.54% BF₁₀
set.seed(277); anovaBF(RT ~ resp + subj, subset(df, align=="M", drop=T), 
                       whichRandom="subj", iterations=1e5, progress=F) # 0.4350769 ±0.58% BF₁₀
set.seed(277); anovaBF(RT ~ align + subj, subset(df, resp=="onehand", drop=T), 
                       whichRandom="subj", iterations=1e5, progress=F) # 9.034816 ±0.51% BF₁₀
set.seed(277); anovaBF(RT ~ align + subj, subset(df, resp=="twohands", drop=T), 
                       whichRandom="subj", iterations=1e5, progress=F) # 2513.757 ±0.52% BF₁₀


#############################
## Including Random Slopes ##
#############################

set.seed(277)
full <- lmBF(RT ~ resp * align + subj + resp:subj + align:subj, df, whichRandom="subj", iterations=1e5, progress=F) # MRE

set.seed(277)
add <- lmBF(RT ~ resp + align + subj + resp:subj + align:subj, df, whichRandom="subj", iterations=1e5, progress=F)

set.seed(277)
omitA <- lmBF(RT ~ align + subj + resp:subj + align:subj, df, whichRandom="subj", iterations=1e5, progress=F)

set.seed(277)
omitB <- lmBF(RT ~ resp + subj + resp:subj + align:subj, df, whichRandom="subj", iterations=1e5, progress=F)

# testing AB
add / full # 0.4879953 ±7.23% BF₀₁

# testing A
omitA / add # 1.891711 ±4.93% BF₀₁

# testing B
omitB / add # 0.0006079976 ±9.26% BF₀₁


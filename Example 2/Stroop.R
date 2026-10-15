## DATA 
## https://github.com/jasp-stats/jasp-desktop/blob/9c22a03b33630dedb073548a4ecbafdae1c1376c/Resources/Data%20Sets/Data%20Library/3.%20ANOVA/Stroop.csv

stroop <- read.csv("~/Downloads/osfstorage-archive/Example 2/Stroop.csv") # 2 x 3 repeated-measures design
colnames(stroop) <- c("a1b1", "a1b2", "a1b3", "a2b1", "a2b2", "a2b3") # wide format
round(cov(stroop[,c(1,4,2,5,3,6)])) # sample covariance
#      a1b1  a2b1 a1b2  a2b2 a1b3  a2b3
# a1b1 5444   106 6688  -103 4904  -363
# a2b1  106  8859 -265 11415 -594  7871
# a1b2 6688  -265 8682  -305 6090  -655
# a2b2 -103 11415 -305 16165 -622 10381
# a1b3 4904  -594 6090  -622 4933  -981       # A: PT (2 levels)
# a2b3 -363  7871 -655 10381 -981  7590       # B: Congruency (3 levels)


df <- data.frame("DV"=c(as.matrix(stroop)),
                 "ID"=rep(paste0("s",1:19), 2*3),
                 "A"=rep(paste0("a",1:2), each=19*3),
                 "B"=rep(paste0("b",1:3), each=19, time=2),
                 stringsAsFactors=T) # long format


#############################
## Repeated-Measures ANOVA ##
#############################

ez::ezANOVA(df, DV, ID, .(A,B), type=3)
# $ANOVA
#   Effect DFn DFd        F            p p<.05         ges
# 2      A   1  18  2.24250 1.515947e-01       0.058252810
# 3      B   2  36 22.15839 5.332514e-07     * 0.034993420
# 4    A:B   2  36  2.29659 1.151568e-01       0.003348043
# 
# $`Mauchly's Test for Sphericity`
#   Effect         W           p p<.05
# 3      B 0.6767876 0.036211394     *
# 4    A:B 0.4987519 0.002704077     *
# 
# $`Sphericity Corrections`
#   Effect       GGe        p[GG] p[GG]<.05       HFe        p[HF] p[HF]<.05
# 3      B 0.7557366 8.938571e-06         * 0.8101994 4.758814e-06         *
# 4    A:B 0.6661124 1.366417e-01           0.6993217 1.344646e-01     

summary(aov(DV ~ (A * B) + Error(ID / (A * B)), df))
# Error: ID
#           Df Sum Sq Mean Sq F value Pr(>F)
# Residuals 18 416417   23134               
# 
# Error: ID:A
#           Df Sum Sq Mean Sq F value Pr(>F)
# A          1  57533   57533   2.242  0.152
# Residuals 18 461801   25656               
# 
# Error: ID:B
#           Df Sum Sq Mean Sq F value   Pr(>F)    
# B          2  33728   16864   22.16 5.33e-07 ***
# Residuals 36  27398     761                     
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Error: ID:A:B
#           Df Sum Sq Mean Sq F value Pr(>F)
# A:B        2   3124  1562.2   2.297  0.115
# Residuals 36  24489   680.2             


################################
## Current anovaBF 0.9.12-4.7 ##
################################
library(BayesFactor)

set.seed(277)
(RIO <- anovaBF(DV ~ A * B + ID, df, whichRandom="ID", iterations=1e5, progress=F)) # RIO
# Bayes factor analysis
# --------------
# [1] A + ID           : 14.78687  ±0.55%
# [2] B + ID           : 0.7821125 ±0.46%
# [3] A + B + ID       : 14.48497  ±0.71%
# [4] A + B + A:B + ID : 2.504408  ±1.42%
# 
# Against denominator:
#   DV ~ ID 
# ---
# Bayes factor type: BFlinearModel, JZS

RIO[3] / RIO[4] # testing A:B      : 5.78379 ±1.59% BF₀₁
RIO[2] / RIO[3] # testing A        : 0.05399476 ±0.84% BF₀₁
RIO[1] / RIO[3] # testing B        : 1.020843 ±0.9% BF₀₁


#############################
## Including Random Slopes ##
#############################

set.seed(277)
full <- lmBF(DV ~ A * B + ID + A:ID + B:ID, df, whichRandom="ID", iterations=1e5, progress=F) # MRE

set.seed(277)
add <- lmBF(DV ~ A + B + ID + A:ID + B:ID, df, whichRandom="ID", iterations=1e5, progress=F)

set.seed(277)
omitA <- lmBF(DV ~ B + ID + A:ID + B:ID, df, whichRandom="ID", iterations=1e5, progress=F)

set.seed(277)
omitB <- lmBF(DV ~ A + ID + A:ID + B:ID, df, whichRandom="ID", iterations=1e5, progress=F)

# testing AB
add / full # 1.009295 ±6.56% BF₀₁

# testing A
omitA / add # 1.468665 ±3.45% BF₀₁

# testing B
omitB / add # 0.0001408403 ±6.04% BF₀₁


setwd("file path")
library(car)
library(vegan)

#----------------------------------------------------------------------
############# PERMANOVA for symbiont community comparisons #########################
# need metadata and need straight count data of symbiont communities

#SymmComm
SymData = read.csv("ADNA.csv")
SymData
SymDataCounts = SymData[,4:7]
SymDataCounts
adonis2(SymDataCounts~SymData$Species*SymData$Treatment, data = SymDataCounts)

#----------------------------------------------------------------------
#############Student and Welch's T-tests #########################
# T-tests were completed in a pariwise fashion within species by treatments
#T-tests were used to compare values for respirometry, PAM, CD, Chl, AFDW, Total lipid, and lipid classes

#load datasheet of interest
Data=read.csv("read in data sheet of interest.csv")
Data
#Subset datasheet for just pairwise comparison of interest
Pair = subset(Mote, Species != "Orbicella faveolata")
Pair = subset (Pair, Treatment != "Nursery2022")

#test assumptions for test
hist(Pair$Measure)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.9898 ->normal
shapiro.test(Pair$PR) 
#test for equal variance, p > 0.05 = equal variance
leveneTest(Pair$PR ~ Treatment, data = Pair)#p > 0.05 = equal variance

#if don't pass assumptions try transformations below and re-test assumptions
Pairlog=log(Pair$PR) #log transformation
PairRecip=(1/Pair$PR) #reciprocal transformation
PairSqrt=sqrt(Pair$FvFm+0.5) #Square root transformation

#if transformations don't work use kruskal-wallis test below, followed by Dunn if significant
#Kruskal- Wallis test how to https://www.statology.org/kruskal-wallis-test/
#checked distribution shapes by species in excel, both left skewed - same shape 
kruskal.test(Pair$FvFm~Pair$Species)
kruskal.test(Pair$FvFm~Pair$Treatment)
#Dunn test, non-parametric Tukey for pairwise
library("FSA")
dunnTest(FvFm~Treatment, data=Pair, method="bonferroni")
#if pass assumptions run paired t-test
t.test(Pair$PR ~ Treatment, data = Pair, var.equal = TRUE)#t = 1.8107, df = 6.3817, p-value = 0.1172




#originally I ran this code with separate sheets of the treatment and species of interest, 
#with values to test as columns as I was generating data. These original tests are shown below, 
#with any attempts at transformation used for that comparison/dataset. The above code for 
#t-tests was designed for use with the now complete dataset.

##Respirometry measurements t-tests
#PR APAM-N21
CDMote=read.csv("MoteOutplantDataAPN21-AM.csv")
CDMote
hist(CDMote$PR)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.9898 ->normal
shapiro.test(CDMote$PR) 
leveneTest(CDMote$PR ~ Treatment, data = CDMote)#p=0.333 -> equal variance
#paired t-test
t.test(CDMote$PR ~ Treatment, data = CDMote, var.equal = TRUE)#t = 1.8107, df = 6.3817, p-value = 0.1172

#PR APAM-N22
CDMote=read.csv("MoteOutplantDataAPAM-N22.csv")
CDMote
hist(CDMote$PR)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.9898 ->normal
shapiro.test(CDMote$PR) 
CDMotelog=log(CDMote$PR)
hist(CDMotelog)
shapiro.test(CDMotelog)
leveneTest(CDMotelog ~ Treatment, data = CDMote)#p=0.333 -> equal variance
#paired t-test
t.test(CDMotelog ~ Treatment, data = CDMote, var.equal = TRUE)#t = 1.8107, df = 6.3817, p-value = 0.1172

#PR APN21-APN22
CDMote=read.csv("MoteOutplantDataAPN21-N22.csv")
CDMote
hist(CDMote$PR)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.00131 ->not normal
shapiro.test(CDMote$PR) 
CDMotelog=log(CDMote$PR)
hist(CDMotelog)
shapiro.test(CDMotelog)#p-value = 0.01696 -> not normal, try reciprocal
CDMoteRecip=(1/CDMote$PR)
hist(CDMoteRecip)
shapiro.test(CDMoteRecip)#p-value = 0.1554 -> normal
leveneTest(CDMoteRecip ~ Treatment, data = CDMote)#p=0.3828 -> equal variance
#paired t-test
t.test(CDMoteRecip ~ Treatment, data = CDMote, var.equal = TRUE)#t = -0.40348, df = 6.7014, p-value = 0.6992

#PR OFN-CI
CDMote=read.csv("MoteOutplantDataOFN-CI.csv")
CDMote
hist(CDMote$PR)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.002572 -> not normal
shapiro.test(CDMote$PR) 
CDMotelog=log(CDMote$PR)
hist(CDMotelog)
shapiro.test(CDMotelog)#p-value =  0.01521 -> not normal, try reciprocal
CDMoteRecip=(1/CDMote$PR)
hist(CDMoteRecip)
shapiro.test(CDMoteRecip)#p-value = 0.04675 -> not normal, but rounds to normal
leveneTest(CDMoteRecip ~ Treatment, data = CDMote)#p = 0.09521 -> equal variance
#paired t-test
t.test(CDMoteRecip ~ Treatment, data = CDMote, var.equal = TRUE)#t = 4.048, df = 3.4709, p-value = 0.02052

#Resp APN21-AM
CDMote=read.csv("MoteOutplantDataAPN21-AM.csv")
CDMote
hist(CDMote$Resp)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.1369 ->normal
shapiro.test(CDMote$Resp) 
leveneTest(CDMote$Resp ~ Treatment, data = CDMote)#p=0.0593 -> equal variance
#paired t-test
t.test(CDMote$Resp ~ Treatment, data = CDMote, var.equal = TRUE)#t = 2.7092, df = 6.6471, p-value = 0.03177

#Resp APN22-AM
CDMote=read.csv("MoteOutplantDataAPAM-N22.csv")
CDMote
hist(CDMote$Resp)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.1369 ->normal
shapiro.test(CDMote$Resp) 
leveneTest(CDMote$Resp ~ Treatment, data = CDMote)#p=0.0593 -> equal variance
#paired t-test
t.test(CDMote$Resp ~ Treatment, data = CDMote, var.equal = TRUE)#t = 2.7092, df = 6.6471, p-value = 0.03177

#Resp APN21-N22
CDMote=read.csv("MoteOutplantDataAPN21-N22.csv")
CDMote
hist(CDMote$Resp)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.4691 ->normal
shapiro.test(CDMote$Resp) 
leveneTest(CDMote$Resp ~ Treatment, data = CDMote)#p=0.1598 -> equal variance
#paired t-test
t.test(CDMote$Resp ~ Treatment, data = CDMote, var.equal = TRUE)#t = 1.812, df = 8.3174, p-value = 0.1061

#Resp OFN-CI
CDMote=read.csv("MoteOutplantDataOFN-CI.csv")
CDMote
hist(CDMote$Resp)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.9583 ->normal
shapiro.test(CDMote$Resp) 
leveneTest(CDMote$Resp ~ Treatment, data = CDMote)#p=0.8398 -> equal variance
#paired t-test
t.test(CDMote$Resp ~ Treatment, data = CDMote, var.equal = TRUE)#t = 1.58, df = 6.4947, p-value = 0.1614

#Pgross APN21-Am
CDMote=read.csv("MoteOutplantDataAPN21-AM.csv")
CDMote
hist(CDMote$Pgross)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.04326 ->not normal
shapiro.test(CDMote$Pgross) 
CDMotelog=log(CDMote$Pgross)
hist(CDMotelog)
shapiro.test(CDMotelog)#p-value = 0.3606 -> normal
leveneTest(CDMotelog ~ Treatment, data = CDMote)#p=0.7015 -> equal variance
#paired t-test
t.test(CDMotelog ~ Treatment, data = CDMote, var.equal = TRUE)#t = 4.3764, df = 9.7493, p-value = 0.001471

#Pgross APN22-Am
CDMote=read.csv("MoteOutplantDataAPAM-N22.csv")
CDMote
hist(CDMote$Pgross)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.04326 ->not normal
shapiro.test(CDMote$Pgross) 
leveneTest(CDMote$Pgross ~ Treatment, data = CDMote)#p=0.7015 -> equal variance
#paired t-test
t.test(CDMotelog ~ Treatment, data = CDMote, var.equal = TRUE)#t = 4.3764, df = 9.7493, p-value = 0.001471

#Pgross APN21-N22
CDMote=read.csv("MoteOutplantDataAPN21-N22.csv")
CDMote
hist(CDMote$Pgross)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.03138 ->not normal
shapiro.test(CDMote$Pgross) 
CDMotelog=log(CDMote$Pgross)
hist(CDMotelog)
shapiro.test(CDMotelog)#p-value = 0.2265 -> normal
leveneTest(CDMotelog ~ Treatment, data = CDMote)#p=0.09405 -> equal variance
#paired t-test
t.test(CDMotelog ~ Treatment, data = CDMote, var.equal = TRUE)#t = 2.4727, df = 7.774, p-value = 0.03938

#Pgross OFN-CI
CDMote=read.csv("MoteOutplantDataOFN-CI.csv")
CDMote
hist(CDMote$Pgross)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.63 -> normal
shapiro.test(CDMote$Pgross) 
leveneTest(CDMote$Pgross ~ Treatment, data = CDMote)#p=0.7391 -> equal variance
#paired t-test
t.test(CDMote$Pgross ~ Treatment, data = CDMote, var.equal = TRUE)#t = -0.72243, df = 7.7048, p-value = 0.4914

#----------------------------------------------------------------------
#CD and Chl measurements t-tests
CDMote=read.csv("MoteOutplantData.csv")
CDMote
#######CD
#CD APAM-N22
CDMote=read.csv("MoteOutplantDataAPAM-N22.csv")
CDMote
hist(CDMote$pgchlacell)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.1689 ->normal
shapiro.test(CDMote$pgchlacell) 
leveneTest(CDMote$pgchlacell ~ Treatment, data = CDMote)#p=0.4971 -> equal variance
#paired t-test
t.test(CDMote$pgchlacell ~ Treatment, data = CDMote, var.equal = TRUE)#t = 2.1276, df = 8.1769, p-value = 0.0653

#CD APAM-N21
CDMote=read.csv("MoteOutplantDataAPN21-AM.csv")
CDMote
hist(CDMote$CD)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.1689 ->normal
shapiro.test(CDMote$CD) 
leveneTest(CDMote$CD ~ Treatment, data = CDMote)#p=0.4971 -> equal variance
#paired t-test
t.test(CDMote$CD ~ Treatment, data = CDMote, var.equal = TRUE)#t = 2.1276, df = 8.1769, p-value = 0.0653

#CD AP N21-N22
CDMote=read.csv("MoteOutplantDataAPN21-N22.csv")
CDMote
hist(CDMote$CD)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.1012 ->normal
shapiro.test(CDMote$CD) 
leveneTest(CDMote$CD ~ Treatment, data = CDMote)#p=0.5268 -> equal variance
#paired t-test
t.test(CDMote$CD ~ Treatment, data = CDMote, var.equal = TRUE)#t = -5.6569, df = 9.0033, p-value = 0.0003104

#CD OF N21-CI
CDMote=read.csv("MoteOutplantDataOFN-CI.csv")
CDMote
hist(CDMote$CD)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.9663 ->normal
shapiro.test(CDMote$CD) 
leveneTest(CDMote$CD ~ Treatment, data = CDMote)#p = 0.7818 -> equal variance
#paired t-test
t.test(CDMote$CD ~ Treatment, data = CDMote, var.equal = TRUE)#t = 1.7397, df = 8.0678, p-value = 0.1198

######Chl
#chlap AP N21-AM
CDMote=read.csv("MoteOutplantDataAPN21-AM.csv")
CDMote
hist(CDMote$pgchlacell)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.003069 -> not normal
shapiro.test(CDMote$pgchlacell) 
CDMotelog = log(CDMote$pgchlacell)
hist(CDMotelog)
shapiro.test(CDMotelog) #p-value = 0.3071 -> normal
leveneTest(CDMotelog ~ Treatment, data = CDMote)#p = 0.05351 -> equal variance
#paired t-test
t.test(CDMotelog ~ Treatment, data = CDMote, var.equal = TRUE)#t = 3.339, df = 6.6295, p-value = 0.01347

#chlap AP N21-N22
CDMote=read.csv("MoteOutplantDataAPN21-N22.csv")
CDMote
hist(CDMote$pgchlacell)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.003352 -> not normal
shapiro.test(CDMote$pgchlacell) 
CDMotelog = log(CDMote$pgchlacell)
hist(CDMotelog)
shapiro.test(CDMotelog) #p-value = 0.313 -> normal
leveneTest(CDMotelog ~ Treatment, data = CDMote)#p = 0.02192 * -> not equal variance
#paired t-test
t.test(CDMotelog ~ Treatment, data = CDMote)#t = 1.3205, df = 5.839, p-value = 0.2361

#chlap OF N21-CI
CDMote=read.csv("MoteOutplantDataOFN-CI.csv")
CDMote
hist(CDMote$pgchlacell)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.7593 -> normal
shapiro.test(CDMote$ugchlacm2) 
leveneTest(CDMote$ugchlacm2 ~ Treatment, data = CDMote)#p = 0.5057 -> equal variance
#paired t-test
t.test(ugchlacm2 ~ Treatment, data = CDMote, var.equal = TRUE)#t = 3.7589, df = 6.9773, p-value = 0.007128

#chlau AP N21-AM
CDMote=read.csv("MoteOutplantDataAPN21-AM.csv")
CDMote
hist(CDMote$ugchlacm2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.009044 -> normal
shapiro.test(CDMote$ugchlacm2) 
CDMotelog = log(CDMote$ugchlacm2)
hist(CDMotelog)
shapiro.test(CDMotelog)#p-value = 0.4899 -> normal
leveneTest(CDMotelog ~ Treatment, data = CDMote)#p-value = 0.4899 -> equal variance
#paired t-test
t.test(CDMotelog ~ Treatment, data = CDMote, var.equal = TRUE)#t = 3.4819, df = 9.6898, p-value = 0.006187

#chlau AP N21-N22
CDMote=read.csv("MoteOutplantDataAPN21-N22.csv")
CDMote
hist(CDMote$ugchlacm2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.06229 -> normal
shapiro.test(CDMote$ugchlacm2) 
leveneTest(CDMote$ugchlacm2 ~ Treatment, data = CDMote)#p = 0.01607 * -> not equal variance
#paired t-test
t.test(CDMote$ugchlacm2 ~ Treatment, data = CDMote)#t = 0.17013, df = 5.1369, p-value = 0.8714

#chlau OF N21-AM
CDMote=read.csv("MoteOutplantDataOFN-CI.csv")
CDMote
hist(CDMote$ugchlacm2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.5576 -> normal
shapiro.test(CDMote$ugchlacm2) 
leveneTest(CDMote$ugchlacm2 ~ Treatment, data = CDMote)#p = 0.2075 -> not equal variance
#paired t-test
t.test(CDMote$ugchlacm2 ~ Treatment, data = CDMote, var.equal = TRUE)#t = 4.967, df = 6.5041, p-value = 0.002011

#----------------------------------------------------------------------
#AFDW t-tests
#AP AM-N21 AFDW
AFDWMote=read.csv("MoteOutplantRsheetsAFDWAPN21-AM.csv")
AFDWMote
#######AFDW
hist(AFDWMote$AFDWmgcm.2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.753 ->  normal 
shapiro.test(AFDWMote$AFDWmgcm.2)
leveneTest(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote)#p = 0.8502 -> not equal variance
#paired t-test
t.test(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote, var.equal = TRUE)#t = -1.5279, df = 9.1702, p-value = 0.1603

#AP AM-N22 AFDW
AFDWMote=read.csv("MoteOutplantRsheetsAFDWAPN22-AM.csv")
AFDWMote
#######AFDW
hist(AFDWMote$AFDWmgcm.2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.175 -> normal
shapiro.test(AFDWMote$AFDWmgcm.2) 
leveneTest(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote)#p = 0.7401 -> not equal variance
#paired t-test
t.test(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote, var.equal = TRUE)#t = 1.0195, df = 8.9934, p-value = 0.3346

#AP N21-N22 AFDW
AFDWMote=read.csv("MoteOutplantRsheetsAFDWAPN21-N22.csv")
AFDWMote
#######AFDW
hist(AFDWMote$AFDWmgcm.2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.7703 -> normal
shapiro.test(AFDWMote$AFDWmgcm.2) 
leveneTest(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote)#p = 0.8308 -> not equal variance
#paired t-test
t.test(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote, var.equal = TRUE)#t = -3.0873, df = 9.9859, p-value = 0.01151

#OF AM-N21 AFDW
AFDWMote=read.csv("MoteOutplantRsheetsAFDWOFN21-AM.csv")
AFDWMote
#######AFDW
hist(AFDWMote$AFDWmgcm.2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.1532 -> normal 
shapiro.test(AFDWMote$AFDWmgcm.2) 
leveneTest(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote)#p = 0.2335 -> not equal variance
#paired t-test
t.test(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote)#t = 2.6285, df = 5.1283, p-value = 0.04548

#OF CI-N21 AFDW
AFDWMote=read.csv("MoteOutplantRsheetsAFDWOFN21-CI.csv")
AFDWMote
#######AFDW
hist(AFDWMote$AFDWmgcm.2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.1583 -> normal
shapiro.test(AFDWMote$AFDWmgcm.2) 
leveneTest(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote)#p = 0.2707 -> not equal variance
#paired t-test
t.test(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote, var.equal = TRUE)#t = 1.6396, df = 5.9172, p-value = 0.1529

#OF CI-AM AFDW
AFDWMote=read.csv("MoteOutplantRsheetsAFDWOFAM-CI.csv")
AFDWMote
#######AFDW
hist(AFDWMote$AFDWmgcm.2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3892 -> normal
shapiro.test(AFDWMote$AFDWmgcm.2) 
leveneTest(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote)#p = 0.2262 -> not equal variance
#paired t-test
t.test(AFDWMote$AFDWmgcm.2 ~ Treatment, data = AFDWMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#----------------------------------------------------------------------
##weight and surface area t-tests
buoyMote=read.csv("MoteOutplantDataBuoyAP.csv")
buoyMote
#Avg AP only area
hist(buoyMote$Change.Total.Area.cm2)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.2094 -> normal
shapiro.test(buoyMote$Change.Total.Area.cm2) 
leveneTest(buoyMote$Change.Total.Area.cm2 ~ Treatment, data=buoyMote) #p=0.8494 -> variance equal
#paired t-test
t.test(buoyMote$Change.Total.Area.cm2 ~ Treatment, data=buoyMote, var.equal = TRUE) #t = 2.5556, df = 5.9104, p-value = 0.04374

#AP only Avg  change weight per change area per day
hist(buoyMote$change.mg.cm.2.day.1)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.2203 -> normal
shapiro.test(buoyMote$change.mg.cm.2.day.1) 
leveneTest(buoyMote$change.mg.cm.2.day.1 ~ Treatment, data=buoyMote) #p=0.6911 -> variance equal
#paired t-test
t.test(buoyMote$change.mg.cm.2.day.1 ~ Treatment, data=buoyMote, var.equal = TRUE) #t = -5.5974, df = 9, p-value = 0.0003354

#----------------------------------------------------------------------
FvFmMote=read.csv("MoteOutplantData.csv")
FvFmMote
#######FvFm measurement t-tests
hist(FvFmMote$FvFm)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 4.583e-05 -> not normal try log 
shapiro.test(FvFmMote$FvFm) 
logFvFmMote=log(1-FvFmMote$FvFm)
hist(logFvFmMote)
shapiro.test(logFvFmMote) #p-value = 0.001229 -> not normal try sqrt
sqrtFvFmMote=sqrt(FvFmMote$FvFm+0.5)
hist(sqrtFvFmMote)
shapiro.test(sqrtFvFmMote)#p-value = 2.109e-05 -> not normal try square
sqrFvFmMote=(FvFmMote$FvFm)^2
hist(sqrFvFmMote)
shapiro.test(sqrFvFmMote)#p-value = 0.0008728 -> not normal try antilog
antiFvFmMote=exp(FvFmMote$FvFm)
hist(antiFvFmMote)
shapiro.test(antiFvFmMote)# p-value = 0.0002163 -> not normal try reciprocal
recipFvFmMote=(1/FvFmMote$FvFm)
hist(recipFvFmMote)
shapiro.test(recipFvFmMote)#p-value = 9.755e-08 -> not normal, transformations are unsuccessful, try non-parametric Krukskal-Wallis
#Kruskal- Wallis test how to https://www.statology.org/kruskal-wallis-test/
#checked distribution shapes by species in excel, both left skewed - same shape 
kruskal.test(FvFmMote$FvFm~FvFmMote$Species)
kruskal.test(FvFmMote$FvFm~FvFmMote$Treatment)
#Dunn test, non-parametric Tukey for pairwise
install.packages("FSA")
library("FSA")
dunnTest(FvFm~Sp_Treat, data=FvFmMote, method="bonferroni")
#______________________________________________________________________
##Lipid class t-tests

#OF  CI-AM - switched through lipid class variables as needed
WaxMote=read.csv("MotePOR_LCforR_OF_CI-Am.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. 
shapiro.test(WaxMote$Phophs) 
logST=log(WaxMote$ST)
logAMPL=log(WaxMote$AMPL)
logTotalLipid=log(WaxMote$TotalLipid) #- not normal
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.04132 -> not equal variance
#paired t-test
t.test(WaxMote$ST ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#OF  Oct21-CI
WaxMote=read.csv("MotePOR_LCforR_OF_Oct21-CI.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$Phophs) 
logTAG=log(WaxMote$TAG) #- not normal or equal variance
logAMPL=log(WaxMote$AMPL) #- not normal or equal variance
logTotalLipid=log(WaxMote$TotalLipid)
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#OF  Oct21-AM
WaxMote=read.csv("MotePOR_LCforR_OF_Oct21-Am.csv")
WaxMote
hist(WaxMote$ST)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$ST) 
logAMPL=log(WaxMote$AMPL)
logPhophs=log(WaxMote$Phophs)
logTotalLipid=log(WaxMote$TotalLipid)
leveneTest(WaxMote$ST ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$ST ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#OF  Oct21-Mar21
WaxMote=read.csv("MotePOR_LCforR_OF_Oct21-Mar21.csv")
WaxMote
#######AFDW
hist(WaxMote$AMPL)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(logPhophs) 
logTAG=log(WaxMote$TAG)
logPhophs=log(WaxMote$Phophs)
leveneTest(logPhophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#t-test
t.test(logPhophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#OF  Mar21-AM
WaxMote=read.csv("MotePOR_LCforR_OF_Mar21-Am.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$Phophs) 
logAMPL=log(WaxMote$AMPL)
logTotalLipid=log(WaxMote$TotalLipid)
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#OF  Mar21-CI
WaxMote=read.csv("MotePOR_LCforR_OF_Mar21-CI.csv")
WaxMote
#######AFDW
hist(WaxMote$AMPL)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(logAMPL) 
logTAG=log(WaxMote$TAG)
logAMPL=log(WaxMote$AMPL)
logTotalLipid=log(WaxMote$TotalLipid)
leveneTest(logAMPL ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(logAMPL ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#AP Oct21-Oct22 - switched through lipid class variables as needed
WaxMote=read.csv("MotePOR_LCforR_AP_Oct21-Oct22.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$Phophs) 
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#AP Oct21-Am
WaxMote=read.csv("MotePOR_LCforR_AP_Oct21-Am.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$Phophs) 
logWAX=log(WaxMote$WAX)
logTAG=log(WaxMote$TAG)
logAMPL=log(WaxMote$AMPL)
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#AP Oct22-AM
WaxMote=read.csv("MotePOR_LCforR_AP_Oct22-Am.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$Phophs) 
logWax=log(WaxMote$WAX)
logST=log(WaxMote$ST)
logAMPL=log(WaxMote$AMPL)
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#AP Mar21-Oct21
WaxMote=read.csv("MotePOR_LCforR_AP_Mar21-Oct21.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$Phophs) 
logWAX=log(WaxMote$WAX)
logTotalLipid=log(WaxMote$TotalLipid)
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#AP Mar21-Oct22
WaxMote=read.csv("MotePOR_LCforR_AP_Mar21-Oct22.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.3478 -> normal
shapiro.test(WaxMote$Phophs) 
logWAX=log(WaxMote$WAX)
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#AP Mar21-Am
WaxMote=read.csv("MotePOR_LCforR_AP_Mar21-Am.csv")
WaxMote
#######AFDW
hist(WaxMote$Phophs)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist.
shapiro.test(WaxMote$Phophs) 
logWAX=log(WaxMote$WAX)
logTAG=log(WaxMote$TAG)
logAMPL=log(WaxMote$AMPL)
leveneTest(WaxMote$Phophs ~ Treatment, data = WaxMote)#pvalue = 0.07081 -> equal variance
#paired t-test
t.test(WaxMote$Phophs ~ Treatment, data = WaxMote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181

#Protein AP Oct22 Nur-Am
Mote=read.csv("MotePOR_LCforR_AP_Oct22-Am.csv")
Mote
#######AFDW
hist(Mote$Protein)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.02738 -> not normal try log
shapiro.test(Mote$Protein) 
logProtein=log(Mote$Protein)
shapiro.test(logProtein) #p-value = 0.5124 -> normal
leveneTest(logProtein ~ Treatment, data = Mote)#pvalue = 0.241 -> equal variance
#paired t-test
t.test(logProtein ~ Treatment, data = Mote, var.equal= TRUE)#t = 0.11371, df = 10, p-value = 0.9117
  
#Carb AP Oct22 Nur-Am
Mote=read.csv("MotePOR_LCforR_AP_Oct22-Am.csv")
Mote
#######AFDW
hist(Mote$Carb)
#test for normality with Shapiro-Wilks test to have official consideration of normal
#p<0.05 = not a normal dist. p-value = 0.6451 -> normal
shapiro.test(Mote$Carb) 
leveneTest(Mote$Carb ~ Treatment, data = Mote)#pvalue = 0.5733 -> equal variance
#paired t-test
t.test(Mote$Carb ~ Treatment, data = Mote, var.equal = TRUE)#t = -2.811, df = 3.7604, p-value = 0.05181





#----------------------------------------------------------------------
############# Isotopes Stats #########################

#removed american shoals OF from dataset before running analyses, due to lack of replication sym n=1 host n=2
Mote=read.csv("IsotopeData_MotePOR_Outplant_NoSkel.csv") #Two-way Anova (uneven N)
Mote=read.csv("SymToHostN_IsotopeData_MotePOR_Outplant.csv")
#Mote = read.csv("IsotopeData_MotePOR_OutplantCN.csv") #C:N ratio , none normal no matter how transformed
#Mote=read.csv("IsotopeData_MotePOR_Outplant_NoSkel_repeatedMeasAnova.csv") #repeated measures ANOVA (equal N)
Mote

MoteHost = subset(Mote, Type == "Host")
MoteHAP = subset(MoteHost, Species == "Acroporapalmata")
MoteHOF = subset(MoteHost, Species == "Orbicellafaveolata")

MoteSym = subset(Mote, Type == "Symbiodiniaceae")
MoteSAP = subset(MoteSym, Species == "Acroporapalmata")
MoteSOF = subset(MoteSym, Species == "Orbicellafaveolata")
#-------------Kruskal and Dunn--------
#CNratios

#Kruskal- Wallis test how to https://www.statology.org/kruskal-wallis-test/
hist(MoteHAP$CNratio)
#checked distribution shapes by species in excel, both left skewed - same shape 
kruskal.test(MoteHOF$CNratio~MoteHOF$Treatment)
kruskal.test(MoteHAP$CNratio~MoteHAP$Treatment)


#Dunn test, non-parametric Tukey for pairwise
library("FSA")
dunnTest(CNratio~Treatment, data=MoteSym, method="bonferroni")


#----Requirements for Anova------
#C, N, Chost-Csym, Nhost -Csym, skeleton C
#does it look normal? can use histogram and/or shapiro test
hist(transHAP)
sqrtSymN = sqrt(MoteSym$δ15NAIR)
sqrtHostN = sqrt(MoteHost$δ15NAIR)
logskel = log(abs(Mote$SkeldC))
transHAP = (MoteHAP$CNratio) # not tranformable to normal - use Kruskal-Wallis

#test for normality - may not be good option for this dataset
shapiro.test(transHAP) # normal if > 0.05

#does it have equal variance?
leveneTest(MoteSAP$CNratio ~ Treatment, data = MoteSAP) #equal variance if > 0.05

#two-way ANOVA
ANOVAmodel <- lm(MoteSym$CNratio ~ Treatment*Species, data = MoteSym)
anova(ANOVAmodel)

#one-way ANOVA
ANOVAmodel <- aov(MoteSAP$CNratio ~ Treatment, data = MoteSAP)
anova(ANOVAmodel)


#Tukey test for significant comparisons
ANOVAmodel <- lm(logskel ~ TreatmentSkel*SpeciesSkel, data=Mote)
SBPaov <- aov(ANOVAmodel)
TukeyHSD(SBPaov)


#----------------------------------------------------------------------
# Graphing Data - exported as 6 x 8 PDFs
#Good website for checking colors to account for colorblindness
#https://davidmathlogic.com/colorblind/#%231572A1-%23b40c17
install.packages("ggplot2")
library(ggplot2)
library(tidyverse)
install.packages("plotrix")
library("plotrix")  
?ggplot
?geom_bar
?geom_errorbar
?geom_point

#PR
Data = read.csv("MoteOutplantData.csv")
Data
p <- ggplot(Data, aes(Species,PR,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Average P:R per cm^2")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks =8, expand = c(0, 0), limits = c(0,8))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#O2 Resp
Data = read.csv("MoteOutplantData.csv")
Data
p <- ggplot(Data, aes(Species,Resp,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Average Absolute O2 Respiration per cm^2")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 4, expand = c(0, 0), limits = c(0,2.0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Pnet
Data = read.csv("MoteOutplantData.csv")
Data
p <- ggplot(Data, aes(Species,Pnet,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Average Net Photosynthesis per cm^2")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 4, expand = c(0, 0), limits = c(0,2.0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Pgross
Data = read.csv("MoteOutplantData.csv")
Data
p <- ggplot(Data, aes(Species,Pgross,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Average Gross Photosynthesis per cm^2")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 4, expand = c(0, 0), limits = c(0,2.0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#CD
CDData = read.csv("MoteOutplantData.csv")
CDData
p <- ggplot(CDData, aes(Species,CD,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Cell cm-2")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 4, expand = c(0, 0), limits = c(0, 2500000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#FvFm
FvFmData = read.csv("MoteOutplantData.csv")
FvFmData
p <- ggplot(FvFmData, aes(Species,FvFm,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Average Fv/Fm")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 1.0)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 1.0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8))

#AFDW
AFDWData = read.csv("MoteOutplantRsheetsAFDW.csv")
AFDWData
p <- ggplot(AFDWData, aes(Species,AFDWmgcm2,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="AFDW mg cm-2")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e"))+
  #scale_y_continuous(n.breaks = 8)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 7)) +
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 12))+
  #expand_limits(y = c(0, 12)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8))
  
#Chlap
ChlapData = read.csv("MoteOutplantData.csv")
ChlapData
p <- ggplot(ChlapData, aes(Species,pgchlacell,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="pg chl a cell-1")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e"))+
  #scale_y_continuous(n.breaks = 8)+
  #scale_x_continuous(expand = c(0, 0), limits = c(0, 7)) +
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 12))+
  #expand_limits(y = c(0, 12)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8))
  
#Chlau
ChlauData = read.csv("MoteOutplantData.csv")
ChlauData
p <- ggplot(ChlauData, aes(Species,ugchlacm2,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+#adjust text sizes
  labs(x="",y="ug chl a cm-2")+
  scale_fill_manual(values=c("#FAC34C", "#A4A19B", "#056273", "#78924e"))+
  #scale_y_continuous(n.breaks = 8)+#number y axis breaks
  #expand_limits(y = c(0, 30)) + #y axis limits
  scale_y_continuous(n.breaks = 5, expand = c(0, 0), limits = c(0, 20))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8))

#Buoyant Weight Data
#change area
Data = read.csv("MoteOutplantDataBuoy.csv")
Data
p <- ggplot(Data, aes(Species,Change.Total.Area.cm2,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Average change cm2")+
  scale_fill_manual(values=c("#FAC34C", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 12, expand = c(0, 0), limits = c(-5, 60))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#change weight per change area per day
Data = read.csv("MoteOutplantDataBuoy.csv")
Data
p <- ggplot(Data, aes(Species,change.mg.cm.2.day.1,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="g wax per g AFDW")+
  scale_fill_manual(values=c("#FAC34C", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(-0.5, 2.0))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Energy Reserves and Lipids
#Total Lipids 
Data = read.csv("MotePOR_LCforR.csv")
Data
p <- ggplot(Data, aes(Species,J.per.g.AFDW,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="J g-1 AFDW")+
  scale_fill_manual(values=c("#ff551f", "#FAC34C","#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 70000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Lipid classes - WAX
Data = read.csv("MotePOR_LCforR.csv")
Data
p <- ggplot(Data, aes(Species,WAX,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="J Wax AFDW g-1")+
  scale_fill_manual(values=c("#FAC34C", "#FAC34C","#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 500))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Lipid classes - TAG
Data = read.csv("MotePOR_LCforR.csv")
Data
p <- ggplot(Data, aes(Species,TAG,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="TAG per AFDW g")+
  scale_fill_manual(values=c("#FAC34C", "#FAC34C","#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 40))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Lipid classes - ST
Data = read.csv("MotePOR_LCforR.csv")
Data
p <- ggplot(Data, aes(Species,ST,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="ST per AFDW g")+
  scale_fill_manual(values=c("#FAC34C", "#FAC34C","#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 40))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Lipid classes - AMPL
Data = read.csv("MotePOR_LCforR.csv")
Data
p <- ggplot(Data, aes(Species,AMPL,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="AMPL per AFDW g")+
  scale_fill_manual(values=c("#FAC34C", "#FAC34C","#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 160))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Lipid classes - Phosphoplipids
Data = read.csv("MotePOR_LCforR.csv")
Data
p <- ggplot(Data, aes(Species,Phophs,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Phophs per AFDW g")+
  scale_fill_manual(values=c("#FAC34C", "#FAC34C","#A4A19B", "#056273", "#78924e")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 120))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Protein
Data = read.csv("MotePOR_LCforR_AP_Oct22-Am.csv")
Data
p <- ggplot(Data, aes(Species,Protein,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Protein J per AFDW g")+
  scale_fill_manual(values=c("#A4A19B", "#056273")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0, 3000))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 

#Carbohydrates
Data = read.csv("MotePOR_LCforR_AP_Oct22-Am.csv")
Data
p <- ggplot(Data, aes(Species,Carb,fill = Treatment))
p + geom_bar(position = 'dodge', stat = 'summary') + #make bar
  theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
  labs(x="",y="Carbs J per AFDW g")+
  scale_fill_manual(values=c("#A4A19B", "#056273")) +
  #scale_y_continuous(n.breaks = 8)+
  #expand_limits(y = c(0, 8000000)) + 
  scale_y_continuous(n.breaks = 8, expand = c(0, 0), limits = c(0,0.25))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
  geom_errorbar(stat = 'summary', position = position_dodge(.9), width = 0.5, ) + #add standard error bars
  #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 


#Isotopes
Data = read.csv("IsotopeData_MotePOR_Outplant_NoSkel.csv")
OF <- subset(Data, Species != "Acroporapalmata")
AP <- subset(Data, Species != "Orbicellafaveolata")
Data 
  
#skeleton C plot
  Data = read.csv("SymToHostN_IsotopeData_MotePOR_Outplant.csv")
  Data
  p <- ggplot(Data, aes(SpeciesSkel,SkeldC,fill = TreatmentSkel))
  p + geom_boxplot(outlier.shape=1) + #make bar
    theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
    labs(x="",y="Skeleton d13C")+
    scale_fill_manual(values=c("#ff551f", "#FAC34C", "#056273", "#78924e")) +
    #scale_y_continuous(trans = "reverse")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
    #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
    geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 
  
#C:N plot
Data = read.csv("IsotopeData_MotePOR_OutplantCN.csv")
Data
OF = subset(Data, Species == "Orbicellafaveolata")
AP = subset(Data, Species == "Acroporapalmata")
p <- ggplot(Data, aes(SpType,CNratio,fill = Treatment))
p + geom_boxplot(outlier.shape=1) + #make boxplot
   theme(legend.text=element_text(size=15), axis.title=element_text(size=12), axis.text=element_text(size=12))+
   labs(x="",y="C:N Ratio")+
   scale_fill_manual(values=c("#ff551f", "#FAC34C", "#056273", "#78924e", "black")) +
   #scale_y_continuous(trans = "reverse")+
   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
         panel.background = element_blank(), axis.line = element_line(colour = "black"))+#remove grid background and make axis black
   #geom_point(aes(x = ID), shape = 21, position = position_dodge(width = 1))#add associated data points
   geom_jitter(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.8)) 
  
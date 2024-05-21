## SIBER Script following Jackson et al. 2011 
## Jackson, A.L., Parnell, A.C., Inger R., & Bearhop, S. 2011. Comparing isotopic niche 
## widths among and within communities: SIBER - Stable Isotope Bayesian Ellipses in R. 
## Journal of Animal Ecology, 80, 595-602. 
## DOI: https://doi.org/10.1111/j.1365-2656.2011.01806.x 
## Example vignette: https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html #####
## Updated by Jackson 28 April 2017 
## Downloaded 02 October 2018
## Used in Conti-Jerpe et al. 2020


#headings of columns in .csv must be "iso1", "iso2", "group", and "community"
#all headings & columns must be present

rm(list=ls())
graphics.off()
set.seed(1)

#had to instal JAG first from this source before could install rjags for use
#https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/
#install.packages("rjags")
#install.packages("SIBER")
library(rjags)
library(SIBER)

setwd("file path")

#practice dataset available see SIBER webpage for details
#https://cran.r-project.org/web/packages/SIBER/vignettes/Introduction-to-SIBER.html

#Full Mote dataset
#Mote <- read.csv("MoteIsotope_TrophicNiche_NoSkel.csv", header=T)

ACRO <- read.csv("MoteIsotope_TrophicNiche_AP_NoSkel_integers.csv", header=T)
FAV <- read.csv("MoteIsotope_TrophicNiche_OF_NoSkel_integers.csv", header=T)

siber.example <- createSiberObject(ACRO)
siber.example <- createSiberObject(FAV)


# calculate metrics
community.hulls.args <- list(col = 1, lty = 1, lwd = 1)
group.ellipses.args  <- list(n = 100, p.interval = 0.40, lty = 1, lwd = 2)
group.hull.args      <- list(lty = 2, col = "grey20")


#Plot data with convex hulls and SEA
palette(c("#ACA4E2","#5CBD92", "#000"))
points <- c(17, 17, 16, 16, 15)


par(mfrow=c(1,1))
plotSiberObject(siber.example,
                ax.pad = 2, 
                points.order = points,
                hulls = F, community.hulls.args, 
                ellipses = T, group.ellipses.args,
                group.hulls = T, group.hull.args,
                bty = "o",
                iso.order = c(2,1),
                xlab = expression({delta}^13*C~'\u2030'),
                ylab = expression({delta}^15*N~'\u2030'),
)
?plotSiberObject

# Calculate summary statistics for each group: TA, SEA and SEAc
group.ML <- groupMetricsML(siber.example)
print(group.ML)
?groupMetricsML

#AP All - results
# 2021Nursery.Symbiodiniaceae 2021Nursery.Host
# TA                      12.29765         1.484200
# SEA                     11.45738         1.064407
# SEAc                    14.32172         1.330508
# 2022Nursery.Symbiodiniaceae 2022Nursery.Host
# TA                      8.591300         4.148000
# SEA                     6.003279         3.761089
# SEAc                    7.504099         4.701361
# AmericanShoals.Symbiodiniaceae AmericanShoals.Host
# TA                         4.224900           0.4003000
# SEA                        3.335359           0.3113045
# SEAc                       4.169199           0.3891307

#OF All - Results
# 2021Nursery.Symbiodiniaceae 2021Nursery.Host
# TA                     10.753800         2.354900
# SEA                     8.050937         1.769295
# SEAc                   10.063671         2.211618
# Cook'sIsland.Symbiodiniaceae Cook'sIsland.Host
# TA                      0.7242500         0.1712000
# SEA                     0.6473551         0.1425284
# SEAc                    0.8631401         0.1900379


# Calculate the various Layman metrics on each of the communities.
community.ML <- communityMetricsML(siber.example) 
print(community.ML)

#All AP-results
#           2021Nursery 2022Nursery AmericanShoals
# dY_range    3.223333   0.1083333      0.4966667
# dX_range    4.096667   3.7883333      1.5550000
# TA          0.000000   0.0000000      0.0000000
# CD          2.606365   1.8949410      0.8161959
# MNND        5.212730   3.7898820      1.6323917
# SDNND       0.000000   0.0000000      0.0000000

#All OF - Results
#            2021Nursery Cook'sIsland
# dY_range    1.010000     0.850000
# dX_range    6.791667     2.010000
# TA          0.000000     0.000000
# CD          3.433178     1.091169
# MNND        6.866355     2.182338
# SDNND       0.000000     0.000000

#Fitting Bayesian models to the data
# options for running jags
parms <- list()
parms$n.iter <- 2 * 10^4   # number of iterations to run the model for
parms$n.burnin <- 1 * 10^3 # discard the first set of values
parms$n.thin <- 10     # thin the posterior by this many
parms$n.chains <- 2        # run this many chains

# define the priors
priors <- list()
priors$R <- 1 * diag(2)
priors$k <- 2
priors$tau.mu <- 1.0E-3

# fit the ellipses which uses an Inverse Wishart prior
# on the covariance matrix Sigma, and a vague normal prior on the 
# means. Fitting is via the JAGS method.
ellipses.posterior <- siberMVN(siber.example, parms, priors)


# The posterior estimates of the ellipses for each group can be used to
# calculate the SEA.B for each group.
SEA.B <- siberEllipses(ellipses.posterior)


siberDensityPlot(SEA.B, xticklabels = colnames(group.ML), 
                 xlab = c("Community | Group"),
                 ylab = expression("Standard Ellipse Area " ('\u2030' ^2) ),
                 bty = "L",
                 las = 1,
                 main = "SIBER ellipses on each group"
)

# Add red x's for the ML estimated SEA-c
points(1:ncol(SEA.B), group.ML[3,], col="red", pch = "x", lwd = 2)


# Calculate some credible intervals 
cr.p <- c(0.95, 0.99) # vector of quantiles


# call to hdrcde:hdr using lapply()
SEA.B.credibles <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$hdr},
  prob = cr.p)


# do similar to get the modes, taking care to pick up multimodal posterior
# distributions if present
SEA.B.modes <- lapply(
  as.data.frame(SEA.B), 
  function(x,...){tmp<-hdrcde::hdr(x)$mode},
  prob = cr.p, all.modes=T)





### Ellipses overlap

# The first ellipse is referenced using a character string representation where 
# in "x.y", "x" is the community, and "y" is the group within that community.
# So in this example: community 1, group 2
ellipse1 <- "1.1" 

# Ellipse two is similarly defined: community 1, group3
ellipse2 <- "1.2"

# Overlap metric (units are per mil^2)
# The overlap of the maximum likelihood fitted standard ellipses are 
# estimated using
sea.overlap <- maxLikOverlap(ellipse1, ellipse2, siber.example, 
                             p.interval = NULL, n = 100)
sea.overlap

#AP results
# area.1    area.2   overlap 
# 14.312111  1.329615  0.000000 

#OF results
# area.1    area.2   overlap 
# 14.312111  1.329615  0.000000 

# overlap as a proportion of group 1 area (host)
prop.sea.over.1 <- sea.overlap[3] / (sea.overlap[1])
prop.sea.over.1

#AP - Results
# overlap 
# 0 

#OF results
# overlap 
# 0 




#________ Stats for Comparing host and symbiont groups_____________

## Following Turner et al. 2010   
## Turner, T.F., Collyer, M.L., & Krabbenhoft, T.J. 2010. A general hypothesis-testing 
## framework for stable isotope ratios in ecological studies. Ecology, 91, 2227-2233.
## DOI: https://doi.org/10.1890/09-1454.1
## Source code and examples: http://www.esapubs.org/archive/ecol/E091/157/suppl-1.htm
## Downloaded 17 June 2016
## Used in Conti-Jerpe et al. 2020


rm(list=ls())
graphics.off()

setwd("/Users/Shelby/Desktop/Kemp Lab/Mote POR Project/Mote POR Outplant Oct2021Oct2022/R code and analyses")

#source('Turner.et.al.ecology.source.r')
source("https://www.esapubs.org/archive/ecol/E091/157/Turner.et.al.ecology.source.r")

#headings of columns in .csv must be "group", "x", and "y"
ex.data<-read.csv("Am_MoteIsotope_TrophicNiche2_AP_NoSkel.csv",header=T)

Y<-as.matrix(ex.data[,2:3])


# Designate groups
group<-as.factor(ex.data[,1])
gp<-length(levels(group)) # number of groups
n.comp<-(gp^2-gp)/2 # number of possible comparisons
rownames(Y)<-group 

lm.gp<-lm(Y~group,x=T,model=T) # for estimating group means
#This outputs the centroid of one group and distance in x and y to the centroid of the second group
res.gp<-resid(lm.gp) # residuals of groups from group means
yhat.gp<-predict(lm.gp) # predicted values

lm.gp.red<-lm(Y~1) # this is the model reduced by the group factor: only estimates an overall mean

#PLOTS
par(mai=c(1,1,0.5,0.5))

Trophic.1<-ex.data[ex.data[,1]=="1",2:3]
Trophic.2<-ex.data[ex.data[,1]=="2",2:3]
All.pts<-rbind(Trophic.1,Trophic.2)


plot(All.pts,type='n',asp=1,xlab=expression(delta^13~C),ylab=expression(delta^15~N),cex.lab=1.7,cex.axis=1.7)
points(Trophic.1, pch=21,bg="#ACA4E2")
points(Trophic.2, pch=21,bg="#5CBD92")
Trophic.1.c<-colMeans(Trophic.1)
Trophic.2.c<-colMeans(Trophic.2)

points(Trophic.1.c[1],Trophic.1.c[2], pch=21,bg="#ACA4E2",cex=1.7,lwd=2.5)
points(Trophic.2.c[1],Trophic.2.c[2], pch=21,bg="#5CBD92",cex=1.7,lwd=2.5)


# DISPERSION MEASURES
ex1<-ds.prep(res.gp,group) # see source file
ex1.ds<-disp.stat(ex1) # see source file
ex1.ds

#Host = 1, sym = 2
#AP 2021 Nursery
#           mdc       mnn       ecc
# Group.1 1.214307 0.8163159 0.9802118
# Group.2 2.475098 1.2704093 0.4408163

#AP 2022 Nursery
#              mdc      mnn       ecc
#  Group.1 1.311354 0.950323 0.7613323
#  Group.2 1.834437 1.707188 0.8556642

#AP AM
#           mdc       mnn       ecc
# Group.1 0.4486287 0.3597586 0.9175850
# Group.2 1.2600451 1.0985450 0.6074886

#OF 2021 Nursery
#             mdc       mnn       ecc
# Group.1 0.9774339 0.7444973 0.5337119
# Group.2 3.3171256 1.3168023 0.9641578

#OF CI
# Group.1 0.3199708 0.2237757 0.8811718
# Group.2 0.6096557 0.6373881 0.8748639



# GROUP MEANS
gp.m<-group.means(Y,group) # finds the group means for the raw data
gp.m

#host =1 sym = 2
#AP 2021 Nursery 
# x         y
# 1 -24.44500 5.0283333
# 2 -21.22167 0.9316667

#AP 2022 nursery
# x        y
# 1 -19.41000 5.980000
# 2 -19.30167 2.191667

#AP Am
# x     y
# 1 -18.37333 3.175
# 2 -17.87667 1.620

#OF 2021 Nursery
# x         y
# 1 -23.84333  3.998333
# 2 -22.83333 -2.793333

#OF CI
# x    y
# 1 -17.67 3.57
# 2 -18.52 1.56


# CONTRASTS
mean.dif<-as.vector(dist(as.matrix(gp.m)))
mean.dif
#AP 2021 Nursery 
# 5.21273

#AP 2022 Nursery
# 3.789882

#AP Am
#1.632392

#OF 2021 Nursery
#6.866355

#OF CI
#2.182338

disp.dif<-ds.dif(ex1.ds) # Dimension names not given by output

disp.dif 
disp.dif$mdc
disp.dif$mnn
disp.dif$ecc
#mdc = mean centroid distance
#mnn = mean nearest neighbor
#ecc = eccentricity

#AP 2021 Nursery host =2 sym = 1
# $mdc
# dis 
# 1.260791 
# 
# $mnn
# dis 
# 0.4540934 
# 
# $ecc
# dis 
# 0.5393956 

#AP 2022 Nursery
# $mdc
# dis 
# 0.5230823 
# 
# $mnn
# dis 
# 0.7568645 
# 
# $ecc
# dis 
# 0.09433198 

#AP Am
# $mdc
# dis 
# 0.8114164 
# 
# $mnn
# dis 
# 0.7387863 
# 
# $ecc
# dis 
# 0.3100964 

#OF 2021 Nursery
# $mdc
# dis 
# 2.339692 
# 
# $mnn
# dis 
# 0.572305 
# 
# $ecc
# dis 
# 0.4304459 

#OF CI
# $mdc
# dis 
# 0.2896849 
# 
# $mnn
# dis 
# 0.4136124 
# 
# $ecc
# dis 
# 0.006307913 





#### PERMUTATION Hotelling's T2 test PROCEDURE ####

# Describe a function that needs
# x = the full linear model
# x2 = the reduced linear model
# y = the raw data
# g = group list
# p = the number of permutations to use (keep in mind that
# the observed values count as 1 random permutation)

# SETUP
permute<-function(x,x2,y,g,p){
  p.table<-NULL
  
  yhat<-predict(x)
  res<-resid(x)
  line<-nrow(res) # this defines the number of permutable objects
  yhat.2<-predict(x2)
  res.2<-resid(x2)
  
  mean.dif<-as.vector(dist(as.matrix(group.means(Y,g))))
  
  
  for(i in 1:p){
    
    #For dispersion tests
    
    line.rand<-sample(line,replace=FALSE) # This creates a line from 1 to the number
    # of permutable objects, then randomizes the order
    res.temp<-cbind(line.rand,res) # attaches the random order to the ordered matrix
    z<-(order(line.rand)) # randomizes matrix as it reorders line
    res.temp2<-as.matrix(res.temp[z,])# removes randomized line
    res.p<-res.temp2[,-1] # Rows of residuals are now randomized
    y.rand<-yhat+res.p # random values created
    
    # The resampling procedure above is the same as randomizing original values
    # (Since the reduced model only contains the overall mean)
    # But using residuals makes it applicable to multi-factor models
    
    lm.r<-lm(y.rand~g,x=T,model=T) # new linear model
    r<-resid(lm.r)
    yhat<-predict(lm.r)
    ex1.r<-ds.prep(r,g) #thanks to the functions, prep takes one step
    ex.r.ds<-disp.stat(ex1.r) # thanks to the function, dispersion stats require one step
    disp.dif.r<-ds.dif(ex.r.ds) # thanks to function, contrasts require one step
    
    
    # For means tests
    
    # uses different linear model
    # the null hypothesis that means are equal
    # means that an intercept model (i.e., defines only the overall mean)
    # is as viable as a group means model.
    # Thus, residuals are calculated from the intercept model (see above)
    # and random means are created despite no mean differences define by the model
    # this creates random distributions of outcomes under the null hypothesis
    
    res.2.temp<-cbind(line.rand,res.2)
    z<-(order(line.rand))
    res.2.temp2<-as.matrix(res.2.temp[z,])
    res.2.p<-res.temp2[,-1] # Rows of residuals are now randomized
    y.rand.2<-yhat.2+res.2.p
    
    
    gm.r<-group.means(y.rand.2,g)
    md.r<-as.vector(dist(as.matrix(gm.r)))
    
    result<-c(i,md.r,disp.dif.r$mdc,disp.dif.r$mnn,disp.dif.r$ecc) # bind all results together
    p.table<-rbind(p.table,result) # add them to a table, row by row
    
  }
  
  head<-NULL # create a header
  line1<-as.vector(c(0,mean.dif,disp.dif$mdc,disp.dif$mnn,disp.dif$ecc))
  
  # The following is a bunch of code simply for generating column names in the output
  
  cn<-length(mean.dif) # cn = column name
  
  test.list<-NULL
  if (cn>1) for(i in 1:cn){
    l1<-rep(i,cn)
    l2<-array(1:cn)
    l12<-cbind(l1,l2)
    test.list<-rbind(test.list,l12)
  }
  
  test.list2<-NULL
  if (cn>1) for(j in 1:nrow(test.list)){
    t<-test.list[j,]
    if(t[2]>t[1]) test.list2<-rbind(test.list2,t)
  }
  
  test.list3<-NULL
  if (cn>1) for(k in 1:nrow(test.list2)){
    t<-test.list2[k,]
    lab<-paste(t[1],t[2],sep="--")
    test.list3<-rbind(test.list3,lab)}
  
  if (cn==1) test.list3<-c("1--2")
  
  lab2<-c(rep("MD",cn),rep("MDC",cn),rep("MNN",cn),rep("ECC",cn))
  test.list4<-paste(lab2,test.list3,sep=".")
  
  head<-c("iteration",test.list4)
  p.table<-rbind(head,line1,p.table)
  
}


test1<-permute(lm.gp,lm.gp.red,Y,group,999) # run the permutation test
test1

# results can be opened in e.g.,
# excel. Find the P-value by measuring the rank of the first value in the
# distribution of values in the same column.  
# p.table.out<-test1
# write.csv(p.table.out,"ACRO_Turner_results.csv") 
# Alternatvely, more R-code can
# be generated to calculate P-values by using a function that ranks observed values
# see below

# P-Values
h<-test1[1,];h<-h[-1]
t<-test1[-1,];t<-t[,-1]
f<-function(t){
  r<-rank(t)
  p<-r[1]
  c<-length(t)
  pv<-(c-p+1)/c}
p.value<-apply(t,2,f)
names(p.value)<-h

p.value

#AP 2021 Nursery host =1 sym = 2
# MD.1--2 MDC.1--2 MNN.1--2 ECC.1--2 
# 0.001    0.017    0.433    0.082 

#AP 2022 Nursery
# MD.1--2 MDC.1--2 MNN.1--2 ECC.1--2 
#0.0020   0.4860   0.3430   0.7915

#AP Am
# MD.1--2 MDC.1--2 MNN.1--2 ECC.1--2 
# 0.007    0.045    0.026    0.216 

#OF 2021 Nursery
# MD.1--2 MDC.1--2 MNN.1--2 ECC.1--2 
# 0.0010   0.0015   0.4060   0.0160 

#OF CI
# MD.1--2 MDC.1--2 MNN.1--2 ECC.1--2 
# 0.002    0.304    0.134    0.958 

#ADDENDUM: HOTELLING'S T2
#Each group comparison one at a time

#Comparing groups 1 (host) and 2 (symbiont) 
gp.m.dif<-gp.m[1,]-gp.m[2,] # vector for difference between means
gn<-tapply(group,group,length) # group sizes
e<-resid(lm.gp)
E<-t(e)%*%e
n<-nrow(e)
k<-lm.gp$rank
V<-(1/(n-k))*E # This is the pooled within-group vcv

d<-gp.m.dif; dim(d)<-c(1,length(d))

D<-d%*%solve(V)%*%t(d) # Squared Mahalanobis Distance

H.T2<-(gn[1]*gn[2])/(gn[1]+gn[2])*D # Hotelling T2

F<-(gn[1]+gn[2]-2-1)/((gn[1]+gn[2])*2)*H.T2 # Convert to an F value

P<-df(F,2,(gn[1]+gn[2]-2-1)) # P-value

P
#AP 2021 Nursery 
#p = 0.001516803, F = 10.1478, T2 = 27.0607

#AP 2022 Nursery
#p = 0.0008330837, F = 11.8338, T2 = 31.5568

#AP Am
#p = 0.005342675, F = 7.1506, T2 = 19.06819

#OF 2021 Nursery
#p = 0.003518716, F = 8.0697, T2 = 21.5192

#OF CI
#p = 9.914311e-06, F = 41.79077, T2 = 119.4022




#_______________________Bayesian Mixing Model________________________________

library(MixSIAR)
#_______________________________
#Example code and datasets available in package library
#find MixSIAR directory
#mixsiar.dir <- find.package("MixSIAR")
#folder with example script files, pull up mantis shrimp R code for practice
#paste0(mixsiar.dir,"/example_scripts")
#_______________________________

#Real data code from mantis shrimp dataset, adapted for our corals

# Load mixture data, i.e. your:
#    Consumer isotope values (trophic ecology / diet)
#    Mixed sediment/water tracer values (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with mix/consumer data
# 'iso_names' - column headings of the tracers/isotopes you'd like to use
# 'random_effects' - column headings of any random effects
# 'cont_effects' - column headings of any continuous effects

# 'iso_names', 'random_effects', and 'cont_effects' can be a subset of your columns
#   i.e. have 3 isotopes in file but only want MixSIAR to use 2,
#   or have data by Region and Pack but only want MixSIAR to use Region

# To run on your data, replace the system.file call with the path to your file
#mix.filename <- system.file("extdata", "mantis_consumer.csv", package = "MixSIAR")#example data
mix.filename <- system.file("extdata", "AP_MotePOR_PecentContribution_consumer.csv", package = "MixSIAR")
mix.filename <- system.file("extdata", "OF_MotePOR_PecentContribution_consumer.csv", package = "MixSIAR")
# I placed my file in the MixSIAR extdata folder as a work around, 
#since read.csv was not loading the correct file type for the mix code

# Load mixture data
?load_mix_data

mix <- load_mix_data(filename = mix.filename, 
                     iso_names = c("d13C","d15N"), 
                     factors = c("Treatment","Symbiont"), 
                     fac_random = c(FALSE,TRUE), 
                     fac_nested = c(FALSE, FALSE), 
                     cont_effects = NULL)

mix <- load_mix_data(filename = mix.filename, 
                     iso_names = c("d13C","d15N"), 
                     factors = c("Treatment"), 
                     fac_random = c(FALSE), 
                     fac_nested = c(FALSE), 
                     cont_effects = NULL)

mix <- load_mix_data(filename = mix.filename, 
                     iso_names = c("d13C","d15N"), 
                     factors = c("Symbiont"), 
                     fac_random = c(FALSE), 
                     fac_nested = c(FALSE), 
                     cont_effects = NULL)

################################################################################
# Load source data, i.e. your:
#    Source isotope values (trophic ecology / diet)
#    Sediment/water source tracer values (sediment/hydrology fingerprinting)

# 'filename': name of the CSV file with source data
# 'source_factors': column headings of random/fixed effects you have source data by
# 'conc_dep': TRUE or FALSE, do you have concentration dependence data in the file?
# 'data_type': "means" or "raw", is your source data as means+SD, or do you have raw data

# To run on your data, replace the system.file call with the path to your file
#source.filename <- system.file("extdata", "mantis_source.csv", package = "MixSIAR")#example
source.filename <- system.file("extdata", "MotePOR_AllsourcesTrea_AP_PercentContribution_sources.csv", package = "MixSIAR")
source.filename <- system.file("extdata", "MotePOR_AllsourcesSym_AP_PercentContribution_sources.csv", package = "MixSIAR")
source.filename <- system.file("extdata", "MotePOR_AllsourcesTrea_OF_PercentContribution_sources.csv", package = "MixSIAR")

# Load source data
source <- load_source_data(filename = source.filename, 
                           source_factors = "Treatment", 
                           conc_dep = FALSE, 
                           data_type = "means", mix)

source <- load_source_data(filename = source.filename, 
                           source_factors = "Symbiont", 
                           conc_dep = FALSE, 
                           data_type = "means", mix)


################################################################################
# Load discrimination data, i.e. your:
#    Trophic Enrichment Factor (TEF) / fractionation values (trophic ecology/diet)
#    xxxxxxxx (sediment/hydrology fingerprinting)

# 'filename' - name of the CSV file with discrimination data

# To run on your data, replace the system.file call with the path to your file
discr.filename <- system.file("extdata", "MotePOR_AP_discrimination.csv", package = "MixSIAR")

# Load discrimination data
discr <- load_discr_data(filename=discr.filename, mix)

#####################################################################################
# Make isospace plot
# Are the data loaded correctly?
# Is your mixture data in the source polygon?
# Are one or more of your sources confounded/hidden?

# 'filename' - name you'd like MixSIAR to save the isospace plot as 
#              (extension will be added automatically)
# 'plot_save_pdf' - TRUE or FALSE, should MixSIAR save the plot as a .pdf?
# 'plot_save_png' - TRUE or FALSE, should MixSIAR save the plot as a .png?

plot_data(filename="isospace_plot", 
          plot_save_pdf=TRUE,
          plot_save_png=FALSE,
          mix,source,discr)

# If 2 isotopes/tracers, calculate normalized surface area of the convex hull polygon(s)
#   *Note 1: discrimination SD is added to the source SD (see calc_area.r for details)
#   *Note 2: If source data are by factor (as in wolf ex), computes area for each polygon
#             (one for each of 3 regions in wolf ex)
# if(mix$n.iso==2) calc_area(source=source,mix=mix,discr=discr)

################################################################################
# Define your prior, and then plot using "plot_prior"
#   RED = your prior
#   DARK GREY = "uninformative"/generalist (alpha = 1)
#   LIGHT GREY = "uninformative" Jeffrey's prior (alpha = 1/n.sources)

# default "UNINFORMATIVE" / GENERALIST prior (alpha = 1)
alpha.unif <- rep(1, source$n.sources)
plot_prior(alpha.prior=alpha.unif, 
           source=source,
           filename="prior_uninf")

################################################################################
# Write JAGS model file (define model structure)
# Model will be saved as 'model_filename' ("MixSIAR_model.txt" is default,
#    but may want to change if in a loop)

# There are 3 error term options available:
#   1. Residual * Process (resid_err = TRUE, process_err = TRUE)
#   2. Residual only (resid_err = TRUE, process_err = FALSE)
#   3. Process only (resid_err = FALSE, process_err = TRUE)

# 'model_filename': don't need to change, unless you create many different models
# 'resid_err': include residual error in the model?
# 'process_err': include process error in the model?

#  *Note: If you have only 1 mix datapoint, you have no information about the 
#         mixture/consumer variability. In this case, we ues the original MixSIR
#         error model (which does not fit a residual error term).
#         This is the same behavior as 'siarsolo' in SIAR.

model_filename <- "MixSIAR_model.txt"
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

################################################################################
# Run model
# JAGS output will be saved as 'jags.1'

# MCMC run options:
# run <- "test"       # chainLength=1000, burn=500, thin=1, chains=3, calcDIC=TRUE
# run <- "very short" # chainLength=10000, burn=5000, thin=5, chains=3, calcDIC=TRUE
# run <- "short"      # chainLength=50000, burn=25000, thin=25, chains=3, calcDIC=TRUE
# run <- "normal"     # chainLength=100000, burn=50000, thin=50, chains=3, calcDIC=TRUE
# run <- "long"       # chainLength=300000, burn=200000, thin=100, chains=3, calcDIC=TRUE
# run <- "very long"  # chainLength=1000000, burn=500000, thin=500, chains=3, calcDIC=TRUE
# run <- "extreme"    # chainLength=3000000, burn=1500000, thin=500, chains=3, calcDIC=TRUE

# Can also set custom MCMC parameters
# run <- list(chainLength=200000, burn=150000, thin=50, chains=3, calcDIC=TRUE)

# Good idea to use 'test' first to check if
#   1) the data are loaded correctly, and 
#   2) the model is specified correctly
# jags.1 <- run_model(run="test", mix, source, discr, model_filename, alpha.prior = alpha.spec)

# After a test run works, increase the MCMC run to a value that may converge
jags.spec <- run_model(run="normal", mix, source, discr, model_filename, alpha.prior = alpha.unif)#changed to alpha.unif for generalist troubleshoot

################################################################################
# Process JAGS output

# Choose output options (see ?output_options for details)
output_options <- list(summary_save = TRUE,                 
                       summary_name = "summary_statistics", 
                       sup_post = FALSE,                    
                       plot_post_save_pdf = TRUE,           
                       plot_post_name = "posterior_density",
                       sup_pairs = FALSE,             
                       plot_pairs_save_pdf = TRUE,    
                       plot_pairs_name = "pairs_plot",
                       sup_xy = TRUE,           
                       plot_xy_save_pdf = FALSE,
                       plot_xy_name = "xy_plot",
                       gelman = TRUE,
                       heidel = FALSE,  
                       geweke = TRUE,   
                       diag_save = TRUE,
                       diag_name = "diagnostics",
                       indiv_effect = FALSE,       
                       plot_post_save_png = FALSE, 
                       plot_pairs_save_png = FALSE,
                       plot_xy_save_png = FALSE,
                       diag_save_ggmcmc = FALSE)

# Create diagnostics, summary statistics, and posterior plots
output_JAGS(jags.spec, mix, source, output_options)

#use combine to calculate Bayesian Credible Intervals
original <- combine_sources(jags.spec, mix, source, alpha.unif, 
                            groups=list(POM="POM",Symbiodiniaceae="Symbiodiniaceae",Zooplankton="Zooplankton"))

#Plot 50% (thin) and 95% (thick) Bayesian Credible Intervals - edit this in illustrator
plot_intervals(
  original,
  toplot = "fac1",
  levels = NULL,
  groupby = "source",
  savepdf = FALSE,
  filename = "post_intervals")

#get numerical representation of Bayesian Credible Intervals - not working
library("bayesplot")
mcmc_intervals_data(
  jags.spec,
  prob = 0.5,
  prob_outer = 0.95,
  point_est = c("mean"),
)


#calculate overlap of proportion of diets for corals by treatment
library(overlapping)
?overlap
set.seed(1)

#Model with ALL SOURCES
OFsym <- list(c(0.407, 0.467, 0.699, 0.836, 0.915,  0.974,  0.983), #nursery
              c(0.572, 0.681, 0.866, 0.933, 0.969,  0.991,  0.994)) #nearshore
OFzoo <- list(c(0.002, 0.004, 0.021, 0.055, 0.122,  0.307,  0.383), #nursery
              c(0.001, 0.001, 0.008, 0.021, 0.057,  0.169,  0.253)) #nearshore
OFpom <- list(c(0.003, 0.005, 0.030, 0.077, 0.172,  0.383,  0.457), #nursery
              c(0.001, 0.002, 0.010, 0.029, 0.072,  0.214,  0.2773)) #nearshore
#AP 2021-2022 Nursery
APsym <-list(c(0.522, 0.576, 0.735, 0.820, 0.891, 0.958, 0.971), #2021
             c(0.197, 0.262, 0.509, 0.695, 0.874, 0.978, 0.988)) #2022
APzoo <-list(c(0.002, 0.005, 0.024, 0.060, 0.122, 0.264, 0.331), #2021
             c(0.001, 0.002, 0.018, 0.064, 0.183, 0.451, 0.576)) #2022
APpom <-list(c(0.004, 0.008, 0.042, 0.092, 0.160, 0.277, 0.320), #2021
             c(0.002, 0.004, 0.037, 0.135, 0.329, 0.592, 0.662)) #2022
#AP 2021-Offshore
APsym <-list(c(0.522, 0.576, 0.735, 0.820, 0.891, 0.958, 0.971), #2021
             c(0.584, 0.669, 0.866, 0.927, 0.964, 0.989, 0.993)) #Offshore
APzoo <-list(c(0.002, 0.005, 0.024, 0.060, 0.122, 0.264, 0.331), #2021
             c(0.001, 0.001, 0.008, 0.022, 0.056, 0.168, 0.227)) #Offshore
APpom <-list(c(0.004, 0.008, 0.042, 0.092, 0.160, 0.277, 0.320), #2021
             c(0.001, 0.002, 0.012, 0.033, 0.076, 0.222, 0.309)) #Offshore
#AP 2022-Offshore
APsym <-list(c(0.197, 0.262, 0.509, 0.695, 0.874, 0.978, 0.988), #2022
             c(0.584, 0.669, 0.866, 0.927, 0.964, 0.989, 0.993)) #Offshore
APzoo <-list(c(0.001, 0.002, 0.018, 0.064, 0.183, 0.451, 0.576), #2022
             c(0.001, 0.001, 0.008, 0.022, 0.056, 0.168, 0.227)) #Offshore
APpom <-list(c(0.002, 0.004, 0.037, 0.135, 0.329, 0.592, 0.662), #2022
             c(0.001, 0.002, 0.012, 0.033, 0.076, 0.222, 0.309)) #Offshore
overlap(APpom, type = "2", pairsOverlap = TRUE, plot = TRUE)


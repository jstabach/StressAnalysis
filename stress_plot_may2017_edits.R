#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Project: SHO Stress and Behavior Analysis
# Date: 7 November 2016 - updated 26 May 2017
# Author: Stephanie Cunningham
# Description: Plotting stress models

#**********************************************************************************************************************************
#**********************************************************************************************************************************

# Remove from memory
rm(list=ls())

library(tidyr)
library(dplyr)
library(sme)

# Set working directory
#setwd("/Volumes/STEPH 1TB/SCBI Work/SHO Behavior and Stress/SCBI/FecalGlucoCorticoids/")
setwd("D:/Jared/Projects/Oryx/StressAnalysis/Stephanie_Analysis/")
stress <- read.csv("stress_8nov.csv")

stress <- gather(stress, "ID", "fGC", 2:14)

stress$ID <- as.character(stress$ID)
Sex <- stress[,2]
Sex[Sex=="Tex" | Sex=="Lewis" | Sex=="DrBob"] <- "Male"
Sex[Sex=="Cat" | Sex=="Chari" | Sex=="Loretta" | Sex=="Niamey" | Sex=="Violet" | Sex=="Jena" | Sex=="Ruby" | Sex=="Shadow" | Sex=="Scout" | Sex=="Bamako"] <- "Female"
stress <- cbind(stress, Sex)
Treatment <- stress[,2]
Treatment[Treatment=="Tex" | Treatment=="Lewis" | Treatment=="Violet" | Treatment=="Jena"] <- "Vectronic" 
Treatment[Treatment=="Niamey" | Treatment=="Cat" | Treatment=="Loretta" | Treatment=="DrBob"  | Treatment=="Chari"] <- "Control"
Treatment[ Treatment=="Ruby" | Treatment=="Shadow" | Treatment=="Scout" | Treatment=="Bamako"] <- "ATS"
stress <- cbind(stress, Treatment)


bam <- subset(stress, ID=="Bamako")  
RelDay <- c(seq(-7,54,1))
bam <- cbind(bam, RelDay)

tex <- subset(stress, ID=="Tex") 
RelDay <- c(seq(-19,42,1))
tex <- cbind(tex, RelDay)

vio <- subset(stress, ID=="Violet") 
RelDay <- c(seq(-7,54,1))
vio <- cbind(vio, RelDay)

lew <- subset(stress, ID=="Lewis") 
RelDay <- c(seq(-19,42,1))
lew <- cbind(lew, RelDay)

drb <- subset(stress, ID=="DrBob") 
RelDay <- c(seq(-19,42,1))
drb <- cbind(drb, RelDay)

jen <- subset(stress, ID=="Jena")
RelDay <- c(seq(-7,54,1))
jen <- cbind(jen, RelDay)

sco <- subset(stress, ID=="Scout")
RelDay <- c(seq(-7,54,1))
sco <- cbind(sco, RelDay)

rub <- subset(stress, ID=="Ruby")
RelDay <- c(seq(-7,54,1))
rub <- cbind(rub, RelDay)

sha <- subset(stress, ID=="Shadow")
RelDay <- c(seq(-7,54,1))
sha <- cbind(sha, RelDay)

nia <- subset(stress, ID=="Niamey")
RelDay <- c(seq(-7,54,1))
nia <- cbind(nia, RelDay)

cha <- subset(stress, ID=="Chari")
RelDay <- c(seq(-7,54,1))
cha <- cbind(cha, RelDay)

cat <- subset(stress, ID=="Cat")
RelDay <- c(seq(-7,54,1))
cat <- cbind(cat, RelDay)

lor <- subset(stress, ID=="Loretta")
RelDay <- c(seq(-7,54,1))
lor <- cbind(lor, RelDay)

stress <- rbind(lor, cat, cha, nia, sha, rub, jen, drb, sco, lew, vio, tex, bam)
stress$Date <- as.character(stress$Date)
stress$Date <- as.POSIXct(stress$Date, format="%d-%b-%y")
N.order <-order(stress$Date,decreasing=FALSE)
stress <- stress[N.order,]

stress <- stress[,c(1,6,2,4,5,3)]

#write.csv(stress, "stress_data_nov8.csv")

# Now that we've got this more organized and clean, can start subsetting and reapplying the analysis.

strsub <- subset(stress, RelDay>=-8  & RelDay<=10)

mcllr <- subset(strsub, Sex=="Male" & Treatment=="Vectronic")
fcllr <- subset(strsub, Sex=="Female" & (Treatment=="Vectronic" | Treatment=="ATS"))
mcntr <- subset(strsub, ID=="DrBob")
fcntr <- subset(strsub, Sex=="Female" & Treatment=="Control")
ccllr <- subset(strsub, Treatment=="Vectronic" | Treatment=="ATS")
ccntr <- subset(strsub, Treatment=="Control")

# quick plot to visualize data
par(mfrow=c(1,1))
plot(mcllr$RelDay, mcllr$fGC, pch=16, ylim=c(25,70), col="blue", ylab="Stress", xlab="Relative Days")
par(new=TRUE)
plot(fcllr$RelDay, fcllr$fGC, pch=16, ylim=c(25,70), col="red", axes=FALSE, ylab=" ", xlab=" ")
par(new=TRUE)
plot(mcntr$RelDay, mcntr$fGC, pch=16, ylim=c(25,70), col="green", axes=FALSE, ylab=" ", xlab=" ")
par(new=TRUE)
plot(fcntr$RelDay, fcntr$fGC, pch=16, ylim=c(25,70), col="purple", axes=FALSE, ylab=" ", xlab=" ")
par(new=FALSE)

# Best-fitting model for each group
# COMBINED COLLARS - MODEL 5........................................................................
x1 <- (ccllr$RelDay %in% 3:10)
x2 <- x1 * ccllr$RelDay
pw1 <- lm(ccllr$fGC~x1 + x2)
(summary(pw1))
AICc(pw1)

# I would have wrote the model statement something like this.....model doesn't end up being as good as your best model (based on AIC)
pw <- lm(ccllr$fGC ~ ccllr$RelDay*(ccllr$RelDay < 3) + ccllr$RelDay*(ccllr$RelDay > 3))
(summary(pw))
AICc(pw)

c1 <- seq(-8,10,1)
cnt1.rev <- rev(c1)
test <- loess(formula=ccllr$fGC~ccllr$RelDay, span=0.75, degree=1)
pr1 <- predict(test, c1, se=T)

up1 <- c(pr1$fit+1.96*pr1$se.fit)
low1 <- c(pr1$fit-1.96*pr1$se.fit)

plot(ccllr$RelDay,ccllr$fGC, type="n", ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
polygon(c(c1,cnt1.rev), c(up1, rev(low1)), col="grey80", border=NA)
lines(c1, pr1$fit, lty=1, col="white", lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(ccllr$RelDay, ccllr$fGC, pch=16, ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2))
axis(2, at=seq(20,70,10))
abline(v=3, lty=3)
curve((coef(pw)[1] + coef(pw)[3]) + (coef(pw)[2]+coef(pw)[5])*x, add=T, from=-8, to=3,lwd=2)
curve((coef(pw)[1] + coef(pw)[4]) + coef(pw)[2]*x, add=T, from=3, to=10,lwd=2)

plot(ccllr$RelDay, ccllr$fGC, type="n", ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
polygon(c(c1,cnt1.rev), c(up1, rev(low1)), col="grey70", border=NA)
lines(c1, pr1$fit, lty=1, col="gray30", lwd=1, xlim=c(-8,10))
par(new=TRUE)
plot(ccllr$RelDay, ccllr$fGC, pch=16, ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2))
axis(2, at=seq(20,70,10))
abline(v=3, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw1)[1], x1=3, y1=coef(pw1)[1], col="black", lwd=2)
#segments(x0=3, y0=coef(pw1)[1]+coef(pw1)[2], x1=10, y1=y2 + 10*coef(pw1)[3], col="black", lwd=2) # This part does not work
y2 <- coef(pw1)[1]+coef(pw1)[2]
segments(x0=3, y0=coef(pw1)[1]+coef(pw1)[2], x1=10, y1=y2 + 10*coef(pw1)[3], col="black", lwd=2) 
par(new=FALSE)

# COMBINED CONTROLS - MODEL 2........................................................................
x1 <- (ccntr$RelDay %in% 1:10)  
pw2 <- lm(ccntr$fGC~x1)
(summary(pw2))
AICc(pw1)

count <- seq(-8,10,1)
cnt.rev <- rev(count)
test1 <- loess(formula=ccntr$fGC~ccntr$RelDay, span=0.75, degree=1)
pr <- predict(test1, count, se=T)
# lines(count, pr$fit, lty=1, col="blue", lwd=2, xlim=c(-8,10))
# lines(pr$fit-1.96*pr$se.fit, lty=2, lwd=2, col="green")
# lines(pr$fit+1.96*pr$se.fit, lty=2, lwd=2, col="green")

up <- c(pr$fit+1.96*pr$se.fit)
low <- c(pr$fit-1.96*pr$se.fit)

plot(ccntr$RelDay, ccntr$fGC,type="n", ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
polygon(c(count,cnt.rev), c(up, rev(low)), col="grey70", border=NA)
lines(count, pr$fit, lty=1, col="gray30", lwd=1, xlim=c(-8,10))
par(new=TRUE)
plot(ccntr$RelDay, ccntr$fGC, pch=16, ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2))
axis(2, at=seq(20,70,10))
abline(v=1, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw2)[1], x1=1, y1=coef(pw2)[1], col="black", lwd=2)
segments(x0=1, y0=coef(pw2)[1]+coef(pw2)[2], x1=10, y1=coef(pw2)[1]+coef(pw2)[2], col="black", lwd=2)
par(new=FALSE)

# COLLARED MALES - MODEL 5........................................................................
x1 <- (mcllr$RelDay %in% 3:10)
x2 <- x1 * mcllr$RelDay
pw3 <- lm(mcllr$fGC~x1 + x2)
(summary(pw3))
AICc(pw3)

mlr <- mcllr[,c(2,6)]
mlr <- filter(mlr, fGC>0)
l3 <- loess.sd(mlr, nsigma=1.96, span=0.75, degree=1)

c2 <- seq(-8,10,1)
cnt2.rev <- rev(c2)
test2 <- loess(formula=mcllr$fGC~mcllr$RelDay, span=0.75, degree=1)
pr2 <- predict(test2, c2, se=T)
# lines(c2, pr2$fit, lty=1, col="blue", lwd=2, xlim=c(-8,10))
# lines(pr2$fit-1.96*pr2$se.fit, lty=2, lwd=2, col="green")
# lines(pr2$fit+1.96*pr2$se.fit, lty=2, lwd=2, col="green")

up2 <- c(pr2$fit+1.96*pr2$se.fit)
low2 <- c(pr2$fit-1.96*pr2$se.fit)

plot(mcllr$RelDay, mcllr$fGC, type="n", ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
polygon(c(c2,cnt2.rev), c(up2, rev(low2)), col="grey70", border=NA)
lines(c2, pr2$fit, lty=1, col="gray30", lwd=1, xlim=c(-8,10))
par(new=TRUE)
plot(mcllr$RelDay, mcllr$fGC, pch=16, ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2))
axis(2, at=seq(20,70,10))
abline(v=3, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw3)[1], x1=3, y1=coef(pw3)[1], col="black", lwd=2)
segments(x0=3, y0=coef(pw3)[1]+coef(pw3)[2], x1=10, y1=coef(pw3)[1]+coef(pw3)[2] + 7*coef(pw3)[3], col="black", lwd=2)
par(new=FALSE)


# CONTROL MALE - MODEL 2........................................................................
x1 <- (mcntr$RelDay %in% 1:10)  
pw4 <- lm(mcntr$fGC~x1)
(summary(pw4))
AICc(pw4)

c3 <- seq(-8,10,1)
cnt3.rev <- rev(c3)
test3 <- loess(formula=mcntr$fGC~mcntr$RelDay, span=0.75, degree=1)
pr3 <- predict(test3, c3, se=T)
# lines(c3, pr2$fit, lty=1, col="blue", lwd=2, xlim=c(-8,10))
# lines(pr3$fit-1.96*pr3$se.fit, lty=2, lwd=2, col="green")
# lines(pr3$fit+1.96*pr3$se.fit, lty=2, lwd=2, col="green")

up3 <- c(pr3$fit+1.96*pr3$se.fit)
low3 <- c(pr3$fit-1.96*pr3$se.fit)

plot(mcntr$RelDay, mcntr$fGC, type="n", ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
polygon(c(c3,cnt3.rev), c(up3, rev(low3)), col="grey70", border=NA)
lines(c3, pr3$fit, lty=1, col="gray30", lwd=1, xlim=c(-8,10))
par(new=TRUE)
plot(mcntr$RelDay, mcntr$fGC, pch=16, ylim=c(20,70), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2))
axis(2, at=seq(20,70,10))
abline(v=1, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw4)[1], x1=1, y1=coef(pw4)[1], col="black", lwd=2)
segments(x0=1, y0=coef(pw4)[1]+coef(pw4)[2], x1=10, y1=coef(pw4)[1]+coef(pw4)[2], col="black", lwd=2)
par(new=FALSE)


# COLLARED FEMALES - MODEL 2........................................................................
x1 <- (fcllr$RelDay %in% 2:10)  
pw5 <- lm(fcllr$fGC~x1)
(summary(pw5))
AICc(pw5)

c4 <- seq(-8,10,1)
cnt4.rev <- rev(c4)
test4 <- loess(formula=fcllr$fGC~fcllr$RelDay, span=0.75, degree=1)
pr4 <- predict(test4, c4, se=T)
# lines(c4, pr2$fit, lty=1, col="blue", lwd=2, xlim=c(-8,10))
# lines(pr4$fit-1.96*pr4$se.fit, lty=2, lwd=2, col="green")
# lines(pr4$fit+1.96*pr4$se.fit, lty=2, lwd=2, col="green")

up4 <- c(pr4$fit+1.96*pr4$se.fit)
low4 <- c(pr4$fit-1.96*pr4$se.fit)

plot(fcllr$RelDay, fcllr$fGC, type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab=" ", xlab=" ")
polygon(c(c4,cnt4.rev), c(up4, rev(low4)), col="grey85", border=NA)
lines(c4, pr4$fit, lty=2, col="gray30", lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(fcllr$RelDay, fcllr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10), cex.lab=1.4, axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2), cex.axis=1.3)
axis(2, at=seq(20,70,10), cex.axis=1.3)
abline(v=1, lty=3, lwd=1)
segments(x0=-6, y0=coef(pw5)[1], x1=1, y1=coef(pw5)[1], col="black", lwd=2)
segments(x0=1, y0=coef(pw5)[1]+coef(pw5)[2], x1=10, y1=coef(pw5)[1]+coef(pw5)[2], col="black", lwd=2)
par(new=FALSE)

# CONTROL FEMALES - MODEL 5........................................................................
x1 <- (fcntr$RelDay %in% 3:10)
x2 <- x1 * fcntr$RelDay
pw6 <- lm(fcntr$fGC~x1 + x2)
(summary(pw6))
AICc(pw6)

#c5 <- seq(-8,10,1)
c5 <- seq(-7,8,1)
cnt5.rev <- rev(c5)
test5 <- loess(formula=fcntr$fGC~fcntr$RelDay, span=0.75, degree=1)
pr5 <- predict(test5, c5, se=T)
# lines(c5, pr5$fit, lty=1, col="blue", lwd=2, xlim=c(-8,10))
# lines(c5,pr5$fit-1.96*pr5$se.fit, lty=2, lwd=2, col="green")
# lines(c5,pr5$fit+1.96*pr5$se.fit, lty=2, lwd=2, col="green")

up5 <- c(pr5$fit+1.96*pr5$se.fit)
low5 <- c(pr5$fit-1.96*pr5$se.fit)

#t <- c(48.560,47.790,47.068,46.502,46.143,45.849,45.859,46.224,46.597,48.012,50.220,49.701,47.869,44.881,39.125,36.908)
#b <- c(25.065,31.581,35.885,39.799,40.493,40.184,38.615,36.433,35.559,35.874,36.221,35.844,35.049,34.223,33.304,32.283)
u <- seq(-7,8,1)
u.rev <- rev(u)
#polygon(c(u, u.rev), c(t,b))
polygon(c(u, u.rev), c(up5,rev(low5)))

plot(fcntr$RelDay, fcntr$fGC, type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
#polygon(c(u, u.rev), c(t,b), col="grey85", border=NA)
polygon(c(u, u.rev), c(up5,rev(low5)), col="grey85", border=NA)
lines(c5, pr5$fit, lty=2,  lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(fcntr$RelDay, fcntr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2))
axis(2, at=seq(20,70,10))
abline(v=3, lty=3, lwd=1)
segments(x0=-7, y0=coef(pw6)[1], x1=3, y1=coef(pw6)[1], col="black", lwd=2)
segments(x0=3, y0=coef(pw6)[1]+coef(pw6)[2], x1=8, y1=coef(pw6)[1]+coef(pw6)[2] + 8*coef(pw6)[3], col="black", lwd=2)
par(new=FALSE)

#*****************************************************************************************************
# Final Plot

par(mfrow=c(3,2), oma=c(0,0,0,0), mar=c(4.5,4.5,0,0))
# Collared Males
plot(mcllr$RelDay, mcllr$fGC, type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab=" ", xlab=" ")
polygon(c(c2,cnt2.rev), c(up2, rev(low2)), col="grey85", border=NA)
lines(c2, pr2$fit, lty=2, col="gray30", lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(mcllr$RelDay, mcllr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10),cex.lab=1.4, axes=FALSE, ylab="Stress (ng/g)", xlab=" ")
axis(1, at=seq(-8,10,2), cex.axis=1.3)
axis(2, at=seq(20,70,10), cex.axis=1.3)
abline(v=3, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw3)[1], x1=3, y1=coef(pw3)[1], col="black", lwd=2)
segments(x0=3, y0=coef(pw3)[1]+coef(pw3)[2], x1=10, y1=coef(pw3)[1]+coef(pw3)[2] + 7*coef(pw3)[3], col="black", lwd=2)
text(-8,68, "A", cex=1.4)
par(new=FALSE)

# Control Male
plot(mcntr$RelDay, mcntr$fGC, type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab=" ", xlab=" ")
polygon(c(c3,cnt3.rev), c(up3, rev(low3)), col="grey85", border=NA)
lines(c3, pr3$fit, lty=2, col="gray30", lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(mcntr$RelDay, mcntr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10), cex.lab=1.4,axes=FALSE, ylab=" ", xlab=" ")
axis(1, at=seq(-8,10,2), cex.axis=1.3)
axis(2, at=seq(20,70,10), cex.axis=1.3)
abline(v=1, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw4)[1], x1=1, y1=coef(pw4)[1], col="black", lwd=2)
segments(x0=1, y0=coef(pw4)[1]+coef(pw4)[2], x1=10, y1=coef(pw4)[1]+coef(pw4)[2], col="black", lwd=2)
text(-8,68, "B", cex=1.4)
par(new=FALSE)

# Collared Females
plot(fcllr$RelDay, fcllr$fGC, type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab=" ", xlab=" ")
polygon(c(c4,cnt4.rev), c(up4, rev(low4)), col="grey85", border=NA)
lines(c4, pr4$fit, lty=2, col="gray30", lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(fcllr$RelDay, fcllr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10),cex.lab=1.4, axes=FALSE, ylab="Stress (ng/g)", xlab=" ")
axis(1, at=seq(-8,10,2), cex.axis=1.3)
axis(2, at=seq(20,70,10), cex.axis=1.3)
abline(v=1, lty=3, lwd=1)
segments(x0=-6, y0=coef(pw5)[1], x1=1, y1=coef(pw5)[1], col="black", lwd=2)
segments(x0=1, y0=coef(pw5)[1]+coef(pw5)[2], x1=10, y1=coef(pw5)[1]+coef(pw5)[2], col="black", lwd=2)
text(-8,68, "C", cex=1.4)
par(new=FALSE)

# Control Females
plot(fcntr$RelDay, fcntr$fGC, type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab=" ", xlab=" ")
#polygon(c(u, u.rev), c(t,b), col="grey85", border=NA)
polygon(c(u, u.rev), c(up5,rev(low5)), col="grey85", border=NA)
lines(c5, pr5$fit, lty=2,  lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(fcntr$RelDay, fcntr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10),cex.lab=1.4, axes=FALSE, ylab=" ", xlab=" ")
axis(1, at=seq(-8,10,2), cex.axis=1.3)
axis(2, at=seq(20,70,10), cex.axis=1.3)
abline(v=3, lty=3, lwd=1)
segments(x0=-7, y0=coef(pw6)[1], x1=3, y1=coef(pw6)[1], col="black", lwd=2)
segments(x0=3, y0=51.9, x1=10, y1=23.2, col="black", lwd=2) # Not sure where this comes from?
text(-8,68, "D", cex=1.4)
par(new=FALSE)

# Combined Collars
y2 <- coef(pw1)[1]+coef(pw1)[2]

plot(ccllr$RelDay, ccllr$fGC, type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab=" ", xlab=" ")
#polygon(c(c1,cnt1.rev), c(up1, rev(low1)), col="grey85", border=NA)
polygon(c(u, u.rev), c(up5,rev(low5)), col="grey85", border=NA)
lines(c1, pr1$fit, lty=2, col="gray30", lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(ccllr$RelDay, ccllr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10),cex.lab=1.4, axes=FALSE, ylab="Stress (ng/g)", xlab="Relative Days")
axis(1, at=seq(-8,10,2), cex.axis=1.3)
axis(2, at=seq(20,70,10), cex.axis=1.3)
abline(v=3, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw1)[1], x1=3, y1=coef(pw1)[1], col="black", lwd=2)
segments(x0=3, y0=coef(pw1)[1]+coef(pw1)[2], x1=10, y1=y2 + 10*coef(pw1)[3], col="black", lwd=2)
text(-8,68, "E", cex=1.4)
par(new=FALSE)

# Combined Controls
plot(ccntr$RelDay, ccntr$fGC,type="n", ylim=c(20,70), xlim=c(-8,10), axes=FALSE, ylab=" ", xlab=" ")
polygon(c(count,cnt.rev), c(up, rev(low)), col="grey85", border=NA)
lines(count, pr$fit, lty=2, col="gray30", lwd=2, xlim=c(-8,10))
par(new=TRUE)
plot(ccntr$RelDay, ccntr$fGC, pch=16, ylim=c(20,70), xlim=c(-8,10),cex.lab=1.4, axes=FALSE, ylab=" ", xlab="Relative Days")
axis(1, at=seq(-8,10,2), cex.axis=1.3)
axis(2, at=seq(20,70,10), cex.axis=1.3)
abline(v=1, lty=3, lwd=1)
segments(x0=-8, y0=coef(pw2)[1], x1=1, y1=coef(pw2)[1], col="black", lwd=2)
segments(x0=1, y0=coef(pw2)[1]+coef(pw2)[2], x1=10, y1=coef(pw2)[1]+coef(pw2)[2], col="black", lwd=2)
text(-8,68, "F", cex=1.4)
par(new=FALSE)









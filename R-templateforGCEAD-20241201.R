# Author: Kelsey J.R.P. Byers (kbyers@alum.mit.edu)
# License: MIT License

# Copyright 2024 Kelsey J.R.P. Byers

# Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the “Software”), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# An R script to enable quantitative data processing of GC-EAD data. See the preprint at Byers & Jacobs (2024) Quantitative analysis of gas chromatography-coupled electroantennographic detection (GC-EAD) of plant volatiles by insects. bioRxiv BIORXIV/2024/626223.


# Import libraries
library(dplyr)

# Import data
insect<-read.csv("") # plug in appropriate GcEAD CSV file export
fid<-read.csv("",header=FALSE) # plug in appropriate OpenLab CDS CSV file export

# Create appropriate vectors for data
fidtime<-fid[,1]
fidamp<-fid[,2]

eadamp<-insect[,2]
maxeadtime<-(length(eadamp)-1)/6000
eadtime<-seq(0,maxeadtime,by=(1/6000))
eadtimetrim<-head(eadtime,-1)

gcfidonsets<-c() # plug in numbers from Perl script output

gcfidoffsets<-c() # plug in numbers from Perl script output

gcfidpeaktimes<-c()
for (i in 1:length(gcfidonsets)) {
    gcfidpeaktimes[i]<-mean(c(gcfidonsets[i],gcfidoffsets[i]))
    }

gcfidpeaknames<-c("") # give names for your peaks

eadcolors<-rainbow(40)
gcnameoffsets<-rep_len(c(700,800,900),length(gcfidpeaknames)) # to stagger peak names in the plot

# Define working window for EAD
eadstartmin<- # specify when in the GcEAD file the sample was injected
eadfinishmin<-eadstartmin+20 # here we have a 20 minute window of interest
eadstartmin<-eadstartmin-0.05 # to correct for slight column length diffs
eadfinishmin<-eadfinishmin-0.05 # to correct for slight column length diffs
eadstart<-eadstartmin*6000
eadfinish<-eadfinishmin*6000
eadonsets<-gcfidonsets+eadstartmin
eadoffsets<-gcfidoffsets+eadstartmin

eadpeaks<-c()
for (i in 1:length(eadonsets)) {
    eadpeaks[i]<-mean(c(eadonsets[i],eadoffsets[i]))
    }

# Take First Forward Difference of EAD data
eaddiff<-diff(eadamp,lag=1,differences=1)
eaddiffsubset<-eaddiff[eadstart:eadfinish]

# Subset EAD data into 5 windows
eaddiffsubset1<-eaddiffsubset[1:24000]
eaddiffsubset2<-eaddiffsubset[24001:48000]
eaddiffsubset3<-eaddiffsubset[48001:72000]
eaddiffsubset4<-eaddiffsubset[72001:96000]
eaddiffsubset5<-eaddiffsubset[96001:120000]

# Calculate 3 sigma (SD) above and below the average for each of 5 windows
threesigmabelow1<-mean(eaddiffsubset1)-3*sd(eaddiffsubset1)
threesigmabelow2<-mean(eaddiffsubset2)-3*sd(eaddiffsubset2)
threesigmabelow3<-mean(eaddiffsubset3)-3*sd(eaddiffsubset3)
threesigmabelow4<-mean(eaddiffsubset4)-3*sd(eaddiffsubset4)
threesigmabelow5<-mean(eaddiffsubset5)-3*sd(eaddiffsubset5)
threesigmaabove1<-mean(eaddiffsubset1)+3*sd(eaddiffsubset1)
threesigmaabove2<-mean(eaddiffsubset2)+3*sd(eaddiffsubset2)
threesigmaabove3<-mean(eaddiffsubset3)+3*sd(eaddiffsubset3)
threesigmaabove4<-mean(eaddiffsubset4)+3*sd(eaddiffsubset4)
threesigmaabove5<-mean(eaddiffsubset5)+3*sd(eaddiffsubset5)

# Here, look at the output of gcfidonsets and determine which peaks fall into the five windows (0-4min, 4-8min, 8-12min, 12-16min, 16-20min here)
#gcfidonsets
# [1]  1.737493  3.350826  4.049993  4.209993  4.525826  5.013326  5.575826
# [8]  5.834993  6.132493  7.070826  7.393326  7.700826  8.095826  8.584993
#[15]  9.246659  9.713326 10.173326 10.388326 10.535826 11.055826 11.276659
#[22] 11.620826 12.594159 13.189159 13.607493 14.024993 14.167493 15.139993
#[29] 16.390826 16.639993 17.094993

# 1:2 in 0-4min
# 3:12 in 4-8min
# 13:22 in 8-12min
# 23:28 in 12-16min
# 29:31 in 16-20min

# Define vectors for output of analyses
eadhitfreq<-c()
eadhitcount<-c()
eadhitsecs<-c()

# First window: determine if a data point is outside the +/- 3 sigma threshold
diffhits1<-c()
for (i in 1:23999) {
  myval<-eaddiffsubset1[i]
   if (myval <= threesigmabelow1 | myval >= threesigmaabove1) {
      timestamp<-(i/6000)+(eadstartmin*1) # 1 should be replaced with which window this is
      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      diffhits1<-append(diffhits1,timestamp)
  }
}

# First window: calculate the number of significant hits/second inside each GC-FID peak within the window
for (i in 1:2) {
    temphits<-length(diffhits1[between(diffhits1,eadonsets[i],eadoffsets[i])])
    tempsec<-60*(eadoffsets[i]-eadonsets[i])
    tempratio<-temphits/tempsec
    print(temphits)
    print(tempsec)
    print(tempratio)
    eadhitcount[i]<-temphits
    eadhitsecs[i]<-tempsec
    eadhitfreq[i]<-tempratio
    }

# Second window: determine if a data point is outside the +/- 3 sigma threshold
diffhits2<-c()
for (i in 1:23999) {
  myval<-eaddiffsubset2[i]
   if (myval <= threesigmabelow2 | myval >= threesigmaabove2) {
      timestamp<-(i/6000)+(eadstartmin+4*1) # 1 should be replaced with which window this is minus 1
      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      diffhits2<-append(diffhits2,timestamp)
  }
}

# Second window: calculate the number of significant hits/second inside each GC-FID peak within the window
for (i in 3:12) {
    temphits<-length(diffhits2[between(diffhits2,eadonsets[i],eadoffsets[i])])
    tempsec<-60*(eadoffsets[i]-eadonsets[i])
    tempratio<-temphits/tempsec
    print(temphits)
    print(tempsec)
    print(tempratio)
    eadhitcount[i]<-temphits
    eadhitsecs[i]<-tempsec
    eadhitfreq[i]<-tempratio
   }

# Third window: determine if a data point is outside the +/- 3 sigma threshold
diffhits3<-c()
for (i in 1:23999) {
  myval<-eaddiffsubset3[i]
   if (myval <= threesigmabelow3 | myval >= threesigmaabove3) {
      timestamp<-(i/6000)+(eadstartmin+4*2) # 1 should be replaced with which window this is minus 1
#      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      diffhits3<-append(diffhits3,timestamp)
  }
}

# Third window: calculate the number of significant hits/second inside each GC-FID peak within the window
for (i in 13:22) {
    temphits<-length(diffhits3[between(diffhits3,eadonsets[i],eadoffsets[i])])
    tempsec<-60*(eadoffsets[i]-eadonsets[i])
    tempratio<-temphits/tempsec
    print(temphits)
    print(tempsec)
    print(tempratio)
    eadhitcount[i]<-temphits
    eadhitsecs[i]<-tempsec
    eadhitfreq[i]<-tempratio
    }

# Fourth window: determine if a data point is outside the +/- 3 sigma threshold
diffhits4<-c()
for (i in 1:23999) {
  myval<-eaddiffsubset4[i]
   if (myval <= threesigmabelow4 | myval >= threesigmaabove4) {
      timestamp<-(i/6000)+(eadstartmin+4*3) # 1 should be replaced with which window this is minus 1
#      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      diffhits4<-append(diffhits4,timestamp)
  }
}

# Fourth window: calculate the number of significant hits/second inside each GC-FID peak within the window
for (i in 23:28) {
    temphits<-length(diffhits4[between(diffhits4,eadonsets[i],eadoffsets[i])])
    tempsec<-60*(eadoffsets[i]-eadonsets[i])
    tempratio<-temphits/tempsec
    print(temphits)
    print(tempsec)
    print(tempratio)
    eadhitcount[i]<-temphits
    eadhitsecs[i]<-tempsec
    eadhitfreq[i]<-tempratio
    }

# Fifth window: determine if a data point is outside the +/- 3 sigma threshold
diffhits5<-c()
for (i in 1:23999) {
  myval<-eaddiffsubset5[i]
   if (myval <= threesigmabelow5 | myval >= threesigmaabove5) {
      timestamp<-(i/6000)+(eadstartmin+4*4) # 1 should be replaced with which window this is minus 1
#      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      diffhits5<-append(diffhits5,timestamp)
  }
}

# Fifth window: calculate the number of significant hits/second inside each GC-FID peak within the window
for (i in 29:31) {
    temphits<-length(diffhits5[between(diffhits5,eadonsets[i],eadoffsets[i])])
    tempsec<-60*(eadoffsets[i]-eadonsets[i])
    tempratio<-temphits/tempsec
    print(temphits)
    print(tempsec)
    print(tempratio)
    eadhitcount[i]<-temphits
    eadhitsecs[i]<-tempsec
    eadhitfreq[i]<-tempratio
    }

# Give each peak a color corresponding to whether its hits/second ratio is above that of the solvent (first) peak
eadhitcols<-c()
for (i in 1:length(eadhitfreq)) {
    if (eadhitfreq[i] > eadhitfreq[1]) {
       eadhitcols[i]<-"red"
       } else {
       eadhitcols[i]<-"black"
       }
       }

# Round significant hit frequencies to 3 decimal places for display
eadhitfreq<-round(eadhitfreq,digits=3)

# Let's plot the plots from the start then - make sure to choose a y-axis limit for the EAD data that is appropriate for your insect
png(file="Insect-10Cmin-00to04min-date.png",width=1500,height=750,units="mm",res=300,pointsize=36,bg="transparent")
par(mfrow=c(2,1))
plot(fidtime,fidamp,type="l",xlim=c(0,4),main="GC-FID Signal",xlab="Time (minutes)",ylab="FID amplitude (pA)",ylim=c(0,1000))
for (i in 1:length(eadonsets)) {
    abline(v=gcfidonsets[i],col=eadcolors[i],lwd=2)
    abline(v=gcfidoffsets[i],col=eadcolors[i],lwd=2)
    text(gcfidpeaktimes[i],gcnameoffsets[i],gcfidpeaknames[i],col=eadhitcols[i])
    }
plot(eadtimetrim,eaddiff,type="l",xlim=c(eadstartmin+(4*0),eadstartmin+(4*1)),main="First forward difference of antennal response",xlab="Time (minutes) (adjusted for GC-FID injection time)",ylab="First forward difference",ylim=c(-0.07,0.07))
abline(h=threesigmabelow1,col="black",lwd=5)
abline(h=threesigmaabove1,col="black",lwd=5)
for (i in 1:length(eadonsets)) {
    abline(v=eadonsets[i],col=eadcolors[i],lwd=2)
    abline(v=eadoffsets[i],col=eadcolors[i],lwd=2)
    text(eadpeaks[i],-0.07,eadhitfreq[i],col=eadhitcols[i])
}
dev.off()

png(file="Insect-10Cmin-04to08min-date.png",width=1500,height=750,units="mm",res=300,pointsize=36,bg="transparent")
par(mfrow=c(2,1))
plot(fidtime,fidamp,type="l",xlim=c(4,8),main="GC-FID Signal",xlab="Time (minutes)",ylab="FID amplitude (pA)",ylim=c(0,1000))
for (i in 1:length(eadonsets)) {
    abline(v=gcfidonsets[i],col=eadcolors[i],lwd=2)
    abline(v=gcfidoffsets[i],col=eadcolors[i],lwd=2)
    text(gcfidpeaktimes[i],gcnameoffsets[i],gcfidpeaknames[i],col=eadhitcols[i])
    }
plot(eadtimetrim,eaddiff,type="l",xlim=c(eadstartmin+(4*1),eadstartmin+(4*2)),main="First forward difference of antennal response",xlab="Time (minutes) (adjusted for GC-FID injection time)",ylab="First forward difference",ylim=c(-0.07,0.07))
abline(h=threesigmabelow2,col="black",lwd=5)
abline(h=threesigmaabove2,col="black",lwd=5)
for (i in 1:length(eadonsets)) {
    abline(v=eadonsets[i],col=eadcolors[i],lwd=2)
    abline(v=eadoffsets[i],col=eadcolors[i],lwd=2)
    text(eadpeaks[i],-0.07,eadhitfreq[i],col=eadhitcols[i])
}
dev.off()

png(file="Insect-10Cmin-08to12min-date.png",width=1500,height=750,units="mm",res=300,pointsize=36,bg="transparent")
par(mfrow=c(2,1))
plot(fidtime,fidamp,type="l",xlim=c(8,12),main="GC-FID Signal",xlab="Time (minutes)",ylab="FID amplitude (pA)",ylim=c(0,1000))
for (i in 1:length(eadonsets)) {
    abline(v=gcfidonsets[i],col=eadcolors[i],lwd=2)
    abline(v=gcfidoffsets[i],col=eadcolors[i],lwd=2)
    text(gcfidpeaktimes[i],gcnameoffsets[i],gcfidpeaknames[i],col=eadhitcols[i])
    }
plot(eadtimetrim,eaddiff,type="l",xlim=c(eadstartmin+(4*2),eadstartmin+(4*3)),main="First forward difference of antennal response",xlab="Time (minutes) (adjusted for GC-FID injection time)",ylab="First forward difference",ylim=c(-0.07,0.07))
abline(h=threesigmabelow3,col="black",lwd=5)
abline(h=threesigmaabove3,col="black",lwd=5)
for (i in 1:length(eadonsets)) {
    abline(v=eadonsets[i],col=eadcolors[i],lwd=2)
    abline(v=eadoffsets[i],col=eadcolors[i],lwd=2)
    text(eadpeaks[i],-0.07,eadhitfreq[i],col=eadhitcols[i])
}
dev.off()

png(file="Insect-10Cmin-12to16min-date.png",width=1500,height=750,units="mm",res=300,pointsize=36,bg="transparent")
par(mfrow=c(2,1))
plot(fidtime,fidamp,type="l",xlim=c(12,16),main="GC-FID Signal",xlab="Time (minutes)",ylab="FID amplitude (pA)",ylim=c(0,1000))
for (i in 1:length(eadonsets)) {
    abline(v=gcfidonsets[i],col=eadcolors[i],lwd=2)
    abline(v=gcfidoffsets[i],col=eadcolors[i],lwd=2)
    text(gcfidpeaktimes[i],gcnameoffsets[i],gcfidpeaknames[i],col=eadhitcols[i])
    }
plot(eadtimetrim,eaddiff,type="l",xlim=c(eadstartmin+(4*3),eadstartmin+(4*4)),main="First forward difference of antennal response",xlab="Time (minutes) (adjusted for GC-FID injection time)",ylab="First forward difference",ylim=c(-0.07,0.07))
abline(h=threesigmabelow4,col="black",lwd=5)
abline(h=threesigmaabove4,col="black",lwd=5)
for (i in 1:length(eadonsets)) {
    abline(v=eadonsets[i],col=eadcolors[i],lwd=2)
    abline(v=eadoffsets[i],col=eadcolors[i],lwd=2)
    text(eadpeaks[i],-0.07,eadhitfreq[i],col=eadhitcols[i])
}
dev.off()

png(file="Insect-10Cmin-16to20min-date.png",width=1500,height=750,units="mm",res=300,pointsize=36,bg="transparent")
par(mfrow=c(2,1))
plot(fidtime,fidamp,type="l",xlim=c(16,20),main="GC-FID Signal",xlab="Time (minutes)",ylab="FID amplitude (pA)",ylim=c(0,1000))
for (i in 1:length(eadonsets)) {
    abline(v=gcfidonsets[i],col=eadcolors[i],lwd=2)
    abline(v=gcfidoffsets[i],col=eadcolors[i],lwd=2)
    text(gcfidpeaktimes[i],gcnameoffsets[i],gcfidpeaknames[i],col=eadhitcols[i])
    }
plot(eadtimetrim,eaddiff,type="l",xlim=c(eadstartmin+(4*4),eadstartmin+(4*5)),main="First forward difference of antennal response",xlab="Time (minutes) (adjusted for GC-FID injection time)",ylab="First forward difference",ylim=c(-0.07,0.07))
abline(h=threesigmabelow5,col="black",lwd=5)
abline(h=threesigmaabove5,col="black",lwd=5)
for (i in 1:length(eadonsets)) {
    abline(v=eadonsets[i],col=eadcolors[i],lwd=2)
    abline(v=eadoffsets[i],col=eadcolors[i],lwd=2)
    text(eadpeaks[i],-0.07,eadhitfreq[i],col=eadhitcols[i])
}
dev.off()

# Write your data to a CSV file
insectout<-cbind(round(gcfidpeaktimes,digits=3),gcfidpeaknames,eadhitfreq,eadhitcols)
write.csv(insectout,file="insect-10cmin-hits-date.csv")
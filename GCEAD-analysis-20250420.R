# load necessary library for 'between' function later
library(dplyr)

# load moth GcEAD software CSV and GC-FID CSV traces
moth1<-read.csv("20240311-MsextaF-EADrun1mixv2-1.ead.csv")
fid<-read.csv("GCFID-bigaverage-10Cmin-2runs.csv")

# define time and amplitude of FID signal
fidtime<-fid[,1]
fidamp<-fid[,2]

# pasted in from the script "gcfid-peaksabove.pl"
gcfidonsets<-c(1.74499375,2.70166041666667,3.43999375,4.15832708333333,4.33582708333333,4.65249375,5.10332708333333,5.64832708333333,5.90249375,6.19666041666667,7.11499375,7.43332708333333,7.73416041666667,8.13166041666667,8.61416041666667,9.27249375,9.73332708333333,10.1891604166667,10.4091604166667,10.5516604166667,11.0691604166667,11.29249375,11.64999375,12.6083270833333,13.19749375,13.61499375,14.0308270833333,14.1783270833333,15.14749375,15.4741604166667,16.39499375,16.64749375,17.10499375)

gcfidoffsets<-c(2.69999375,2.71082708333333,3.77332708333333,4.31749375,4.48582708333333,4.78832708333333,5.20916041666667,5.73832708333333,6.00666041666667,6.31499375,7.22499375,7.53082708333333,7.87166041666667,8.24582708333333,8.71332708333333,9.37999375,9.84082708333333,10.30249375,10.5483270833333,10.6483270833333,11.2691604166667,11.4883270833333,11.80499375,12.7158270833333,13.2758270833333,13.7308270833333,14.1358270833333,14.30749375,15.3516604166667,15.47749375,16.5033270833333,16.7641604166667,17.2141604166667)

# get average time between an onset and offset for plotting purposes
gcfidpeaktimes<-c()
for (i in 1:length(gcfidonsets)) {
    gcfidpeaktimes[i]<-mean(c(gcfidonsets[i],gcfidoffsets[i]))
    }

# define peak names for plotting purposes
gcfidpeaknames<-c("Hexane","Hexane","EthylLactate","unknown","Z3Hexen1Ol","IsoamylAcetate","EthylValerate","APinene","Camphene","Benzaldehyde","Z3HexenylAcetate","Limonene","Phenylacetaldehyde","Acetophenone","MethylBenzoate","Nopinone","BenzylAcetate","1PhenylethylAcetate","Verbenone","unknown","Geraniol","Cinnamaldehyde","Indole","3PhenylpropylAcetate","Caryophyllene","AHumulene","Phenethyl2Methylbutyrate","Methylisoeugenol","CaryophylleneOxide*","unknown","Farnesol1","Farnesol2","BenzylBenzoate")

# define amplitude of EAD signal
m1eadamp<-moth1[,2]

# median filter the EAD amplitude signal to remove electrical noise
m1eadmedfilt<-runmed(m1eadamp,17) # we chose 17 as it's the width in samples of the electrical noise spikes we've seen, this can be adjusted to suit your data

# define time of EAD signal since not provided in the file - n.b. GcEAD samples at 100Hz = 100 times/second so need to multiply the start time by 6000 
m1maxeadtime<-(length(m1eadamp)-1)/6000
m1eadtime<-seq(0,m1maxeadtime,by=(1/6000))
m1eadtimetrim<-head(m1eadtime,-1)

# define EAD start time, i.e. when the GC injection happened according to the GcEAD software
m1eadstartmin<-48.87
m1eadfinishmin<-m1eadstartmin+20 # this run has no peaks after ca. 18 minutes so cut off at 20 minutes
m1eadstartmin<-m1eadstartmin-0.05 # to correct for slight delay between softwares - needs to be checked each time! 
m1eadfinishmin<-m1eadfinishmin-0.05 # to correct for slight delay between softwares - needs to be checked each time!
m1eadstart<-m1eadstartmin*6000
m1eadfinish<-m1eadfinishmin*6000
# define when the GC-FID peaks are in EAD time
m1eadonsets<-gcfidonsets+m1eadstartmin
m1eadoffsets<-gcfidoffsets+m1eadstartmin

# get average time of EAD for peak onset/offset for plotting purposes
m1eadpeaks<-c()
for (i in 1:length(m1eadonsets)) {
    m1eadpeaks[i]<-mean(c(m1eadonsets[i],m1eadoffsets[i]))
    }

# take First Forward Difference (FFD) of median-filtered EAD amplitude signal to remove baseline drift
m1eaddiff<-diff(m1eadmedfilt,lag=1,differences=1)
m1eaddiffsubset<-m1eaddiff[m1eadstart:m1eadfinish]

# divide the FFD into 5 non-overlapping windows of 4 minutes (= 24,000 samples) each
m1eaddiffsubset1<-m1eaddiffsubset[1:24000]
m1eaddiffsubset2<-m1eaddiffsubset[24001:48000]
m1eaddiffsubset3<-m1eaddiffsubset[48001:72000]
m1eaddiffsubset4<-m1eaddiffsubset[72001:96000]
m1eaddiffsubset5<-m1eaddiffsubset[96001:120000]

# calculate 3*sigma above and below the mean of each window
m13sigmabelow1<-mean(m1eaddiffsubset1)-3*sd(m1eaddiffsubset1)
m13sigmabelow2<-mean(m1eaddiffsubset2)-3*sd(m1eaddiffsubset2)
m13sigmabelow3<-mean(m1eaddiffsubset3)-3*sd(m1eaddiffsubset3)
m13sigmabelow4<-mean(m1eaddiffsubset4)-3*sd(m1eaddiffsubset4)
m13sigmabelow5<-mean(m1eaddiffsubset5)-3*sd(m1eaddiffsubset5)
m13sigmaabove1<-mean(m1eaddiffsubset1)+3*sd(m1eaddiffsubset1)
m13sigmaabove2<-mean(m1eaddiffsubset2)+3*sd(m1eaddiffsubset2)
m13sigmaabove3<-mean(m1eaddiffsubset3)+3*sd(m1eaddiffsubset3)
m13sigmaabove4<-mean(m1eaddiffsubset4)+3*sd(m1eaddiffsubset4)
m13sigmaabove5<-mean(m1eaddiffsubset5)+3*sd(m1eaddiffsubset5)

# define placeholder variables for the loops
m1eadhitfreq<-c()
m1eadhitcount<-c()
m1eadhitsecs<-c()

# for window 1, determine which FFD values fall outside the 3*sigma threshold
m1diffhits1<-c()
for (i in 1:23999) {
  myval<-m1eaddiffsubset1[i]
   if (myval <= m13sigmabelow1 | myval >= m13sigmaabove1) {
      timestamp<-(i/6000)+(m1eadstartmin*1) 
      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      m1diffhits1<-append(m1diffhits1,timestamp)
  }
}

# repeat for windows 2-5 with their respective GC-FID peaks
m1diffhits2<-c()
for (i in 1:23999) {
  myval<-m1eaddiffsubset2[i]
   if (myval <= m13sigmabelow2 | myval >= m13sigmaabove2) {
      timestamp<-(i/6000)+(m1eadstartmin+4*1) # 1 should be replaced with which window this is minus 1
      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      m1diffhits2<-append(m1diffhits2,timestamp)
  }
}

m1diffhits3<-c()
for (i in 1:23999) {
  myval<-m1eaddiffsubset3[i]
   if (myval <= m13sigmabelow3 | myval >= m13sigmaabove3) {
      timestamp<-(i/6000)+(m1eadstartmin+4*2) # 2 should be replaced with which window this is minus 1
#      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      m1diffhits3<-append(m1diffhits3,timestamp)
  }
}

fm1diffhits4<-c()
for (i in 1:23999) {
  myval<-m1eaddiffsubset4[i]
   if (myval <= m13sigmabelow4 | myval >= m13sigmaabove4) {
      timestamp<-(i/6000)+(m1eadstartmin+4*3) # 3 should be replaced with which window this is minus 1
#      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      m1diffhits4<-append(m1diffhits4,timestamp)
  }
}

m1diffhits5<-c()
for (i in 1:23999) {
  myval<-m1eaddiffsubset5[i]
   if (myval <= m13sigmabelow5 | myval >= m13sigmaabove5) {
      timestamp<-(i/6000)+(m1eadstartmin+4*4) # 4 should be replaced with which window this is minus 1
#      print(paste("found a hit which was",myval,"at timestamp",timestamp,sep=" "))
      m1diffhits5<-append(m1diffhits5,timestamp)
  }
}

# for window 1, which contains the first 3 GC-FID peaks, calculate the number of significant FFD values within those peak windows (this uses dplyr::between), the length of each peak in seconds, and the ratio of the two 

# calculate the background hit frequency using the above diffhits arrays
m1diffhits1bg<-c()
counter<-0
for (i in 1:length(m1diffhits1)) {
    print(paste("searching",i,sep=" "))
    counter<-0
    for (j in 1:length(m1eadonsets)) {
        # did we see this diffhits data point during an FID window?
        if(m1diffhits1[i] %in% m1diffhits1[between(m1diffhits1,m1eadonsets[j],m1eadoffsets[j])]) {
            counter<-counter+1
            print(paste("found",m1diffhits1[i],"in FID peak window",j,sep=" "))
        }
    }
    if (counter == 0) {
        # no, we didn't see this diffhits data point in an FID window so it is part of the background rate, so keep it!
        print(paste("didn't see",m1diffhits1[i],"in any FID peak window",sep=" "))
        m1diffhits1bg<-c(m1diffhits1bg,m1diffhits1[i])
    }
}

# repeat for windows 2-5
m1diffhits2bg<-c()
counter<-0
for (i in 1:length(m1diffhits2)) {
    print(paste("searching",i,sep=" "))
    counter<-0
    for (j in 1:length(m1eadonsets)) {
        if(m1diffhits2[i] %in% m1diffhits2[between(m1diffhits2,m1eadonsets[j],m1eadoffsets[j])]) {
            counter<-counter+1
            print(paste("found",m1diffhits2[i],"in FID peak window",j,sep=" "))
        }
    }
    if (counter == 0) {
        print(paste("didn't see",m1diffhits2[i],"in any FID peak window",sep=" "))
        m1diffhits2bg<-c(m1diffhits2bg,m1diffhits2[i])
    }
}

m1diffhits3bg<-c()
counter<-0
for (i in 1:length(m1diffhits3)) {
    print(paste("searching",i,sep=" "))
    counter<-0
    for (j in 1:length(m1eadonsets)) {
        if(m1diffhits3[i] %in% m1diffhits3[between(m1diffhits3,m1eadonsets[j],m1eadoffsets[j])]) {
            counter<-counter+1
            print(paste("found",m1diffhits3[i],"in FID peak window",j,sep=" "))
        }
    }
    if (counter == 0) {
        print(paste("didn't see",m1diffhits3[i],"in any FID peak window",sep=" "))
        m1diffhits3bg<-c(m1diffhits3bg,m1diffhits3[i])
    }
}

m1diffhits4bg<-c()
counter<-0
for (i in 1:length(m1diffhits4)) {
    print(paste("searching",i,sep=" "))
    counter<-0
    for (j in 1:length(m1eadonsets)) {
        if(m1diffhits4[i] %in% m1diffhits4[between(m1diffhits4,m1eadonsets[j],m1eadoffsets[j])]) {
            counter<-counter+1
            print(paste("found",m1diffhits4[i],"in FID peak window",j,sep=" "))
        }
    }
    if (counter == 0) {
        print(paste("didn't see",m1diffhits4[i],"in any FID peak window",sep=" "))
        m1diffhits4bg<-c(m1diffhits4bg,m1diffhits4[i])
    }
}

m1diffhits5bg<-c()
counter<-0
for (i in 1:length(m1diffhits5)) {
    print(paste("searching",i,sep=" "))
    counter<-0
    for (j in 1:length(m1eadonsets)) {
        if(m1diffhits5[i] %in% m1diffhits5[between(m1diffhits5,m1eadonsets[j],m1eadoffsets[j])]) {
            counter<-counter+1
            print(paste("found",m1diffhits5[i],"in FID peak window",j,sep=" "))
        }
    }
    if (counter == 0) {
        print(paste("didn't see",m1diffhits5[i],"in any FID peak window",sep=" "))
        m1diffhits5bg<-c(m1diffhits5bg,m1diffhits5[i])
    }
}

m3window1length<-4*60 # the window is 4 minutes long * 60 seconds per minute
m3window2length<-4*60
m3window3length<-4*60
m3window4length<-4*60
m3window5length<-4*60

m3sig1bg<-length(m3diffhits1bg)/m3window1length # what is the background significant spike rate for window 1?
m3sig2bg<-length(m3diffhits2bg)/m3window2length
m3sig3bg<-length(m3diffhits3bg)/m3window3length
m3sig4bg<-length(m3diffhits4bg)/m3window4length
m3sig5bg<-length(m3diffhits5bg)/m3window5length

# for each window's set of FID peaks, calculate the significant spike frequency in spikes/second
for (i in 1:3) {
    temphits<-length(m3diffhits1[between(m3diffhits1,m3eadonsets[i],m3eadoffsets[i])]) # this is how many spikes were in the FID peak
    tempsec<-60*(m3eadoffsets[i]-m3eadonsets[i]) # this is the length of the FID peak
    tempratio<-temphits/tempsec # this is the spikes/second ratio
    m3eadhitcount[i]<-temphits
    m3eadhitsecs[i]<-tempsec
    m3eadhitfreq[i]<-tempratio
    if(m3eadhitfreq[i] > m3sig1bg) { # here we calculate whether this ratio is above the background for this window
       m3eadhitcols[i]<-"red" # if yes, color it red
       } else {
       m3eadhitcols[i]<-"black" # if no, black
       }
}
# repeat for windows 2-5
for (i in 4:13) {
    temphits<-length(m3diffhits2[between(m3diffhits2,m3eadonsets[i],m3eadoffsets[i])])
    tempsec<-60*(m3eadoffsets[i]-m3eadonsets[i])
    tempratio<-temphits/tempsec
    m3eadhitcount[i]<-temphits
    m3eadhitsecs[i]<-tempsec
    m3eadhitfreq[i]<-tempratio
    if(m3eadhitfreq[i] > m3sig2bg) {
       m3eadhitcols[i]<-"red"
       } else {
       m3eadhitcols[i]<-"black"
       }
   }
for (i in 14:23) {
    temphits<-length(m3diffhits3[between(m3diffhits3,m3eadonsets[i],m3eadoffsets[i])])
    tempsec<-60*(m3eadoffsets[i]-m3eadonsets[i])
    tempratio<-temphits/tempsec
    m3eadhitcount[i]<-temphits
    m3eadhitsecs[i]<-tempsec
    m3eadhitfreq[i]<-tempratio
    if(m3eadhitfreq[i] > m3sig3bg) {
       m3eadhitcols[i]<-"red"
       } else {
       m3eadhitcols[i]<-"black"
       }
    }
for (i in 24:30) {
    temphits<-length(m3diffhits4[between(m3diffhits4,m3eadonsets[i],m3eadoffsets[i])])
    tempsec<-60*(m3eadoffsets[i]-m3eadonsets[i])
    tempratio<-temphits/tempsec
    m3eadhitcount[i]<-temphits
    m3eadhitsecs[i]<-tempsec
    m3eadhitfreq[i]<-tempratio
    if(m3eadhitfreq[i] > m3sig4bg) {
       m3eadhitcols[i]<-"red"
       } else {
       m3eadhitcols[i]<-"black"
       }
    }
for (i in 31:33) {
    temphits<-length(m3diffhits5[between(m3diffhits5,m3eadonsets[i],m3eadoffsets[i])])
    tempsec<-60*(m3eadoffsets[i]-m3eadonsets[i])
    tempratio<-temphits/tempsec
    m3eadhitcount[i]<-temphits
    m3eadhitsecs[i]<-tempsec
    m3eadhitfreq[i]<-tempratio
    if(m3eadhitfreq[i] > m3sig5bg) {
       m3eadhitcols[i]<-"red"
       } else {
       m3eadhitcols[i]<-"black"
       }
    }

# the variable 'm1eadhitfreq' now contains the results (example output below) indicating how many significant FFD hits/second were during each GC-FID peak
m1eadhitfreq
# [1]  3.3856894  0.0000000  3.9000000  2.7225131  3.0000000  5.5214724
# [7] 10.8661417  2.9629630  0.0000000  0.0000000 10.3030303  7.5213675
#[13]  6.0606061  2.1897810  3.0252101  3.5658915  5.5813953  3.0882353
#[19]  4.6706587  2.4137931  5.8333333  6.6382979  2.2580645  0.9302326
#[25]  5.9574468  2.3021583  0.0000000  0.2580645  2.2040816  5.0000000
#[31]  1.8461538  1.4285714  1.8320611


# round to 3 decimel places for plotting purposes
m1eadhitfreq<-round(m1eadhitfreq,digits=3)

# stagger GC-FID peak names for clarity when plotting
gcnameoffsets<-rep_len(c(700,800,900),33)

# make a nice rainbow spectrum to line peaks up
eadcolors<-rainbow(40)

# write the results to an external CSV file (with GC-FID peak retention time, GC-FID peak name, significant hit frequency during each peak, and whether it was significant (red/black color)
moth1out<-cbind(gcfidpeaktimes,gcfidpeaknames,m1eadhitfreq,m1eadhitcols)
write.csv(moth1out,file="some-useful-filename-date.csv")

# plot comparing EAD with FFD data - here we are plotting Window 2 (4-8 minutes of elapsed GC-FID time)
png(file="some-pithy-filename-date.png",width=1000,height=700,units="mm",res=300,pointsize=48,bg="transparent")
par(mfrow=c(3,1))
plot(fidtime,fidamp,type="l",xlim=c(4,8),main="GC-FID Signal",xlab="Time (minutes)",ylab="FID amplitude (pA)",ylim=c(0,1000))
for (i in 1:length(gcfidonsets)) {
    abline(v=gcfidonsets[i],col=eadcolors[i],lwd=7)
    abline(v=gcfidoffsets[i],col=eadcolors[i],lwd=7)
    text(gcfidpeaktimes[i],gcnameoffsets[i],gcfidpeaknames[i],col="black")
    }
plot(m1eadtime,m1eadamp,type="l",xlim=c(m1eadstartmin+(4*1),m1eadstartmin+(4*2)),main="Moth 1 (female): raw EAD amplitude",xlab="Time (minutes) (adjusted for GC-FID injection time)",ylab="Raw EAD amplitude (mV)",ylim=c(-2.5,2.5))
for (i in 1:length(m1eadonsets)) {
    abline(v=m1eadonsets[i],col=eadcolors[i],lwd=7)
    abline(v=m1eadoffsets[i],col=eadcolors[i],lwd=7)
}
plot(m1eadtimetrim,m1eaddiff,type="l",xlim=c(m1eadstartmin+(4*1),m1eadstartmin+(4*2)),main="Moth 1 (female): first forward difference of median-filtered antennal response",xlab="Time (minutes) (adjusted for GC-FID injection time)",ylab="First forward difference",ylim=c(-0.3,0.3))
abline(h=m13sigmabelow2,col="black",lwd=5)
abline(h=m13sigmaabove2,col="black",lwd=5)
for (i in 1:length(m1eadonsets)) {
    abline(v=m1eadonsets[i],col=eadcolors[i],lwd=7)
    abline(v=m1eadoffsets[i],col=eadcolors[i],lwd=7)
}
dev.off()

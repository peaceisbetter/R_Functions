#calculate linear analysis of COP data in AP/ML directions
#created 3/4/2022 by Jack Manning
#Packages: astsa, seewave, pracma, zoo, rlist
#C:\Users\jackpmanning\OneDrive - Texas A&M University\Documents\Projects\Manuscript JM_ES\Manning DT Data
#Purpose: The code sets the working directory, loops through subfolders, reads in data from text files, 
#creates new data frames, downsamples the data, calculates center of pressure (COP) positions, and calls the 
#function COPAnalysis.

#load in appropriate packeges
library(astsa)
library(seewave)
library(pracma)
library(zoo)
library(rlist)
library(writexl)
#library(sazedR)

#COP Analysis code
COPAnalysis <- function(sampfreq, COPx, COPy) {

  # calculate number of data points
  n <- length(COPx)
  
  # create time vector
  t <- (0:(n-1)) * (1/sampfreq)
  
  # get current date and time for filename
  datetime <- format(Sys.time(), "%Y_%m_%d_%H-%M")
  
  ## Begin analysis of the segments
  
  # Normalizing x and y data.
  COPxnorm <- COPx - mean(COPx)
  COPynorm <- COPy - mean(COPy)
  
  # Finds radius of COP path.
  COPd <- sqrt(COPxnorm^2 + COPynorm^2)
  
  ## Linear analyses
  
  # Root Mean Square (x,y,d)
  rms_x <- sqrt((1/length(COPxnorm)) * sum(COPxnorm^2))
  rms_y <- sqrt((1/length(COPynorm)) * sum(COPynorm^2))
  rms_d <- sqrt((1/length(COPd)) * sum(COPd^2))
  
  # Range of mediolateral and anteroposterior motion
  rangeml <- max(COPxnorm) - min(COPxnorm)
  rangeap <- max(COPynorm) - min(COPynorm)
  
  # Sway Path
  swaypathd <- sum(abs(diff(COPd)))
  swaypathx <- sum(abs(diff(COPxnorm)))
  swaypathy <- sum(abs(diff(COPynorm)))
  
  # This version of the calculation provides the total tangential distance.
  i <- 1:(n-1)
  swaypathtangential <- sum(sqrt((COPxnorm[i+1] - COPxnorm[i])^2 + 
                                   (COPynorm[i+1] - COPynorm[i])^2))
  
  # Area of 95% Confidence Circle
  # (Notes) Finds the area of a circle that includes 95% of radii(COPd).
  # This is a one sided test so the z score has a value of 1.645. It
  # assumes a normal distribution of radii, as does Prieto. However the
  # distribution of radii is not normally distributed. For now use 1.645
  # from the normal pdf. Chi-square pdf needed?
  
  confcir95 <- mean(COPd) + 1.645 * sd(COPd)
  Acir <- pi * confcir95^2
  
  # Area of 95% Confidence Ellipse
  # (Notes) This is from Prieto (1996) and from Sokal and Rohlf
  # (1995) Biometry p589, also cited by Prieto.
  
  sML <- sqrt((1/length(COPxnorm)) * sum(COPxnorm^2))  # Prieto
  sAP <- sqrt((1/length(COPynorm)) * sum(COPynorm^2))  # Prieto
  sAPML <- (1/length(COPxnorm)) * sum(COPxnorm * COPynorm)  # Prieto
  
  prieto <- list()
  
  prieto$f <- 3
  
  # (Notes) Prieto has missed squaring the sum of the first two terms in his
  # equation(16) page 959. This equation is actually from Sokal and Rohlf
  # (1995).
  
  prieto$d <- sqrt((sAP^2 + sML^2)^2 - 4 * (sAP^2 * sML^2 - sAPML^2))
  
  prieto$ellipRA <- sqrt(prieto$f * (sAP^2 + sML^2 + prieto$d))  # Prieto eq 14, (1995)
  prieto$ellipRB <- sqrt(prieto$f * (sAP^2 + sML^2 - prieto$d))  # Prieto eq 15
  prieto$ellipA <- 2 * pi * prieto$f * sqrt(sAP^2 * sML^2 - sAPML^2)  # Prieto eq 18
  
  prieto$lambda <- (sAP^2 + sML^2 + prieto$d) / 2  # from Sokal and Rohlf
  prieto$slope <- prieto$ellipRA / prieto$ellipRB
  prieto$slope <- sAPML / (prieto$lambda - sML^2) # from Sokal and Rohlf
  prieto$angle <- atan(prieto$slope)
  
  # Plotting 95% circle and ellipse over data
  
  conffig <- plot(COPxnorm, COPynorm, type="l", col="green")
  
  theta <- (pi/180) * seq(0, 360)
  
  # Add 95% circle to plot
  x1circ <- confcir95 * cos(theta)
  y1circ <- confcir95 * sin(theta)
  lines(x1circ, y1circ, col="blue")
  
  # Add ellipse to plot
  
  x2 <- prieto$ellipRA * cos(theta)
  y2 <- prieto$ellipRB * sin(theta)
  x2ellip <- x2 * cos(prieto$angle) - y2 * sin(prieto$angle)
  y2ellip <- x2 * sin(prieto$angle) + y2 * cos(prieto$angle)
  lines(x2ellip, y2ellip, col="red")
  
  title(paste(FileName, ": Green = data, Blue =  95% circle, Red = 95% Ellipse",
              xlab = "COPx (mm)",
              ylab = "COPy (mm)"))
  
  # Add ellipse to plot
  x2 <- prieto$ellipRA*cos(theta)
  y2 <- prieto$ellipRB*sin(theta)
  x2ellip <- x2*cos(prieto$angle)-y2*sin(prieto$angle)
  y2ellip <- x2*sin(prieto$angle)+y2*cos(prieto$angle)
  
  plot(x2ellip, y2ellip, col="red", type = "l")
  title(paste(FileName, ': Green = data, Blue = 95% circle, Red = 95% Ellipse'))
  xlab = 'COPx (mm)'
  ylab = 'COPy (mm)'
  
  # saveas(conffig,fullfile(PathName,[FileName '_95%_Confidence_Circle_and_Ellipse']),'jpg')
  
  # Frequency Domain Analyses
  # (Notes) Prieto, et al (1996) used a sinusiodal multitaper method with eight
  # tapers for their spectral analyses. They cite Minimum Bias Multiple Taper
  # Spectral Emission by Riedel and Sidorenko(1995). It was published in the
  # IEEE Transactions on Signal Processing 43(1) 188-195. This implementation
  # follows the equation on page 188, in paragraph 4 of Riedel and Sidorenko.
  
  N <- length(COPd)
  n <- 1:N
  TaperWindow <- 0.5*(1-cos(2*pi*((n-1)/(N-1))))
  WindowedData <- (COPd*TaperWindow)
  
  TransformFFT <- fft(WindowedData) # Get fourier transform of data.
  
  # Check for an even number of data points, else toss last data point.
  if (length(TransformFFT) %% 2) {
    TransformFFT[length(TransformFFT)] <- NULL
  }
  
  Spectrum <- TransformFFT[1:length(TransformFFT)/2]*Conj(TransformFFT[1:length(TransformFFT)/2])
  
  # (Notes) This is why the index starts at 3 (so points 1 and 2 are removed)
  # The big difference is here a Hanning window is used,
  # whereas Prieto used a multitaper method.
  # (BS) To make sense of this comment look ahead to the Median
  # Frequency calculation about 40 lines down. That is where the
  # index starts at 3.
  
  f <- (.5/length(Spectrum))*sampfreq*(1:length(Spectrum))
  
  powfig <- figure()
  #par(powfig,'Position',[20 40 screen(3)/2-20 screen(4)/2-120],'Visible','off')
  # create a new figure window
  powfig <- plot.new()
   
  # set the position and size of the figure window
  win.width <- (dev.size("px")[1] - 20)/2
  win.height <- (dev.size("px")[2] - 120)/2
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4), xpd = TRUE)
  plot.new()
  grconvertX(c(20, 20+win.width), from = "px", to = "npc")
  grconvertY(c(120+win.height, 120), from = "px", to = "npc")
  par(fig = c(0, 1, 0, 1), mar = c(5, 4, 4, 2) + 0.1, new = TRUE)
   
  # # make the figure window invisible
  # dev.control(displaylist="enable")
  # dev.cur(dev.next())
  # dev.control("off")
  

  
  # Plot entire spectrum.
  subplot(2,1,1)
  plot(f,Spectrum)
  title(paste(FileName, ': Power Spectrum'))
  xlabel('Frequency (Hz)')
  ylabel('Amplitude')
  
  # Plot spectrum where most biomechanical data appears.
  # subplot(2,1,2)
  # plot(f,Spectrum)
  # title(paste(FileName, ': Power Spectrum (to 5 Hz)'))
  # xlabel('Frequency (Hz)')
  # ylabel('Arbitrary Units')
  # # Errors sometime occur here. It seems to be dependent on the
  # # length of the segment.
  # PlotMax <- max(Spectrum[f>0 & f<10])
  # axis(c(0, 5, 0, PlotMax))
  
  #stop sending graphical output to pdf
  dev.off(paste(outputs, "\\", FileName, ".pdf", sep=""))
  
  # Save figure
  # saveas(powfig,fullfile(PathName,[FileName '_Power_Spectrum']),'jpg')
  
  # Median Frequency
  # (Notes) Prieto discards the first two points in the spectra
  # and only uses data to 5 Hz. We have frequencies up at 7, so
  # cutoff here is at 10 Hz - to avoid 60 Hz noise
  # when doing power spectral densities. See left

  # freqindex <- 1
  # while (f[freqindex] < 10 && freqindex != length(f)){
  #   freqindex <- freqindex + 1
  # }
  # analysisspectrum <- Spectrum[3:freqindex]
  # analysisfrequency <- f[3:freqindex]
  # 
  # CumSumPower <- cumsum(analysisspectrum)
  # FindMedian <- which(CumSumPower > 0.5 * sum(analysisspectrum))
  # MedianIndex <- min(FindMedian)
  # 
  # medianfreq <- (0.5 / length(Spectrum)) * sampfreq * MedianIndex

  #Frequency dispersion
  # Mu0 <- (1 / length(analysisspectrum)) * sum(analysisspectrum)
  # Mu1 <- (1 / length(analysisspectrum)) * sum(analysisspectrum * analysisfrequency)
  # Mu2 <- (1 / length(analysisspectrum)) * sum(analysisspectrum * analysisfrequency^2)
  # 
  # freqdisp <- sqrt(1 - Mu1^2 / (Mu0 * Mu2))
  
  #Add new variables to output  
  # copoutput <- data.frame(rms_x, rms_y, rms_d, 
  #                         rangeml, rangeap, 
  #                         swaypathd, swaypathx, 
  #                         swaypathy, swaypathtangential)
  copoutput <- data.frame(
    'rms_x' = rms_x,
    'rms_y' = rms_y,
    'rms_d' = rms_d,
    'rangeap' = rangeap,
    'rangeml' = rangeml,
    'swaypathd' = swaypathd,
    'swaypathx' = swaypathx,
    'swaypathy' = swaypathy,
    'swaypathtangential' = swaypathtangential,
    'Acir' = Acir,
    'prieto.ellipA' = prieto$ellipA)
  
    # 'medianfreq' = medianfreq,
    # 'freqdisp' = freqdisp
  
  return(copoutput)
}



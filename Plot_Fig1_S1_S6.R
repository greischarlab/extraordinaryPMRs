################################
# Code to generate schematic & #
#   calculate PMRs from data.  #
################################

options(warn=1)

# Code to generate schematic from Archer et al. 2018 model:

library(deSolve)

# betaOffsetStart returns a vector of the initial iRBC cohort
#   sorted into age compartments according to a Beta distribution
#   where the median age of infection is adjusted from the default
#   (halfway through IED) according to the offset parameter that can
#   take on any value between 0 and 1 (note the function uses the CDF
#   to sort the initial cohort into compartments)
betaOffsetStart = function(compartments, shape, offset, initialI){
  times = seq(0, 1, length.out = round(compartments + 1)) + offset
  correctedTimes = times
  correctedTimes[times > 1] = times[times > 1] - 1
  pbetaValues = pbeta(correctedTimes, shape1 = shape, shape2 = shape)
  correctedpbetaValues = pbetaValues
  correctedpbetaValues[times > 1] = pbetaValues[times > 1] + 1
  iVector = initialI * diff(correctedpbetaValues)
  if(round(sum(iVector)) != round(initialI)){
    warning(paste("betaOffsetStart magnitude error (",
                  round(sum(iVector)) - round(initialI), ", offset = ", 
                  round(offset, digits = 2), ", shape = ", 
                  round(shape, digits = 2), sep = ''), immediate. = T)}
  if(round(length(iVector)) != round(compartments)){
    warning(paste("betaOffsetStart length error, offset = ", 
                  round(offset, digits = 2), ", shape = ", 
                  shape, sep=''), immediate. = T)}
  return(iVector)
}

# function to convert offset value to median parasite age 
# (as fraction through IED)
# multiply by cycle length in hours to get age in hours
offsetToPeak = function(offset){
  res = rep(NA, length(offset))
  for (i in 1:length(offset)){
    res[i] = 1.5 - offset[i]
    if(offset[i]<0.5){res[i] = 0.5 - offset[i]}
  }
  return(res)
}

# function parameterized from Kriek et al. 2003
# for sequestration as a function of hpi
gfx = function(age, p1 = 11.3869/467.6209, 
               p2 = 1,
               p3 = 18.5802, 
               p4 = 0.2242){
  gVal = p1+((p2-p1)/(1+10^(p4*(p3-age))))
  return(gVal)
} # end gfx function

# generate vector of ages through IED in ten minute increments:
ageVals = seq(1/6,48,by=1/6) 
gValues = gfx(ageVals)
yValues = 1-gValues
qValues = c(0, 1 - yValues[2:length(yValues)]/yValues[1:(length(yValues)-1)])

constantPMR.gamma=function(t, y, parms){
  # model from Archer et al. 2018 PNAS (constant PMR)
  cycleLength = parms[1] # in hours
  mu = parms[2] # mortality of non-sequestered iRBCs
  museq = parms[3] # mortality of sequestered iRBCs
  R = parms[4] # parasite multiplication rate (constant)
  
  lambdaN = 6
  lambdaS = 6
  m = round(6*cycleLength)
  
  N = y[1:m]
  S = y[(m+1):(2*m)]
  
  dN = rep(0, m)
  dS = rep(0, m)
  dN[1] = R*(lambdaN * N[m] + lambdaS * S[m]) - (lambdaN + mu) * N[1]
  dS[1] = 0
  for (i in 2:m){
    dN[i] = (1-qValues[(i-1)]) * lambdaN * N[(i-1)] - (lambdaS + mu) * N[i]
    dS[i] = qValues[(i-1)] * lambdaN * N[(i-1)] + lambdaS * S[(i-1)] - (lambdaS + museq) * S[i]
  } # end for loop through n compartments
  
  res=c(dN,dS)
  list(res)
} # end simplified gamma chain model

pfCycleLength = 48 # hours

mu = 0
museq = 0
R = 6
parmValues = c(pfCycleLength, mu, museq, R)

# choose offset value equivalent to initial 
#   median age of 16 hours
startI0All = betaOffsetStart(pfCycleLength*6,100,0.5-16/48,1000)
startI0 = c((1-qValues)*startI0All,qValues*startI0All)

timeValues = seq(0, pfCycleLength*4, by = 1/6) # time vector in hours
archer.output = as.data.frame(lsoda(startI0, times = timeValues,
                              func = constantPMR.gamma, 
                              parms = parmValues))

# calculate circulating iRBC abundance
circ.iRBC = apply(archer.output[ , 2:(pfCycleLength*6 + 1)],1,sum)
# calculate sequestered iRBC abundance
seq.iRBC = apply(archer.output[ , (pfCycleLength*6 + 2):(pfCycleLength*6*2 + 1)],1,sum)

# Code to analyze Wockner et al. 2020 data:

# load package to read in Excel data 
  # from Wockner et al. 2020 J Inf Dis
library("readxl")
# to use Arial fonts in figures:
library("extrafont")
#font_import() # takes awhile only run when needed
#font_import(recursive=FALSE) # takes awhile only run when needed
quartzFonts(avenir=c("Avenir Book", "Avenir Black", "Avenir Book Oblique", "Avenir Black Oblique"),Arial=c("Arial Book", "Arial Black", "Arial Book Oblique", "Arial Black Oblique"))
loadfonts(device = "postscript")

## @knitr pmrfx
# fixed bug with incorrect indexing with short time series that have missing values
pmrfx=function(timeSeries, alpha){
  # function to calculate parasite multiplication rates
  
  # 1. any zeroes should be coded as missing values rather than true zeroes
  # but retain and output the original counts as para.obs column
  obsPara = timeSeries$Asex
  
  # 2. all days encompassed in the experiment, even the days on which samples were not
  # taken while ensuring that the first time point contains a parasite count
  # and the last time point contains a gametocyte count
  allTime = seq(min(timeSeries$day), 
                max(timeSeries$day), 
                by = min(diff(timeSeries$day)))
  
  # 3. build a matrix to store results (and missing values)
  # missing values might be recorded as NAs, or rows could simply be missing
  # both possibilities dealt with using the matrix pmrRes
  pmrRes = as.data.frame(matrix(NA, nrow = length(allTime), ncol = 3))
  colnames(pmrRes) <- c("day", "pmr.obs", "para.obs")

  # 3a. store all time points encompassed in experiment, 
  # even ones on which no samples were recorded
  pmrRes$day = allTime
  # 3b. store all parasite counts recorded
  # if any days were unsampled (i.e., rows are missing from the time series),
  # the parasite counts will be recorded here as NA
  pmrRes$para.obs[which(signif(allTime) %in% signif(timeSeries$day))] <-     
    obsPara[which(signif(timeSeries$day) %in% signif(allTime))]
  
  # 4. calculate parasite multiplication rate (PMR) from observed time series
  alphaIndex = alpha/(pmrRes$day[2]-pmrRes$day[1])
  estPMR = pmrRes$para.obs[(alphaIndex+1):(length(pmrRes$para.obs))] / 
    pmrRes$para.obs[1:(length(pmrRes$para.obs) - alphaIndex)]
  pmrRes$pmr.obs[1:(length(pmrRes[,1]) - alphaIndex)] <- estPMR[1:(length(pmrRes[,1]) - alphaIndex)]
  
  return(pmrRes)
} # end pmrfx function

# need to download supplementary table 2 from Wockner et al.
  # can be found at https://academic.oup.com/jid/article/221/6/963/5611305?login=false#supplementary-data
  # add file to the working directory before running next line:
chiData = read_xlsx("jiz557_suppl_supplementary_table_2.xlsx", range = "R12C1:R1142C7", col_names = T)

chiPMR = data.frame(id = rep(NA, length(unique(chiData$Inoculum.Coh.Subject.ID))), 
                       pmrObs = rep(NA, length(unique(chiData$Inoculum.Coh.Subject.ID))),
                       pmrReg = rep(NA, length(unique(chiData$Inoculum.Coh.Subject.ID))))

chiTimes = seq(min(chiData$Time), max(chiData$Time), by = 0.5)
chiParaCounts = matrix(NA, nrow = length(unique(chiData$Inoculum.Coh.Subject.ID)), 
                       ncol = length(chiTimes))

if(exists("chiPMRThruTime")){rm(chiPMRThruTime)}
for (pindex in 1:length(unique(chiData$Inoculum.Coh.Subject.ID))){
  id = unique(chiData$Inoculum.Coh.Subject.ID)[pindex]
  subData = chiData[which(chiData$Inoculum.Coh.Subject.ID==id),]
  colnames(subData) <- c("day", "Treatment.Day", "log10.av.ppml", "Asex", "Time.centered",
                         "Inoculum.Cohort.ID", "Inoculum.Coh.Subject.ID")
  pmrEstSubData = pmrfx(timeSeries = subData, alpha=2)
  yValues = pmrEstSubData$para.obs
  xValues = pmrEstSubData$day
  
  chiThruTime = data.frame(id = rep(pindex, length(which(is.finite(pmrEstSubData$pmr.obs)))),
                           day = xValues[which(is.finite(pmrEstSubData$pmr.obs))], 
                           obsPMR = pmrEstSubData$pmr.obs[which(is.finite(pmrEstSubData$pmr.obs))])
  if(exists("chiPMRThruTime")){newRes = rbind(chiPMRThruTime, chiThruTime, deparse.level = 0)}
  if(!exists("chiPMRThruTime")){newRes = chiThruTime}
  chiPMRThruTime = newRes
  
  linr = lm(log10(yValues)~xValues)
  pmrReg.chi = 10^(2*linr$coefficients[[2]])
  chiPMR[pindex,] = c(id, max(pmrEstSubData$pmr.obs, na.rm = T), pmrReg.chi)
  chiParaCounts[pindex, which(chiTimes %in% xValues)] = yValues
  }

# Calculate PMRs for MT data

asexLimit = 10 # detection limit based on Simpson et al. 2002 Parasitology

# read in table of days listing when drugs were administered:
drugdays=read.csv("DrugDays.csv", header=T)
# read in table of circulating iRBC counts
asex=read.csv("Asexual.csv",header=F)

mtPMR = data.frame(id = rep(NA, 322), # patient id
                   pmrObs = rep(NA, 322), # max PMR_obs over first seven days
                   pmrReg = rep(NA, 322), # PMR_reg over first seven days
                   toPlot = rep(NA, 322)) # set to T if both PMRs are defined
                                            # and population is not declining
                                            # (which would render regression inappropriate)

mtParaCounts = matrix(NA, nrow = 322, ncol = 7)
mtParaCountsExcl = matrix(NA, nrow = 322, ncol = 7)

if(exists("mtPMRThruTime")){rm(mtPMRThruTime)}
for (patientIndex in 1:322){
  drugDay = NA 
  if(any(!is.na(unlist(drugdays[patientIndex,])))){
    firstDrugDay = min(unlist(drugdays[patientIndex,]), na.rm=T)
    if(firstDrugDay<=7){drugDay = firstDrugDay}
    } # end conditional are any drugdays listed for this patient
  
  if(is.na(drugDay)){
    workingData=as.data.frame(cbind(1:480, asex[, patientIndex]))
    colnames(workingData)<-c("day", "Asex")
    
    questionablePts = which(workingData$Asex<asexLimit)
    indexForReg = c(1:7)
    # exclude any values that coincide with iRBC counts below detection limit
    if(any(indexForReg %in% questionablePts)){
      indexForReg = c(1:7)[-which(c(1:7)%in%questionablePts)]
      }
    # index appropriate iRBC and day values for regression
    paraForLM = workingData$Asex[indexForReg]
    dayForLM = workingData$day[indexForReg]
    
    # calculate PMR_reg
    linr.mt = lm(log10(paraForLM)~dayForLM)
    pmrReg.mt = 10^(2*linr.mt$coefficients[[2]])
    
    # calculate PMRobs from same iRBC counts used in regression
    pmrEstMTData = pmrfx(timeSeries = workingData[indexForReg,], alpha = 2)
    pmrObs.mt = NA
    # if iRBC counts contain missing values, it may not be possible 
      # to calculate any PMRobs values, much less a maximum
    if(any(!is.na(pmrEstMTData$pmr.obs))){
      pmrObs.mt = max(pmrEstMTData$pmr.obs, na.rm = T)
      }
    
    # identify early peaks by checking for early declines (PMR<1)
    declPMR = F
    # set declPMR to T if any PMR falls below 1
    if(any(pmrEstMTData$pmr.obs<1, na.rm = T)){declPMR = T}
    # include in plot if both pmrObs and pmrReg are defined,
      # and PMR is not declining
    toPlot = is.finite(pmrObs.mt) & is.finite(pmrReg.mt) & declPMR==F
    if(toPlot==F){
      mtParaCountsExcl[patientIndex, indexForReg] = paraForLM
      } # end conditional: should data for this patient be excluded?
    if(toPlot){
      mtPMR[patientIndex,] = c(patientIndex, pmrObs.mt, pmrReg.mt, toPlot)
      mtParaCounts[patientIndex, indexForReg] = paraForLM
      # calculate observed PMRs from entire infections to plot later
      # remove `questionable points` below detection threshold
      allWorkingData = workingData
      if(length(questionablePts)>0){
        allWorkingData = workingData[-questionablePts,]}
      pmrAllMTData = pmrfx(timeSeries = allWorkingData, alpha = 2)
      # check whether any PMRs were defined and save those
      if(any(!is.na(pmrAllMTData$pmr.obs))){
        mtThruTime = data.frame(id = rep(patientIndex, length(which(is.finite(pmrAllMTData$pmr.obs)))),
                                day = pmrAllMTData$day[which(is.finite(pmrAllMTData$pmr.obs))], 
                                obsPMR = pmrAllMTData$pmr.obs[which(is.finite(pmrAllMTData$pmr.obs))])
        if(exists("mtPMRThruTime")){newRes = rbind(mtPMRThruTime, mtThruTime, deparse.level = 0)}
        if(!exists("mtPMRThruTime")){newRes = mtThruTime}
        mtPMRThruTime = newRes
        } # end conditional: any PMRs calculable from this patient?
      } # end conditional: should data for this patient be plotted?
    } # end conditional: drugs were not used in first seven days
  } # end patientIndex loop
# R automatically converts numeric values to strings in data.frame
  # because id's are strings
  # switch them back:
chiPMR$pmrObs = as.numeric(chiPMR$pmrObs)
chiPMR$pmrReg = as.numeric(chiPMR$pmrReg)

# order values from CHI data
pmrRegOrder = order(chiPMR$pmrReg)

cexpt = 0.4
cexText = 0.83
lineCol = "gray70"
maxLogPMR = 1.04*log10(max(c(chiPMR$pmrObs,chiPMR$pmrReg),na.rm=T))
max.iRBC = max(c(chiParaCounts,mtParaCounts),na.rm=T)

postscript("Fig1.eps", height = 7.25, width=6.5, paper="special", horizontal=F, family="ArialMT")

par(mfrow = c(3,2), mar = c(1.5,3,0,0), bty='n', oma = c(1.5,0,0,0), xpd = F)

plot(archer.output$time/24, log10(circ.iRBC+seq.iRBC), type='l', col = 'black', 
     xaxt = 'n', yaxt = 'n',
     xlim = c(0, 9), ylim = c(0, 6.5), xaxs = 'i', yaxs = 'i', 
     xlab = "",
     ylab = "", lwd = 1)
axis(1, at = 0:9, labels = c(0,"",2,"",4,"",6,"",8, ""), tcl = -0.1, mgp = c(3,0.1,0))
axis(2, at = log10(c(1,10,100,1e3,1e4,1e5,1e6)), 
     labels = c(1,10,expression(10^2),expression(10^3),expression(10^4),
                expression(10^5),expression(10^6)), las=2, tcl = -0.1, mgp = c(3,0.2,0))
mtext(expression("iRBC abundance"), side = 2, line = 1.6, cex = 0.66)
points(archer.output$time/24, log10(circ.iRBC), type='l', col = lineCol, lwd = 2)
# calculate PMRs
# index appropriate iRBC and day values
indexForReg = which((archer.output$time/24)%in%(1:7))
paraForLM = circ.iRBC[indexForReg]
dayForLM = (archer.output$time/24)[indexForReg]
# calculate PMR_obs
workingData = data.frame(day = dayForLM, Asex = paraForLM)
pmrEstSimData = pmrfx(timeSeries = workingData, alpha = 2)
# calculate PMR_reg
linr.sim = lm(log10(paraForLM)~dayForLM)
points(dayForLM[c(1,3,5,7)], log10(paraForLM[c(1,3,5,7)]),
       type = 'l', col = 'red', lwd = 0.5)
points(dayForLM[c(2,4,6)], log10(paraForLM[c(2,4,6)]),
       type = 'l', col = 'red', lwd = 0.5)
# bold first obs PMR segment
points(dayForLM[c(1,3)], log10(paraForLM[c(1,3)]), type = 'l',
       col = 'red', lwd = 2)
points(dayForLM, coef(linr.sim)[[1]] + dayForLM*coef(linr.sim)[[2]], 
       type = 'l', col = "blue", lty = 2, lwd = 2)
points((archer.output$time/24)[indexForReg], 
       log10(circ.iRBC[indexForReg]), col = lineCol,
       pch = 21, bg = "white", cex = 1.3)
legend(x = "bottomright", 
       legend = c("Total", "Circulating", "Observed (max. bolded)", "From regression"),
       col = c('black', lineCol, "red", "blue"), lty = c(1,1,1,2), 
       lwd = c(1,2,0.5,2), bty='n', y.intersp = 1.2)
text(x=0.05*9, y=0.95*6.5, "(A) Simulated data", adj = 0)

par(mar = c(1.5,2.75,0,0.25))

plot(c(0,10), y = c(0,1), ylim = c(0,maxLogPMR),
     type = 'n', xlab = "", ylab = "", xaxt = 'n', yaxt = 'n', yaxs = 'i')
rect(xleft = -20, xright = 400, ybottom = log10(1), ytop = log10(16), col = "gray80", border=NA)
rect(xleft = -20, xright = 400, ybottom = log10(16), ytop = log10(32), col = "gray90", border=NA)
axis(2,las=2,at = log10(c(1,2,4,8,16,32,64,128,256,512,1024,2048)), 
     labels = c(1,2,4,8,16,32,64,128,256,512,1024,2048),
     tcl = -0.1, mgp = c(3,0.2,0))
mtext("Estimated PMR", side=2, line=1.6, cex = 0.66)
points(x = 0, y = log10(max(pmrEstSimData$pmr.obs, na.rm = T)), 
       pch=21, col = "red", bg = "red", cex = 1.2*cexpt)
text(x = 0.3, y = log10(max(pmrEstSimData$pmr.obs, na.rm = T)),
     bquote("Max. PMR"[italic("obs")]~"="~.(round(max(pmrEstSimData$pmr.obs, na.rm = T)))),
     col = "red", adj = c(0,0.5))
points(x = 0, y = log10(round(10^(2*coef(linr.sim)[[2]]))), 
       pch=23,col = "blue", bg = "white",cex = 1.5*cexpt)
text(x = 0.3, y = log10(round(10^(2*coef(linr.sim)[[2]]))),
     bquote("PMR"[italic("reg")]~"="~.(round(10^(2*coef(linr.sim)[[2]])))),
     col = "blue", adj = c(0,0.5))
abline(h=log10(R), lty=1)
rect(xleft = 0, xright = 3.2, ybottom = log10(4), ytop = log10(8), 
     col = "gray80", border=NA)
text(x = 0.3, y = log10(R), bquote("True PMR ="~.(R)), adj = c(0,0.5))
text(x=0.1, y = 0.97*maxLogPMR, "(B) Inference from\n    simulated data", 
     adj = c(0,1))

par(mar = c(1.5,3,0,0))

plot(x = range(chiTimes), y = range(0:6), type = 'n', 
     xaxt = 'n', yaxt = 'n', xlab = "", ylab = "", 
     xlim = c(0,9), ylim = c(0,log10(max.iRBC)),
     xaxs = 'i', yaxs = 'i')
axis(1, at = 0:9, labels = c(0,"",2,"",4,"",6,"",8, ""), tcl = -0.1, mgp = c(3,0.1,0))
axis(2, tcl = -0.1, mgp = c(3,0.2,0), las = 2, at = c(0:5), labels = c(1, 10, expression(10^2),
                                                                       expression(10^3),
                                                                       expression(10^4),
                                                                       expression(10^5)))
mtext("iRBC abundance", side=2, line=1.6, cex = 0.66)

for (id in c(1:length(chiParaCounts[,1]))){
  subDataParaToPlot = chiParaCounts[id, ]
  subDataToPlot = data.frame(day = chiTimes[!is.na(subDataParaToPlot)], 
                             para = subDataParaToPlot[!is.na(subDataParaToPlot)])
  if(id==which.max(chiPMR$pmrReg)){maxDataToPlot = subDataToPlot}
  points(subDataToPlot$day, log10(subDataToPlot$para), type = 'l', col = lineCol)
  } # end for loop through patient IDs
# overplot trajectory corresponding to max PMR from regression
points(maxDataToPlot$day, log10(maxDataToPlot$para), type = 'l', col = "blue", lwd = 3)
legend(x = 0.2, y = 5, 
       legend = c(expression("Individual trajectory"~phantom("PMR"[italic("reg")])),
                  expression("associated with"~phantom("PMR"[italic("reg")])),
                  expression("max."~"PMR"[italic("reg")]),
                  expression("Individual trajectory")),
       col = c("blue", "white", "white", lineCol), lty = c(1,1,1,1), lwd = c(3,1,1,1), 
       bty = 'n', adj = c(0,0.5), y.intersp = 1.1)
text(x = 0.05*9, y = 0.95*log10(max.iRBC), "(C) CHI data", adj=0)

par(mar = c(1.5,2.75,0,0.25))

plot(log10(as.numeric(chiPMR$pmrObs[pmrRegOrder])),type = 'n',yaxt='n',xaxt='n',
     ylab = "", xlab = "", ylim = c(0,maxLogPMR), yaxs = 'i')
mtext("Estimated PMR", side=2, line=1.6, cex = 0.66)
rect(xleft = -20, xright = 400, ybottom = log10(1), ytop = log10(16), col = "gray80", border=NA)
rect(xleft = -20, xright = 400, ybottom = log10(16), ytop = log10(32), col = "gray90", border=NA)
points(log10(as.numeric(chiPMR$pmrObs[pmrRegOrder])), pch=21, col = "red", bg = "red", cex = 1.2*cexpt)
axis(2,las=2,at = log10(c(1,2,4,8,16,32,64,128,256,512,1024,2048)), 
     labels = c(1,2,4,8,16,32,64,128,256,512,1024,2048),
     tcl = -0.1, mgp = c(3,0.2,0), cex=cexText)
arrows(x0 = 1.035*length(pmrRegOrder),
       x1 = 1.01*length(pmrRegOrder),
       y0 = 0.99*max(log10(as.numeric(chiPMR$pmrReg[pmrRegOrder]))) - 0.05*maxLogPMR,
       y1 = 0.99*max(log10(as.numeric(chiPMR$pmrReg[pmrRegOrder]))), 
       col = "blue", length = 0.05, angle = 20)
points(log10(as.numeric(chiPMR$pmrReg[pmrRegOrder])),pch=23,col = "blue", bg = "white", cex=1.5*cexpt)
text(x = 0.01*length(pmrRegOrder), y = 0.95*maxLogPMR, 
     "(D) Inference from CHI data", adj=0)

par(mar = c(1.5,3,0,0))

plot(x = range(1:7), y = range(0:6), type = 'n', 
     xaxt = 'n', yaxt = 'n', 
     xlim = c(0,9), ylim = c(0,log10(max.iRBC)),
     xaxs = 'i', yaxs = 'i')
axis(1, at = 0:9, labels = c(0,"",2,"",4,"",6,"",8, ""), tcl = -0.1, mgp = c(3,0.1,0))
axis(2, tcl = -0.1, mgp = c(3,0.2,0), las = 2, at = c(0:5), cex = cexText,
     labels = c(1, 10, expression(10^2),
                       expression(10^3),
                       expression(10^4),
                       expression(10^5)))
mtext("iRBC abundance", side=2, line=1.6, cex = 0.66)
mtext("Day of blood-stage infection", side=1, line=0, cex = 0.66, 
      adj = 0.23, outer=T)

for (id in c(1:length(mtParaCounts[,1]))){
  subDataParaToPlot = mtParaCounts[id, ]
  subDataToPlot = data.frame(day = (1:7)[!is.na(subDataParaToPlot)], 
                             para = subDataParaToPlot[!is.na(subDataParaToPlot)])
  if(id==which.max(mtPMR$pmrReg)){maxDataToPlot = subDataToPlot}
  points(subDataToPlot$day, log10(subDataToPlot$para), type = 'l', col = lineCol)
  } # end loop through patient IDs
# overplot trajectory corresponding to max PMR from regression
points(maxDataToPlot$day, log10(maxDataToPlot$para), type = 'l', col = "blue", lwd = 3)

text(x = 0.05*9, y = 0.95*log10(max.iRBC), "(E) MT data", adj=0)

par(mar = c(1.5,2.75,0,0.25))
# order values for MT data 
  # (only consider values for which both PMRobs and PMRreg are defined,
  # and for which there is no evidence of decline)
mtToPlot = mtPMR[!is.na(mtPMR$toPlot),]
regOrdering = order(mtToPlot$pmrReg)

plot(log10(as.numeric(mtToPlot$pmrObs[regOrdering])),type='n',yaxt='n',xaxt='n',
     ylab = "", xlab = "", ylim = c(0,maxLogPMR), yaxs = 'i')
mtext("Individual", side=1, line=0, cex = 0.66, outer=T, adj = 0.79)
mtext("Estimated PMR", side=2, line=1.6, cex = 0.66)
rect(xleft = -20, xright = 400, ybottom = log10(1), ytop = log10(16), col = "gray80", border=NA)
rect(xleft = -20, xright = 400, ybottom = log10(16), ytop = log10(32), col = "gray90", border=NA)
points(log10(as.numeric(mtToPlot$pmrObs[regOrdering])),pch=21, col = "red", bg = "red", cex = 1.2*cexpt)
axis(2,las=2,at = log10(c(1,2,4,8,16,32,64,128,256,512,1024,2048)), 
     labels = c(1,2,4,8,16,32,64,128,256,512,1024,2048),
     tcl = -0.1, mgp = c(3,0.2,0),cex=cexText)
arrows(x0 = 1.035*length(regOrdering), 
       x1 = 1.01*length(regOrdering), 
       y0 = 0.99*max(log10(mtPMR$pmrReg), na.rm=T) - 0.05*maxLogPMR,
       y1 = 0.98*max(log10(mtPMR$pmrReg), na.rm=T), 
       col = "blue", length = 0.05, angle = 20)
points(log10(as.numeric(mtToPlot$pmrReg[regOrdering])),pch=23,col = "blue", bg = "white", cex = 1.5*cexpt)

text(x = 0.01*length(regOrdering), y = 0.95*maxLogPMR, 
     "(F) Inference from MT data", adj=0)

dev.off()

c(max(chiPMR$pmrReg), max(mtToPlot$pmrReg[regOrdering]))

c(max(chiPMR$pmrObs), max(mtToPlot$pmrObs[regOrdering]))

# calculate fraction of extraordinary and borderline PMRreg's from each dataset
length(which(as.numeric(chiPMR$pmrReg[pmrRegOrder])>=32))/length(pmrRegOrder)
length(which(as.numeric(chiPMR$pmrReg[pmrRegOrder])>=16))/length(pmrRegOrder)

length(which(as.numeric(mtToPlot$pmrReg[regOrdering])>=32))/length(regOrdering)
length(which(as.numeric(mtToPlot$pmrReg[regOrdering])>=16))/length(regOrdering)

# calculate max PMR & associated day for each infection
maxDayCHI = data.frame(id = unique(chiPMRThruTime$id),
                       day = rep(NA, length(unique(chiPMRThruTime$id))),
                       obsPMR = rep(NA, length(unique(chiPMRThruTime$id))))
for (pid in chiPMRThruTime$id){
  subData = data.frame(day = chiPMRThruTime$day[chiPMRThruTime$id==pid],
                       obsPMR = chiPMRThruTime$obsPMR[chiPMRThruTime$id==pid])
  maxDayCHI$day[which(maxDayCHI$id==pid)] = subData$day[which.max(subData$obsPMR)]
  maxDayCHI$obsPMR[which(maxDayCHI$id==pid)] = max(subData$obsPMR)
  }

maxDayMT = data.frame(id = unique(mtPMRThruTime$id),
                       day = rep(NA, length(unique(mtPMRThruTime$id))),
                       obsPMR = rep(NA, length(unique(mtPMRThruTime$id))))
for (pid in mtPMRThruTime$id){
  subData = mtPMRThruTime[mtPMRThruTime$id==pid & mtPMRThruTime$day<=5,]
  maxDayMT$day[which(maxDayMT$id==pid)] = subData$day[which.max(subData$obsPMR)]
  maxDayMT$obsPMR[which(maxDayMT$id==pid)] = max(subData$obsPMR)
  }

# plot for supplement: observed PMRs through time
cexpt = 0.6
postscript("FigS6.eps", width = 6.5, height = 7.5, paper = "special", horizontal = F)

par(mfrow = c(3,1), mar = c(1,3,1,0), bty='n', oma = c(1,0,0,0.5), xpd = NA)

plot(c(0.5,7.5), c(0, log10(2048)), type = 'n', xlab = '', ylab = '',
     xaxt='n', yaxt='n', ylim = c(0, log10(2048+1)), xaxs = 'i', yaxs = 'i')
rect(xleft = 0.5, xright = 7.5, ybottom = log10(1), ytop = log10(16+1), col = "gray80", border=NA)
rect(xleft = 0.5, xright = 7.5, ybottom = log10(16+1), ytop = log10(32+1), col = "gray90", border=NA)
axis(1, tcl = -0.1, mgp = c(0,0.1,0), cex = 0.83, cex.axis = 0.83)
axis(2,las=2,at = log10(c(0,1,2,4,8,16,32,64,128,256,512,1024,2048)+1), 
     labels = c(0,1,2,4,8,16,32,64,128,256,512,1024,2048),
     tcl = -0.1, mgp = c(3,0.2,0), cex=0.83, cex.axis = 0.83)
for (pid in chiPMRThruTime$id){
  subData = chiPMRThruTime[chiPMRThruTime$id==pid,]
  points(subData$day, 
         log10(subData$obsPMR+1), type = 'l', lwd = 0.5, col = 'darkred')
  }
points(chiPMRThruTime$day, 
       log10(chiPMRThruTime$obsPMR+1), pch = 4, col = 'darkred', cex = 0.8*cexpt)
points(maxDayCHI$day, log10(maxDayCHI$obsPMR+1), pch = 16, col = 'red', cex = 1.4*cexpt)
text(x = 0.5 + 0.02*7, y = 0.99*log10(2048), "(A) CHI data", adj = 0, cex = 0.83)
legend(x = 0.5 + 0.02*7, y = log10(512+1), legend = c(expression(italic("PMR"["obs"])), expression("Max."~italic("PMR"["obs"]))),
       col = c("darkred", "red"), pch = c(4,16), pt.cex = c(0.8,1.2), y.intersp = 1.2, cex = cexText, bty = 'n')

plot(c(0.5,7.5), c(0, log10(2048)), type = 'n',
     xaxt='n', yaxt='n', ylim = c(0, log10(2048+1)), xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
rect(xleft = 0.5, xright = 7.5, ybottom = log10(1), ytop = log10(16+1), col = "gray80", border=NA)
rect(xleft = 0.5, xright = 7.5, ybottom = log10(16+1), ytop = log10(32+1), col = "gray90", border=NA)
axis(1, tcl = -0.1, mgp = c(0,0.1,0), cex = 0.83, cex.axis = 0.83)
axis(2,las=2,at = log10(c(0,1,2,4,8,16,32,64,128,256,512,1024,2048)+1), 
     labels = c(0,1,2,4,8,16,32,64,128,256,512,1024,2048),
     tcl = -0.1, mgp = c(3,0.2,0), cex=0.83, cex.axis = 0.83)
mtext("Observed PMR", side=2, line=-1, outer=T, cex = 0.66)
for (pid in mtPMRThruTime$id){
  subData = mtPMRThruTime[mtPMRThruTime$id==pid & mtPMRThruTime$day<=5,]
  points(subData$day, 
         log10(subData$obsPMR+1), type = 'l', lwd = 0.5, col = 'darkred')
  }
points(mtPMRThruTime$day[mtPMRThruTime$day<=5], 
       log10(mtPMRThruTime$obsPMR[mtPMRThruTime$day<=5]+1), pch = 4, col = 'darkred', cex = 0.8*cexpt)
points(maxDayMT$day, log10(maxDayMT$obsPMR+1), pch = 16, col = 'red', cex = 1.4*cexpt)
text(x = 0.5 + 0.02*7, y = 0.99*log10(2048+1), "(B) MT data, first 5 days", adj = 0, cex = 0.83)

plot(c(0,max(mtPMRThruTime$day)), c(0, log10(2048)), type = 'n',
     xaxt='n', yaxt='n', ylim = c(0, log10(2048+1)), xaxs = 'i', yaxs = 'i', xlab = '', ylab = '')
rect(xleft = 0, xright = max(mtPMRThruTime$day), ybottom = log10(1), ytop = log10(16+1), col = "gray80", border=NA)
rect(xleft = 0, xright = max(mtPMRThruTime$day), ybottom = log10(16+1), ytop = log10(32+1), col = "gray90", border=NA)
axis(1, tcl = -0.1, mgp = c(0,0.1,0), cex = 0.83, cex.axis = 0.83)
axis(2,las=2,at = log10(c(0,1,2,4,8,16,32,64,128,256,512,1024,2048)+1), 
     labels = c(0,1,2,4,8,16,32,64,128,256,512,1024,2048),
     tcl = -0.1, mgp = c(3,0.2,0), cex=0.83, cex.axis = 0.83)
mtext("Day of blood-stage infection", side=1, line=1, cex = 0.66, adj = 0.5)
points(mtPMRThruTime$day, 
       log10(mtPMRThruTime$obsPMR+1), pch = 4, col = 'darkred', cex = 0.8*cexpt)
points(maxDayMT$day, log10(maxDayMT$obsPMR+1), pch = 16, col = 'red', cex = 1.4*cexpt)
text(x = 0.02*max(mtPMRThruTime$day), y = 0.99*log10(2048+1), "(C) MT data (all)", adj = 0, cex = 0.83)

dev.off()

# plot for supplement: patients excluded from MT data:
postscript("FigS1.eps", height = 4, width=4, paper="special", horizontal=F, family="ArialMT")

par(mfrow = c(1,1), mar = c(3,3,0,0), bty='n', oma = c(0,0,0,0))

plot(x = range(1:7), y = range(0:6), type = 'n', 
     xaxt = 'n', yaxt = 'n')
axis(1, tcl = -0.1, mgp = c(3,0.1,0))
axis(2, tcl = -0.1, mgp = c(3,0.2,0), las = 2, at = c(0:5), labels = c(1, 10, expression(10^2),
                                                                       expression(10^3),
                                                                       expression(10^4),
                                                                       expression(10^5)))
mtext("iRBC abundance", side=2, line=1.6, cex = cexText)
mtext("Day of blood-stage infection", side=1, line=1.2, cex=cexText)

for (id in c(1:length(mtParaCountsExcl[,1]))){
  subDataParaToPlot = mtParaCountsExcl[id, ]
  subDataToPlot = data.frame(day = (1:7)[!is.na(subDataParaToPlot)], 
                             para = subDataParaToPlot[!is.na(subDataParaToPlot)])
  points(subDataToPlot$day, log10(subDataToPlot$para), type = 'l', col = lineCol)
}
text(x = 1, y = 6, "MT data (excluded patients)", adj=0)

dev.off()

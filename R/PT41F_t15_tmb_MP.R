#------------------------------------------------------------------------------
# Pella-Tomlinson 40:10-type MPs (details implemented below)
# highest level function sets options including priors and TAC change constraints
# mid-level function does TMB interface, minimization (and grid-search if req'd)
# low level function is written in TMB / C++
#------------------------------------------------------------------------------
# F-based hockey-stick 15% change constraint

PT41F.t15.tmb<-function(pset, BLower=0.1,BUpper=0.4,CMaxProp=1., useF=1, deltaTACLimUp=0.15, deltaTACLimDown=0.15)
{
  return (PT41F.base.tmb(pset,
                         BLower,
                         BUpper,
                         CMaxProp,
                         useF,
                         deltaTACLimUp,
                         deltaTACLimDown))
}

class(PT41F.t15.tmb)<-"IO_MP_tune"

#------------------------------------------------------------------------------

PT41F.base.tmb<-function(pset, BLower=0.1,BUpper=0.4,CMaxProp=1., useF=1, deltaTACLimUp=0.25, deltaTACLimDown=0.25, MShape=TRUE){
  # define priors
  initDepMode <- 1. # (initial depletion = B/K)
  #set MSY prior Mode (and initial value) to 80% of mean of 5 highest catches
  MSYMode    <- 0.8*mean(sort(pset$Cobs, decreasing=T)[1:5]) 
  #set k prior mode relative to MSY
  kMode      <- 20*MSYMode
  #logit transform bounds for k
  kBounds <- c(MSYMode*5, MSYMode*40)
    
  sigmaPMode <- 0.15 #sigma productivity
  sigmaIMode <- 0.15 #sigma Index
  shapeMode  <- -0.16 #shape #-0.16 # do not expect this to be estimable

  Data = list("I_t"=pset$Iobs, "c_t"=pset$Cobs,  
              "priorMode" =c(initDepMode,MSYMode, kMode, shapeMode, sigmaIMode, sigmaPMode), 
              "priorLogCV"=c(0.05,       .25,    5.,      0.05,      0.1,        0.2),
              "log_kBounds"=log(kBounds))

  Binit <- rep(kMode,pset$y)
  
  Params = list("log_MSY"=log(MSYMode), "log_k"=log(kMode), "shape"=shapeMode, "log_sigmaP"=log(sigmaPMode), "log_sigmaI"=log(sigmaIMode),"log_B_t"=log(Binit)) #,  "log_sigmac"=log(0.01), "logit_exploit_t"=rep(0,n_years-1) )
  
  Random = c("log_B_t")
  
  Map = list() # fix parameters at correct simulation values
  if (MShape)
  {
    Map[["shape"]] = factor(NA)
  }

  Map[["log_sigmaI"]] = factor(NA)
  tmbList <- list(Data, Params, Random, Map) 
  names(tmbList) <- c('Data', 'Params', 'Random', 'Map')
    return (PT4010tmb(pset, BLower=BLower,BUpper=BUpper,CMaxProp=pset$tune * CMaxProp, deltaTACLimUp=deltaTACLimUp, deltaTACLimDown=deltaTACLimDown, useF=1,
                   gridSearch = 0, tmbList = tmbList))
}

#------------------------------------------------------------------------------

shouldLogPerformance <- function(pset)
{
  return ((!is.null(pset$MP_environment)                  & 
           exists("TAC",       envir=pset$MP_environment) &
           exists("B",         envir=pset$MP_environment) &
           exists("Depletion", envir=pset$MP_environment) &
           exists("q",         envir=pset$MP_environment)))
}
   
#------------------------------------------------------------------------------
# Function for logging MP performance data
#------------------------------------------------------------------------------
logPerformance <- function(pset, Report, TAC, BY, plots=NA)
{
  if (!is.null(pset$MP_environment)                  & 
      exists("TAC",       envir=pset$MP_environment) &
      exists("B",         envir=pset$MP_environment) &
      exists("Depletion", envir=pset$MP_environment) &
      exists("q",         envir=pset$MP_environment))
  {
    pset$MP_environment$TAC        <- c(pset$MP_environment$TAC, as.double(TAC))
    pset$MP_environment$B          <- c(pset$MP_environment$B, as.double(Report$B_t[length(Report$B_t)]))
    pset$MP_environment$Depletion  <- c(pset$MP_environment$Depletion, as.double(Report$Depletion_t[length(Report$Depletion_t)]))
    pset$MP_environment$q          <- c(pset$MP_environment$q, as.double(Report$q[length(Report$q)]))
    pset$MP_environment$plots      <- plots

    pset$MP_environment$ModelData  <- list(
      BY          = as.double(BY),
      B_t         = as.double(Report$B_t),
      Bpred_t     = as.double(Report$Bpred_t),
      recDev      = as.double(Report$recDev),
      msy         = exp(as.double(Report$log_MSY)),
      k           = as.double(Report$k),
      r           = as.double(Report$r),
      shape       = as.double(Report$shape),
      q           = as.double(Report$q),
      Depletion_t = as.double(Report$Depletion_t),
      log_sigmaI  = as.double(Report$log_sigmaI), 
      log_sigmaP  = as.double(Report$log_sigmaP)
    )
  }
}


#------------------------------------------------------------------------------
# Pella Tomlinson Production model with generic 40-10 type rule - MPs are defined with tuning parameters above
# useF option uses the 40:10 rule for F rather than C, in which case FMax = FMSY*CMaxProp
# positive gridSearch value is preferred at this time
# tmb indicates joint process and observation error model implemented with TMB
#------------------------------------------------------------------------------
PT4010tmb<-function(pset, BLower=0.1,BUpper=0.4,CMaxProp=1.0, deltaTACLimUp=0.9, deltaTACLimDown=0.9, shockAbsorber=1, useF=1, gridSearch=0, tmbList =list){

  diagnose <- 0.01 #0.05 # 0 = skip, >0 is T and rnd pproportion to plot
  
  # bypass PT fitting if CPUE depletion indicates collapse
  I_hist <- pset$Iobs
  Iy <- I_hist[!is.na(I_hist)]
  Idepletion <- mean(Iy[(length(Iy)-1):length(Iy)])/mean(Iy[1:5])
  print("CPUE % depletion: " %&% floor(Idepletion*100))
  if (Idepletion < 0.1){
    print("skip PT fitting - stock clearly needs rebuilding")
    newTAC <- 1.    #i.e. shutdown fishery as fast as possible
  } else {
    
  obj = TMB::MakeADFun( data=tmbList$Data, parameters=tmbList$Params, random=tmbList$Random, map=tmbList$Map, 
                       DLL="PT41F.t15.tmb.MP") #, silent=T)

  #suppress output
  obj$env$tracemgc <- FALSE
  obj$env$inner.control$trace <- FALSE
  obj$env$silent <- TRUE  
  
  # Optimize - no bounds
  Opt = stats::nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr, control=list("trace"=1, "eval.max"=1e4, "iter.max"=1e4))
    
  Opt[["diagnostics"]] = data.frame( "Est"=Opt$par, "final_gradient"=obj$gr(Opt$par) )
  Report = obj$report()
  SD = try(TMB::sdreport( obj ))

  if(!all(is.finite(Report$nll_comp))) browser()   
  if(!all(is.finite(Opt$objective))) browser()   

  if(sum(Report$B_t<0)>0){  #negative biomass
    graphics::plot(Report$B_t)
    print("")
    print("Negative biomass PT fitting failure")
  }
  
  # visualize and save results
  if( all(abs(Opt$diagnostics$final_gradient)<0.01) ){
    if( !inherits(SD,"try-error")) {
      print("convergence okay")
      
    } else {
      print("converged, but some other error likely...")
      Plot_Fn( report=Report, sdsummary=summary(SD), tmbList = tmbList, OMMSY=pset$MSY)
      if(diagnose) browser()
    }
  } else {
    print("")
    print("gradient convergence failure")
    Plot_Fn( report=Report, sdsummary=summary(SD), tmbList = tmbList, OMMSY=pset$MSY)
  }  
    msy <- exp(Report$log_MSY)
    k   <- Report$k
    Y   <- length(Report$B_t)
    SD1 <- summary(SD)
    
    # Current biomass point estimate
    # Current biomass probably simplistic lower confidence bound
    BCI <- 0.4   # 0.25 = lower 25th assuming normal
    BY  <- SD1[rownames(SD1)=="B_t",][Y,1] + stats::qnorm(BCI)* SD1[rownames(SD1)=="B_t",][Y,2]
    
    #Apply the 40:10 rule to F ...
    if(useF){
      FMSY = -log(1-msy/k)
      FMult = CMaxProp # maximum F relative to FMSY
      if(BY / k <= BLower) TACF <- 0.0001    #i.e. (almost) shutdown fishery
      if(BY / k  >  BLower & BY / k  <= BUpper) TACF <- FMult*FMSY*((BY / k )/(BUpper-BLower) + ( 1 - (BUpper/(BUpper-BLower))))^shockAbsorber
      if(BY / k  >  BUpper) TACF <- FMult*FMSY
      newTAC <-  BY*(1-exp(-TACF))
    }

    names(newTAC) <- "TAC"
  } # fit the PT model
  
      
  lastTAC  <- pset$prevTACE$TAC
  deltaTAC <- newTAC/lastTAC - 1

  if(deltaTAC >  deltaTACLimUp)   deltaTAC =  deltaTACLimUp
  if(deltaTAC < -deltaTACLimDown) deltaTAC = -deltaTACLimDown
  newTAC <- lastTAC*(1+deltaTAC)
  if(newTAC<9) newTAC <- 9 #shut the fishery down, except collect some data
  TAEbyF <- 0.0 * pset$prevTACE$TAEbyF #TAE by fishery

  if (min(TAEbyF) < 0)
  {
    print("MP TAEbyF < 0")
  }

  if (min(newTAC) < 0)
  {
    print("MP TAC < 0")
  }

  if (shouldLogPerformance(pset))
  {
    tmbObj <- if (!is.null(pset$MP_environment$calcLikelihoodProfiles) && pset$MP_environment$calcLikelihoodProfiles) obj else NULL
    plots  <- reportPlots(report=Report, sdsummary=summary(SD), tmbList = tmbList, firstYr=pset$firstYr, obj=tmbObj)

    logPerformance(pset, Report, newTAC, BY, plots)
  }

  rm(SD, Report, obj)

  return (list(TAEbyF=pset$prevTACE$TAEbyF,TAC=newTAC))
}

#------------------------------------------------------------------------------

# MP plot func for testing
Plot_Fn <- function( report, sdsummary, tmbList, OMMSY){
  graphics::par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
  Y <- length(report$B_t)
  data <- tmbList$Data
  # Biomass
  graphics::plot(c(report$B_t/1000), type="l", 
           ylim=c(0,report$k*1.5/1000), ylab="", xlab="year", main="Biomass & CPUE")
  Mat = cbind( report$B_t, sdsummary[which(rownames(sdsummary)=="B_t"),"Std. Error"])
  graphics::polygon( x=c(1:Y,Y:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2]))/1000, col=grDevices::rgb(1,0,0,0.2), border=NA)
  graphics::points(data$I_t/report$q/1000, col=1, pch=15)

  # Depletion
  graphics::plot(report$Depletion_t, type="l", ylim=c(-.4,1.5), ylab="", xlab="year", main="B/K & Prod devs")
  Mat = cbind( report$Depletion_t, sdsummary[which(rownames(sdsummary)=="Depletion_t"),"Std. Error"])
  graphics::polygon( x=c(1:Y,Y:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2])), col=grDevices::rgb(1,0,0,0.2), border=NA)

  Mat2 = cbind(report$recDev, sdsummary[which(rownames(sdsummary)=="recDev"),"Std. Error"])
  graphics::polygon( x=c(1:Y,Y:1), y=c(Mat2[,1]+Mat2[,2],rev(Mat2[,1]-Mat2[,2])), col=grDevices::rgb(0,0,1,0.2), border=NA)
  
  graphics::lines(report$recDev)
  graphics::abline(h=0., lty=2)
  graphics::abline(h=1., lty=2)
  
  # Catch
  graphics::plot(1:Y, data$c_t/1000, type="l", ylab="", xlab="year", main="Catch")
  
  #production function
  B   <- seq(0.01, report$k, report$k/100)
  RG1 <- ((report$shape+1)/report$shape)*report$r*B*(1-(abs(B/report$k))^report$shape) 
  BMSYoK <- floor(100*B[RG1==max(RG1)]/report$k)
  msy <- floor(exp(report$log_MSY)/1000)
  k   <- floor(report$k/1000)
  msyok  <- floor(1000*msy/k)/10
  OMMSY  <- floor(OMMSY/1000)
  graphics::plot(B/1000,RG1/1000, type='l', cex.main=0.7,
       main = "MSY: " %&% msy %&% "    k: " %&% k %&% "  OM-MSY: " %&% OMMSY
           %&% "\nBMSY/K: " %&% BMSYoK %&% "%  MSY/k: " %&% msyok %&% "%")
}

#------------------------------------------------------------------------------

Plot_FnProj = function( report, reportProj, sdsummary, sdsummaryProj, tmbList, OMMSY){

  graphics::par( mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
  Y <- length(report$B_t)
  data <- tmbList$Data
  # Biomass
  graphics::plot(c(reportProj$B_t/1000,reportProj$BProj_t/1000), type="l", 
       ylim=c(0,report$k*1.5/1000), ylab="", xlab="year", main="Biomass & CPUE")
  Mat = cbind( report$B_t, sdsummary[which(rownames(sdsummary)=="B_t"),"Std. Error"])
  graphics::polygon( x=c(1:Y,Y:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2]))/1000, col=grDevices::rgb(1,0,0,0.2), border=NA)
  graphics::points(data$I_t/report$q/1000, col=1, pch=15)
  
  # Depletion
  graphics::plot(c(report$Depletion_t, reportProj$BProj_t/report$k), type="l", ylim=c(-.4,1.5), ylab="", xlab="year", main="B/K & Prod devs")
  Mat = cbind( report$Depletion_t, sdsummary[which(rownames(sdsummary)=="Depletion_t"),"Std. Error"])
  graphics::polygon( x=c(1:Y,Y:1), y=c(Mat[,1]+Mat[,2],rev(Mat[,1]-Mat[,2])), col=grDevices::rgb(1,0,0,0.2), border=NA)
  
  Mat2 = cbind(report$recDev, sdsummary[which(rownames(sdsummary)=="recDev"),"Std. Error"])
  graphics::polygon( x=c(1:Y,Y:1), y=c(Mat2[,1]+Mat2[,2],rev(Mat2[,1]-Mat2[,2])), col=grDevices::rgb(0,0,1,0.2), border=NA)
  
  graphics::lines(report$recDev)
  graphics::abline(h=0., lty=2)
  graphics::abline(h=1., lty=2)
  
  # Catch
  graphics::plot(c(data$c_t/1000, rep(reportProj$newTAC/1000,length(reportProj$BProj_t))), type="l", ylab="", xlab="year", main="Catch", col=2)
  graphics::lines(data$c_t/1000)
  
  #production function
  B   <- seq(0.01, report$k, report$k/100)
  RG1 <- ((report$shape+1)/report$shape)*report$r*B*(1-(abs(B/report$k))^report$shape) 
  BMSYoK <- floor(100*B[RG1==max(RG1)]/report$k)
  msy <- floor(exp(report$log_MSY)/1000)
  k   <- floor(report$k/1000)
  msyok  <- floor(1000*msy/k)/10
  OMMSY  <- floor(OMMSY/1000)
  graphics::plot(B/1000,RG1/1000, type='l', cex.main=0.7,
       main = "MSY: " %&% msy %&% "    k: " %&% k %&% "  OM-MSY: " %&% OMMSY
       %&% "\nBMSY/K: " %&% BMSYoK %&% "%  MSY/k: " %&% msyok %&% "%")
}

#------------------------------------------------------------------------------

reportPlots <- function(report, sdsummary, tmbList, firstYr, obj=NULL)
{
  Y       <- firstYr - 1 + (1:length(report$B_t))
  data    <- tmbList$Data
  colors  <- c("Catch"="#001D34", "TAC"="#00A9CE")

  # CPUE
  cpue_serr  <- as.double(sdsummary[which(rownames(sdsummary)=="B_t"), "Std. Error"]) * report$q
  cpue       <- as.double(report$B_t * report$q)
  cpue_lower <- cpue - cpue_serr
  cpue_upper <- cpue + cpue_serr
  cpue_data  <- data.frame(t=Y, cpue_t=cpue, lower=cpue_lower, upper=cpue_upper, cpue=data$I_t)

  cpue_plot  <- ggplot2::ggplot(data=cpue_data, ggplot2::aes(x=cpue_data$t, y=cpue_data$cpue_t)) + 
                                ggplot2::geom_line(colour=colors[1], size=2, alpha=0.5) +
                                ggplot2::geom_point(color="black", shape=1, size=5, mapping=ggplot2::aes(x=cpue_data$t, y=cpue_data$cpue)) + 
                                ggplot2::geom_ribbon(data=cpue_data, ggplot2::aes(x=cpue_data$t, ymin=cpue_data$lower, ymax=cpue_data$upper), alpha=0.07) + 
                                ggplot2::ggtitle("CPUE") +
                                ggplot2::xlab('\n Year') +
                                ggplot2::ylab('CPUE\n ') +      
                                ggplot2::theme_bw()

  # Biomass
  biomass_serr  <- as.double(sdsummary[which(rownames(sdsummary)=="B_t"), "Std. Error"]) / 1000
  biomass       <- as.double(report$B_t) / 1000
  biomass_lower <- biomass - biomass_serr
  biomass_upper <- biomass + biomass_serr
  biomass_cpue  <- as.double(data$I_t / report$q) / 1000
  biomass_data  <- data.frame(t=Y, B_t=biomass, lower=biomass_lower, upper=biomass_upper, B_cpue=biomass_cpue)
  biomass_plot  <- ggplot2::ggplot(data=biomass_data, ggplot2::aes(x=biomass_data$t, y=biomass_data$B_t)) + 
                                   ggplot2::geom_line(colour=colors[1], size=2, alpha=0.5) +
                                   ggplot2::geom_point(color="black", shape=1, size=5, mapping=ggplot2::aes(x=biomass_data$t, y=biomass_data$B_cpue)) + 
                                   ggplot2::geom_ribbon(data=biomass_data, ggplot2::aes(x=biomass_data$t, ymin=biomass_data$lower, ymax=biomass_data$upper), alpha=0.07) + 
                                   ggplot2::ggtitle("Biomass") + 
                                   ggplot2::xlab('\n Year') + 
                                   ggplot2::ylab('Biomass\n ') +  
                                   ggplot2::theme_bw()

  # Depletion
  depletion       <- as.double(report$Depletion_t)
  depletion_serr  <- as.double(sdsummary[which(rownames(sdsummary)=="Depletion_t"), "Std. Error"])
  depletion_lower <- depletion - depletion_serr
  depletion_upper <- depletion + depletion_serr
  depletion_data  <- data.frame(t=Y, depletion_t=depletion, lower=depletion_lower, upper=depletion_upper)
  depletion_plot  <- ggplot2::ggplot(data=depletion_data, ggplot2::aes(x=depletion_data$t, y=depletion_data$depletion_t)) + 
                                     ggplot2::geom_line(colour=colors[1], size=2, alpha=0.5) +
                                     ggplot2::geom_ribbon(data=depletion_data, ggplot2::aes(x=depletion_data$t, ymin=depletion_data$lower, ymax=depletion_data$upper), alpha=0.07) + 
                                     ggplot2::ggtitle("Depletion") + 
                                     ggplot2::xlab('\n Year') + 
                                     ggplot2::ylab('Depletion\n ') +  
                                     ggplot2::theme_bw()

  # Recruitment deviations
  recDev       <- as.double(report$recDev)
  recDev_serr  <- as.double(sdsummary[which(rownames(sdsummary)=="recDev"), "Std. Error"])
  recDev_lower <- recDev - recDev_serr
  recDev_upper <- recDev + recDev_serr
  recDev_data  <- data.frame(t=Y, recDev_t=recDev, lower=recDev_lower, upper=recDev_upper)
  recDev_plot  <- ggplot2::ggplot(data=recDev_data, ggplot2::aes(x=recDev_data$t, y=recDev_data$recDev_t)) + 
                                  ggplot2::geom_line(colour=colors[1], size=2, alpha=0.5) +
                                  ggplot2::geom_ribbon(data=recDev_data, ggplot2::aes(x=recDev_data$t, ymin=recDev_data$lower, ymax=recDev_data$upper), alpha=0.07) + 
                                  ggplot2::ggtitle("Recuitment Deviations") + 
                                  ggplot2::xlab('\n Year') +
                                  ggplot2::ylab('Recruitment Deviation\n ') +      
                                  ggplot2::theme_bw()

  # Production function
  B         <- seq(0.01, as.double(report$k), as.double(report$k) / 100)
  RG1       <- as.double((report$shape + 1)/ report$shape) * report$r * B * (1.0 - (abs(B / report$k)) ^ report$shape)
  BMSYoK    <- as.double(floor(100 * B[RG1==max(RG1)] / report$k))
  msy       <- as.double(floor(exp(report$log_MSY) / 1000))
  k         <- as.double(floor(report$k / 1000))
  msyok     <- as.double(floor(1000 * msy / k) / 10)
  Title     <- paste("MSY:", msy, " k:", k, "BMSY/K: ", BMSYoK, "%  MSY/k:", msyok, "%")
  prod_data <- data.frame(B=B / 1000, RG1=RG1 / 1000)

  prod_plot <- ggplot2::ggplot(data=prod_data, ggplot2::aes(x=B, y=RG1)) + 
                               ggplot2::geom_line(colour=colors[1], size=2, alpha=0.5) +
                               ggplot2::ggtitle(Title) + 
                               ggplot2::theme_bw()

  plots <- list(cpue_plot=cpue_plot, biomass_plot=biomass_plot, depletion_plot=depletion_plot, recDev_plot=recDev_plot, prod_plot=prod_plot)

  if (!is.null(obj))                 
  {
    profile_plots <- list()
    bestNames     <- names(obj$env$last.par.best)
    params        <- unique(bestNames)

    for (param in params)
    {
      if (sum(param == bestNames) == 1)
      {
        profile   <- TMB::tmbprofile(obj, param)
        prof_data <- data.frame(x=as.double(profile[[param]]), value=profile$value)
        prof_plot <- ggplot2::ggplot(data=prof_data, ggplot2::aes(x=x, y=value)) + 
                                     ggplot2::geom_line(colour=colors[1], size=2, alpha=0.5) + 
                                     ggplot2::xlab(param) + 
                                     ggplot2::theme_bw()

        profile_plots[[param]] <- prof_plot
      }
    }

    plots[["profiles"]] <- profile_plots
  }

  return (plots)
}


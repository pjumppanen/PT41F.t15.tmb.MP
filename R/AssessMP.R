# -----------------------------------------------------------------------------
# named cache for the loaded compiled MP libraries. We require this because 
# when run in a shiny context because of dangling pointers left by TMB if the
# library is unloaded it crashes the R session. To avoid this we never unload
# the library but we need to keep track of when it is loaded so we don't load
# multiple times which would result in the previous instance being unloaded. 
# Hopefully TMB does not have memory leaks as this could potentially result in
# memory exhaustion if that is the case.
# -----------------------------------------------------------------------------
MP_cache         <- rlang::new_environment()
MP_cache$sources <- list()


# -----------------------------------------------------------------------------
# assessMP
# -----------------------------------------------------------------------------
# this function applies the designated MP to the catch and cpue series data
# in the supplied csv file and measures performance statistics which are 
# returned in an environment. The reported statistics are:
#
#   TAC
#   B
#   Depletion
#   q
# 
# Typical example usage:
#
# source("./Source/AssessMP.R")
#
# theta           <- 1.0 #tuning parameter
# MP_Interval     <- 3
# LastTAC         <- 100000
# CatchAndCPUEcsv <- "./Source/test.csv"
# MP_SourcePath   <- "./MPs/PTTMB/MPs_TMBMSY_tidied.R"
# MP_Name         <- "PTBoB0Targ.t15.pr15"
# 
# results <- assessMP(MP_Name, MP_SourcePath, CatchAndCPUEcsv, MP_Interval, theta)
#
# print(results$TAC)
# print(results$B)
# print(results$Depletion)
# print(results$q)
#
# -----------------------------------------------------------------------------
assessMP <- function(CatchAndCPUEcsv, calcLikelihoodProfiles=FALSE)
{
  CatchAndCPUE <- read.csv(CatchAndCPUEcsv)
  Names        <- names(CatchAndCPUE)
  theta        <- 3.718027
  MP_Interval  <- 1

  # Check for y column
  if (!any(Names == "Year"))
  {
    stop("Year Column missing (year)")
  }

  # Check for Catch column
  if (!any(Names == "Catch"))
  {
    stop("Catch Column missing (observed catch)")
  }

  # Check for CPUE column
  if (!any(Names == "CPUE"))
  {
    stop("CPUE Column missing (observed CPUE)")
  }

  # Check for TAC column
  if (!any(Names == "TAC"))
  {
    return("TAC Column missing (historic TAC)")
  }

  # Check for ascending contiguous years
  MinYear <- min(CatchAndCPUE$Year)
  MaxYear <- max(CatchAndCPUE$Year)
  y       <- MinYear

  for (idx in 1:length(CatchAndCPUE$Year))
  {
    if (CatchAndCPUE$Year[idx] != y)
    {
      stop("y must be in contiguous ascending years")
    }

    y <- y + 1
  }

  # Check MP_interval
  if (MP_Interval < 1)
  {
    stop("MP_Interval must be >= 1")
  }

  # Find LastTAC
  LastTAC <- CatchAndCPUE$TAC[length(CatchAndCPUE$TAC)]

  if (is.null(LastTAC) || is.na(LastTAC))
  {
    stop("MP requires a TAC value in the last recorded year of data")
  }

  # Evaluate selected MP over all available MP years
  nfleets                               <- 1
  y                                     <- MaxYear
  MP_environment                        <- rlang::new_environment()
  MP_environment$TAC                    <- c()
  MP_environment$B                      <- c()
  MP_environment$Depletion              <- c()
  MP_environment$q                      <- c()
  MP_environment$plots                  <- c()
  MP_environment$calcLikelihoodProfiles <- calcLikelihoodProfiles
  TAE                                   <- keep::karray(rep(0, nfleets), dim=c(nfleets))
  TACE                                  <- list(TAEbyF=TAE, TAC=LastTAC)

  pset <- list(y              = y - MinYear + 1,
               firstYr        = MinYear,
               Cobs           = CatchAndCPUE$Catch,
               Iobs           = CatchAndCPUE$CPUE,
               tune           = theta,
               MSY            = NA,
               prevTACE       = TACE,
               interval       = MP_Interval,
               MP_environment = MP_environment)

  # Call MP with pset set of data
  TACE <- PT41F.t15.tmb(pset)

  return (MP_environment)
}

\name{assessMP}
\alias{assessMP}
\title{
  Use the MP to assess the future TAC for a given historic catch and CPUE time series
}
\description{
  Calling \code{assessMP} applies the \emph{PT41F.t15.tmb} MP to the supplied Catch, CPUE and TAC data to produce a new TAC recommendation for
  future catches. 
}
\usage{
  assessMP(CatchAndCPUEcsv, calcLikelihoodProfiles=FALSE)
}
\arguments{
  \item{CatchAndCPUEcsv}{
    \code{CatchAndCPUEcsv} is the file name and path of the csv file containing Catch, CPUE and TAC series data with the expected column
    names \code{Year}, \code{CPUE}, \code{Catch} and \code{TAC}.
  }
  \item{calcLikelihoodProfiles}{
    \code{Force} is a boolean that when \code{TRUE} calculates likelihood profiles for the parameters in the Pella-Tomlinson model fit.
  }
}
\details{
  Calling \code{assessMP} applies the \emph{PT41F.t15.tmb} MP to the supplied Catch, CPUE and TAC data to produce a new TAC recommendation for
  future catches. The data needed for the MP evaulation is supplied in the \code{CatchAndCPUEcsv} file and must have the columns \code{Year}, 
  \code{CPUE}, \code{Catch} and \code{TAC}. The \code{CPUE} column data is obtained from the output of the \code{\link{assessmentCPUE}} function.
  For the \code{TAC} series, only the last year TAC figure for historic catches is required and is needed to apply the +/- 15\% TAC change 
  restriction to the recommended future TAC. The MP results are returned in a enviroment with \code{B}, \code{Depletion}, \code{q}, \code{TAC} 
  members outlining the recommendation and \code{ModelData} showing the entire TMB Pella-Tomlinson model fit and \code{plots} showing \code{ggplot} 
  objects for pertinant results. 
}
\value{
  Returns an environment containing the MP assessment results.
}
\author{
  \packageAuthor{PT41F.t15.tmb.MP}

  Maintainer: \packageMaintainer{PT41F.t15.tmb.MP}
}

\seealso{
  \code{\link{PT41F.t15.tmb.MP-package}}
}

\note{
  see \code{\link{PT41F.t15.tmb.MP-package}} for examples of use.
}

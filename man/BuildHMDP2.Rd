\name{BuildHMDP2}
\alias{BuildHMDP2}
\title{Build the HMDP (2 levels with shared linking) using the C++ binary writer.}
\usage{
BuildHMDP2(filePrefix, param, dlms, dglmParam)
}
\arguments{
  \item{filePrefix}{Prefix used by the binary files storing
  the MDP model.}

  \item{param}{Model parameters a list created using
  \code{\link{setParameters}}.}

  \item{dlms}{A list a long as the number of rations where
  each entry is a DLM created using \code{\link{iniDLM}}
  and \code{\link{buildDLM}}.}

  \item{dglmParam}{Model parameters a list created using
  \code{\link{iniDGLM}}.}
}
\value{
Build log (character).
}
\description{
Build the HMDP (2 levels with shared linking) using the C++
binary writer.
}
\author{
Lars Relund \email{lars@relund.dk}
}


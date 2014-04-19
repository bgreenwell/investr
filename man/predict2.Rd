\name{predict2}
\alias{predict2}
\alias{predict2.lm}
\alias{predict2.nls}
\title{Predict method for (Single-Regressor) Linear and Nonlinear Model Fits}
\usage{
predict2(object, ...)

\method{predict2}{lm}(object, newdata, interval = c("none", "confidence",
  "prediction"), level = 0.95, adjust = c("none", "Bonferroni", "Scheffe"),
  k, ...)

\method{predict2}{nls}(object, newdata, interval = c("none", "confidence",
  "prediction"), level = 0.95, adjust = c("none", "Bonferroni", "Scheffe"),
  k, ...)
}
\description{
Convenience function to be called by \code{plotFit}. It is
not for routine use.
}
\keyword{internal}


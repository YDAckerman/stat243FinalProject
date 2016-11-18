print.mfa <- function(x){
  cat('object: "mfa"\n\n')
  cat('number of tables:', length(mfa1$sets), '\n')
  cat('number of observations:', length(x$eigenValues), '\n')
  cat('number of random variables:', length(mfa1$a), '\n')
  cat('First Two Eigenvalues:', x$eigenValues[1:2], '\n')
  cat(sprintf('Variance Explained by first two eigenvalues: (%.1f, %.1f)'
              , 100*x$eigenValues[1] / sum(x$eigenValues)
              , 100*x$eigenValues[2] / sum(x$eigenValues))
      , '\n\n')
  cat('First Two Factor Scores\n')
  print(x$factorScores[,1:2])
  invisible(x)
}
#' plot_variable_loadings
#'
#' Plot the compromise scores for the first two extracted dimensions
#' @param mfa - the mfa object
#' @param title - chart title
#' @param dims - Numeric vector of two values signifying dimensions to plot
#' @return a plot with the variable loadings for each table
#' @export

plot_variable_loadings <- function(mfa, dims = 1:2, title = 'Variable Loadings'){
  
  ## perform checking
  check_mfa(mfa)
  check_plot_dims(length(attributes(mfa)$sets), dims)
  
  # Set plotting parameters
  tables <- length(mfa$partialFactorScores)
  cols <- 5 # Number of columns in the grid
  rows <- ceiling(tables / cols) # Number of rows in the grid
  par(mfrow = c(rows, cols)
      , oma = c(4, 2, 4, 2) # Outer margin (b, r, u, l)
      , mar = c(2, 2, 2, 2) ) # Inner margins (b, r, u, l)
  
  # Min and max plotting dimensions
  x_range <- c(min(mfa$Q[,dims[1]]), max(mfa$Q[,dims[1]]))
  y_range <- c(min(mfa$Q[,dims[2]]), max(mfa$Q[,dims[2]]))
  
  # Loop through and plot each table
  for (k in 1:tables){
    
    # Subset the relevant matrix
    loadings_k <- mfa$Q[attributes(mfa)$sets[[k]],]

    
    # Plot the data
    plot(loadings_k[,dims[1]]
         , loadings_k[,dims[2]]
         , pch = 1
         , cex = 2
         , xlim = x_range
         , ylim = y_range
         , xlab = ''
         , ylab = '')
    
    # Text to identify the table number
    text(x = x_range[2] * .8
         , y = y_range[2]
         , labels = paste0('k=', k))
    
    # Text to identify the variable if it exists
    if(!is.null(attributes(mfa)$var_names)){
        text(x = loadings_k[,dims[1]]
             , y = loadings_k[,dims[2]]
             , labels = attributes(mfa)$var_names[[k]]
             , pos = 1
        )
    }
    
    # Plot center lines
    abline(v = 0, h = 0)
    }
  mtext(title, outer = TRUE, cex = 1.5)
}

# plot_variable_loadings(mfa1)
# plot_variable_loadings(mfa2)

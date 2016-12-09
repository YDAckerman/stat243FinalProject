# Multiple Factor Analysis in R

> Multiple factor analysis (MFA, also called multiple factorial analysis) is an extension
> of principal component analysis (PCA) tailored to handle multiple data tables that
> measure sets of variables collected on the same observations…

This description is from the paper that this project is based upon: [Multiple factor analysis: principal component analysis for multitable and multiblock data sets](https://www.utdallas.edu/~herve/abdi-WiresCS-mfa-2013.pdf) by Hervé Abdi, Lynne J. Williams and Domininique Valentin (2013)  

### Contents

* [R Code](mfa/R)
* [Unit Tests](mfa/tests)
* [Slides](slides/slides.md)
* [Vignettes](/mfa/vignettes)
* [Shiny](https://mfashinyapp.shinyapps.io/MFA_Shiny_App/)
* [License](./LICENSE.txt)

### Quick Start Guide

```R
# Install and load package through devtools
library(devtools)
devtools::install_github("fussballball/stat243FinalProject/mfa", 
                         force_deps = FALSE)
library(mfa)

# Fit an MFA model on the wines dataset
data(wine)                                    
i <- grep("V", colnames(wine))                
wine_ratings <- wine[,i]                              
expert_sets <- list(1:6, 7:12, 13:18, 19:23, 24:29,  
                    30:34, 35:38, 39:44, 45:49, 50:53)
wine_region <- substr(wine$ID, 1, 2)
wine_names <- as.character(wine$ID)
MFA <- mfa(data = wine_ratings, sets = expert_sets
          , color = wine_region, ids = wine_names)

# Get a summary of the MFA object
print(MFA)

# Print out the eigenvalue summary
summary_eigenvalues(MFA)

# Plot the compromise scores, partial factor scores, and variable loadings
plot(MFA)
```

### Package Developers

* Dario Cantore
* Josiah Davis
* Yanli Fan
* Yoni Ackerman


### References

* [Multiple factor analysis: principal component analysis for multitable and multiblock data sets](https://www.utdallas.edu/~herve/abdi-WiresCS-mfa-2013.pdf) by Hervé Abdi, Lynne J. Williams and Domininique Valentin (2013)  
* [Singular Value Decomposition (SVD) and Generalized Singular Value Decomposition (GSVD)](http://www.cimat.mx/~alram/met_num/clases/Abdi-SVD2007-pretty.pdf) by Hervé Abdi (2007)
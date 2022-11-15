#  multiLODmice

Custom functions to extend mice for left-censored predictors with varying LODs.

To install:
```R
install.packages('devtools')
devtools::install_github('glenmcgee/multiLODmice')
```

Based on functions MImpute_lcens, mice.impute.cens, and .cens.draw3 from doMIsaul package by Faucheux et al (https://github.com/LilithF/doMIsaul).

The method is based on extended version of Lapidus et al (2014; Statistics in Medicine) extended to allow different LODs for different observations (e.g. if LODs vary over time).

GNU General Public License v3.0 (https://github.com/LilithF/doMIsaul/blob/main/LICENSE.md)




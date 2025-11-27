# Bit scale
A package used to convert the theta scale of an IRT model (fitted using the [mirt package](https://cran.r-project.org/web/packages/mirt/index.html)) into a bit scale. Bit scores can be computed directly from estimated thetas.

The package can be installed directly from github using the devtools package as follows:
```R
devtools::install_github("joakimwallmark/bitscale")
```
If you don't have devtools installed yet, you can install it via:
```R
install.packages("devtools")
```
After installing, refer to the package help files for more information on how to use the package:
```R
library(bitscale)
?bit_scores
?bit_plot_item
help(package = "bitscale")
```

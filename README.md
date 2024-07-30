---
output:
  pdf_document: default
  html_document: default
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version-last-release/ModTools)](https://CRAN.R-project.org/package=ModTools) [![downloads](https://cranlogs.r-pkg.org/badges/grand-total/ModTools)](https://CRAN.R-project.org/package=ModTools) [![downloads](http://cranlogs.r-pkg.org/badges/last-week/ModTools)](https://CRAN.R-project.org/package=ModTools) [![License: GPL v2+](https://img.shields.io/badge/License-GPL%20v2+-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html) [![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html) [![R build status](https://github.com/AndriSignorell/ModTools/workflows/R-CMD-check/badge.svg)](https://github.com/AndriSignorell/ModTools/actions) [![pkgdown](https://github.com/AndriSignorell/ModTools/workflows/pkgdown/badge.svg)](https://andrisignorell.github.io/ModTools/)

<!-- badges: end -->

# Tools for Building Regression and Classification Models

There is a rich selection of R packages implementing algorithms for classification and regression tasks out there. The authors legitimately take the liberty to tailor the function interfaces according to their own taste and needs. For us other users, however, this often results in struggling with user interfaces, some of which are rather weird - to put it mildly - and almost always different in terms of arguments and result structures. ModTools pursues the goal of offering uniform handling for the most important regression and classification models in applied data analyses. The function FitMod() is designed as a simple and consistent interface to these original functions while maintaining the flexibility to pass on all possible arguments. print, plot, summary and predict operations can so be carried out following the same logic. The results will again be reshaped to a reasonable standard. For all the functions of this package Google styleguides are used as naming rules (in absence of convincing alternatives). The ’BigCamelCase’ style has been consequently applied to functions borrowed from contributed R packages as well.

As always: Feedback, feature requests, bug reports and other suggestions are welcome! Please report problems to [GitHub issues tracker](https://github.com/AndriSignorell/ModTools/issues) (preferred) or directly to the maintainer.

## Installation

You can install the released version of **ModTools** from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("ModTools")
```

And the development version from GitHub with:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("AndriSignorell/ModTools")
```

# Warning

**Warning:** This package is still under development. Although the code seems meanwhile quite stable, until release of version 1.0 you should be aware that everything in the package might be subject to change. Backward compatibility is not yet guaranteed. Functions may be deleted or renamed and new syntax may be inconsistent with earlier versions. By release of version 1.0 the “deprecated-defunct process” will be installed.

# Authors

Andri Signorell\
Helsana Versicherungen AG, Health Sciences, Zurich\
HWZ University of Applied Sciences in Business Administration Zurich.

R is a community project. This can be seen from the fact that this package includes R source code and/or documentation previously published by [various authors and contributors](https://andrisignorell.github.io/ModTools/authors.html). Special thanks go to Beat Bruengger, Mathias Frueh, Daniel Wollschlaeger for their valuable contributions and testing. The good things come from all these guys, any problems are likely due to my tweaking. Thank you all!

**Maintainer:** Andri Signorell

# Examples

``` r
library(ModTools)
```

<!-- ## Demo "describe" -->

``` r
demo(describe, package = "ModTools")
```

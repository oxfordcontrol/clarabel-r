# clarabel <img src="man/figures/logo.png" width="100" align="right" />

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/clarabel)](https://cran.r-project.org/package=clarabel)

R interface to the
[Clarabel](https://oxfordcontrol.github.io/ClarabelDocs/stable/)
interior point conic optimization solver from the [Oxford Control
Group](https://github.com/oxfordcontrol).

Stable versions can be installed from CRAN as usual. Development
versions from this repo can be installed via:

```
## Install remotes packages if not available
if (! "remotes" %in% installed.packages()[, 1] ) {
	install.packages("remotes", repository = "https://cran.r-project.org")
}
remotes::install_github("bnaras/clarabel")
```

The above code assumes availability of the Cargo/Rust infrastructure
and R compilation tools to install from source, all freely available.

Vignettes are provided and may be perused at the [clarabel package
site](https://bnaras.github.io/clarabel/articles/clarabel.html).


# clarabel <img src="man/figures/logo.png" width="100" align="right" />

R interface to the
[Clarabel](https://oxfordcontrol.github.io/ClarabelDocs/stable/)
interior point conic optimization solver from the [Oxford Control
Group](https://github.com/oxfordcontrol).

Until this gets on CRAN, where binaries will be provided for
various platforms, one can install via:

```
## Install remotes packages if not available
if (! "remotes" %in% installed.packages()[, 1] ) {
	install.packages("remotes", repository = "https://cran.r-project.org")
}
remotes::install_github("bnaras/clarabel")
```

The above code assumes availability of the Cargo/Rust infrastructure
and R compilation tools to install from source, all freely available.

Examples code to run (for now) may be found
[here](https://bnaras.github.io/clarabel/articles/clarabel.html).


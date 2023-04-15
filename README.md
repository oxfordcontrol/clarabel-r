# Clarabel

R interface to the Interior Point Conic Optimization Solver from the
[Oxford Control
Group](https://oxfordcontrol.github.io/ClarabelDocs/stable/).

Until this gets on CRAN, you can install via:

```
## Install remotes packages if not available
if (! "remotes" %in% installed.packages()[, 1] ) {
	install.packages("remotes", repository = "https://cran.r-project.org")
}
remotes::install_github("bnaras/clarabel")
```

Examples code to run (for now) may be found
[here](https://bnaras.github.io/clarabel/articles/clarabel.html).


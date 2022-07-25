# Eutropia

Agent-based metabolic modelling of microbial communities in time and continuous µ-meter-scale space

<!-- badges: start -->
  [![R-CMD-check](https://github.com/Waschina/Eutropia/workflows/R-CMD-check/badge.svg)](https://github.com/Waschina/Eutropia/actions)
<!-- badges: end -->

### What is *Eutropia*

*Eutropia* is a R-package for cell agent-based metabolic modelling of microbial communities. It allows dynamic simulations of two-dimensional surface-attached cell communities. A few features of Eutropia, that you might find interesting:

- Complex polygons (also non-convex) as growth environment
- Extracellular enzymes
- Nutrient regimes
- Chemotaxis (attracting and repelling)
- Direct calculations of community assortment and segregation (see e.g. [Yanni *et al. (2015) Current Biology*](https://doi.org/10.1016/j.cub.2019.03.068))

### Installation

*Eutropia* is in its development phase. The current development version can be installed using: 

```R
# install.packages("devtools")
devtools::install_github("Waschina/Eutropia")
```
If you have not installed `devtools` yet, just uncomment the first line.

### Getting started

For an example simulation, please have a look at the package's vignette.

### Some tricks to increase performance

For the individual Flux-Balance-Analysis steps, Eutropia requires a LP-solver. Ideally, the GLPK (GNU Linear Programming Kit) together with the R-interface-package `glpkAPI` should work. If you like to speed up the simulations, you could consider installing IBM's *ILOG CPLEX Optimization Studio* along with the R-interface-package `cplexAPI`. If you are in academia, make sure to check out [IBM's academic initiative](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students).  Installing the R-package `cplexAPI` can sometimes be a bit tricky, so we'd suggest to get the latest version [here](https://gitlab.cs.uni-duesseldorf.de/general/ccb/cplexAPI) and following the [install instructions](https://gitlab.cs.uni-duesseldorf.de/general/ccb/cplexAPI/-/blob/master/inst/INSTALL).

------

### A big Thank You :bouquet:

- To the developers and maintainers of the *sybil*-universe (`sybil`, `glpkAPI`, `cplexAPI`, and `sybilSBML`). R is a better place with those packages. If you use Eutropia, please make sure to give credits to *sybil* by citing:
  Gelius-Dietrich, G., Desouki, A.A., Fritzemeier, C.J., Lercher, M.J. sybil – Efficient constraint-based modelling in R.                    *BMC Syst Biol* **7,** 125 (2013). https://doi.org/10.1186/1752-0509-7-125
- To everyone behind [GLPK](https://www.gnu.org/software/glpk/) (GNU Linear Programming Kit)


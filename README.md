# Eutropia

Agent-based metabolic modelling of microbial communities in time and continuous µ-meter-scale space

## What is *Eutropia*

*Eutropia* is an R-package for cell agent-based metabolic modelling of microbial communities. It allows dynamic simulations of two-dimensional surface-attached cell communities. A few features of Eutropia, that you might find interesting:

- Complex polygons (also non-convex) as growth environment
- Extracellular enzymes
- Nutrient regimes
- Chemotaxis (attracting and repelling)
- Direct calculations of community assortment and segragation (see e.g. [Yanni *et al. (2015) Current Biology*](https://doi.org/10.1016/j.cub.2019.03.068))

## Getting started

> TODO

## Some tricks to increase performance

For the individual Flux-Balance-Analysis steps, Eutropia requires a LP-solver. Ideally, the GLPK (GNU Linear Programming Kit) together with the R-interface-package `glpkAPI` should work. If you like to speed up the simulations, you could consider installing IBM's *ILOG CPLEX Optimization Studio* along with the R-interface-package `cplexAPI`. If you are in academia, make sure to check out [IBM's academic initiative](https://community.ibm.com/community/user/datascience/blogs/xavier-nodet1/2020/07/09/cplex-free-for-students).  Installing the R-package `cplexAPI` can sometimes be a bit tricky, so we'd suggest to get the latest version [here](https://gitlab.cs.uni-duesseldorf.de/general/ccb/cplexAPI) and following the [install instructions](https://gitlab.cs.uni-duesseldorf.de/general/ccb/cplexAPI/-/blob/master/inst/INSTALL).

Eutropia uses the R-package `particles` as physics engine (e.g. cell collision). The package is available from [CRAN](https://cran.r-project.org/web/packages/particles/index.html). However, it is recommended to install the latest development version directly from github. You can do so with the following commands (if you haven't installed `devtools` yet, simply un-comment the first line):

```R
# install.packages("devtools")
devtools::install_github("thomasp85/particles")
```



------

## A big "Thank you" :bouquet:

- To the developers and maintainers of the *sybil*-universe (`sybil`, `glpkAPI`, `cplexAPI`, and `sybilSBML`). R is a better place with those packages. If you use Eutropia, please make sure to give credits to *sybil* by citing:
  Gelius-Dietrich, G., Desouki, A.A., Fritzemeier, C.J., Lercher, M.J. sybil – Efficient constraint-based modelling in R.                    *BMC Syst Biol* **7,** 125 (2013). https://doi.org/10.1186/1752-0509-7-125
- To everyone behind [GLPK](https://www.gnu.org/software/glpk/) (GNU Linear Programming Kit)


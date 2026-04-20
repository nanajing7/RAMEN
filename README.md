# RAMEN

**RAMEN** implements dynamic bipartite additive and multiplicative effects (AME) models using a temporal alternating least squares (ALS) algorithm. It is designed for bipartite panel data in which relationships between two sets of nodes evolve over time. It supports:

- node- and dyad-level covariates
- nodal heterogeneity and higher-order dependence
- temporal smoothing of coefficients, effects, and latent factors
- penalized ALS estimation in place of MCMC

Typical applications include:

- country–industry trade
- legislator–bill co-sponsorship
- other bipartite relational panel data

### Paper replication

The `paper_example_simulation/` folder contains the simulation scripts used to reproduce the results reported in the accompanying paper.

### Acknowledgements

I thank Prof. Jared Edgerton and Prof. Ryan Kennedy for their invaluable guidance and support in the development of this package.

---

# Installation

Currently the package can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("nanajing7/RAMEN")
```

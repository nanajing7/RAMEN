# RAMEN

**RAMEN** implements dynamic bipartite additive and multiplicative effects (AME) models using a temporal alternating least squares (ALS) algorithm.

The package is designed for bipartite panel data where relationships between two sets of nodes evolve over time. It supports:

- latent factor models for bipartite networks
- temporal smoothing of latent positions
- row, column, and dyadic covariates
- multi-start initialization to reduce local optima

Typical applications include:

- supply chain networks
- trade networks
- district–industry employment panels
- other bipartite relational panel data

### Acknowledgements

I thank Prof. Jared Edgerton and Prof. Ryan Kennedy for their invaluable guidance and support in the development of this package.

---

# Installation

Currently the package can be installed from GitHub:

```r
# install.packages("remotes")
remotes::install_github("nanajing7/RAMEN")
```

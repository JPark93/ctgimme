# ctgimme

Pre-Package Version of `ct-gimme`, a package for the analysis of time-series data with continuous-time vector autoregressive models (CT-VAR)
If making use of this code, please cite the following manuscript:

> Park, J. J., Fisher, Z. F., Hunter, M. D., Shenk, C., Russell, M., Molenaar, P. C., & Chow, S. M. (2024). Unsupervised Model Construction in Continuous-Time. Structural Equation Modeling: A Multidisciplinary Journal, 1-23.

Or use the `BibTeX` version below:

```
@article{park2024unsupervised,
  title={Unsupervised Model Construction in Continuous-Time},
  author={Park, Jonathan J and Fisher, Zachary F and Hunter, Michael D and Shenk, Chad and Russell, Michael and Molenaar, Peter CM and Chow, Sy-Miin},
  journal={Structural Equation Modeling: A Multidisciplinary Journal},
  pages={1--23},
  year={2024},
  publisher={Taylor \& Francis}
}
```
## Usage

`ct-gimme` is primarily defined by the function `ctsgimme()` which accepts several arguments as inputs including:

- `varnames`:
  - vector of variable names in dataset used for analysis
- `dataframe`:
  - the dataset as a data.frame object or as a matrix
- `id`:
  - Character string denoting the ID variable; ideally, set to "id"
- `time`:
  - Character string denoting the time index. If you encounter errors, rename the time-column to "Time". This is a bug I'm working on.
  - For each subject, time should begin at t = 0. Each successive time-point encapsulates time since 0
- `ME.var`:
  - Measurement error variances
  - Define a p x p matrix of the expected ME variances
  - Defaults to a diagonal matrix with near-zero variances. Denotes a model with almost no measurement errors
- `PE.var`:
  - Process noise variances
  - Define a p x p matrix of the expected process noise variances
  - Defaults to an Identity matrix
- `cores`:
  - Number of computational cores to use for parallelization; Defaults to 4
  - WARNING: Do not use more cores than you have
- `directory`:
  - Existing folder to save results to
- `ben.hoch`:
  - Alpha-correction for multiple testing; Defaults to TRUE
- `Galpha`:
  - p-value for group-level effects
  - When `ben.hoch = TRUE`, this is the starting value for the first test all other successive tests use a more stringent p-value
- `Ialpha`:
  - p-value for individual-level effects
- `sig.thrsh`:
  - Proportion of sample that must have a statistically significant effect for a path to be added to the "group"-level.
  - Classic GIMME uses 75% threshold; testing for continuous-time suggests that this can be lower for CT-VAR


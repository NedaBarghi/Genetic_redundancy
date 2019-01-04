# Genetic Redundancy fuels polygenic adaptation in Drosophila

This folder contains all the scripts needed to reproduce the results in Barghi N, Tobler R, Nolte V, Jaksic AM, Mallard F, Otte KA, Dolezal M, Taus T, Kofler R & Schlötterer C. (2019) PLoS Biology.

## SNP_Calling_CMH_FET.sh

Command lines for SNP calling, SNP filteration steps, performing CMH and Fisher's exact tests are provided. The provided script for this step: 

```
call_SNP.py 
```

## estimate_Ne.sh

Command lines needed for estimation on effective population size (Ne). The provided scripts for this step:

```
compute_AF_cov_Ne.py
Ne_estimation.R
```
## neutral_simulations_FDR.sh

Command lines required for performing neutral simulations and estimation of false discovery rate corrected q-values of the CMH and Fisher's exact tests. The provided scripts for this step:

```
neutral_simulation_CMH.R
FDR_CMH.R

neutral_simulation_FET.R
FDR_FET.R
```

## Authors
Neda Barghi



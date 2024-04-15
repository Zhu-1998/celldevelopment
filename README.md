## Uncovering underlying physical principles and driving forces of cell differentiation & reprogramming using single-cell transcriptomics

![image](https://github.com/Zhu-1998/celldevelopment/blob/main/Graph.jpg)

# Analysis tutorials
## 1. Downloading and processing the single-cell transcriptomic data
The scRNA-seq raw data can be downloaded from GEO with an accession number.  

## 2. Estimating RNA velocity from scRNA-seq data
For the scRNA-seq for ~1k U2OS-FUCCI data, we can estimate the RNA velocity by [scvelo](https://github.com/theislab/scvelo) or [dynamo](https://github.com/aristoteleo/dynamo-release). 

For the scEU-seq data for ~3k RPE1-FUCCI cells, `dyn.sample_data.scEU_seq_rpe1()` to acquire the processed data, which includes cell cycle clustering and RNA velocity.


## 3. Reconstructing vector field of cell cycle dynamics
We can reconstruct the vector field based RNA velocity with `dyn.vf.VectorField()` by using [dynamo](https://github.com/aristoteleo/dynamo-release). Then, we can calculate the divergence, curl, acceleration and curvature to perform the differential geometry analysis. We can also calculate the jacobian to perform genetic perturbation and inference gene regulatory interaction.

## 4. Quantifying landscape-flux of cell cycle global dynamics and thermodynamics
We can learn an analytical function of vector field from sparse single cell samples on the entire space robustly by `vector_field_function`. Then, we could simulate stochastic dynamics by solving the Langevin equation based analytical function and quantify the non-equilibrium landscape-flux of the cell cycle.

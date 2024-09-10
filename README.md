## Uncovering underlying physical principles and driving forces of cell differentiation & reprogramming using single-cell transcriptomics

![image](https://github.com/Zhu-1998/celldevelopment/blob/main/Graph.jpg)

# Analysis tutorials
## 1. Downloading and processing the single-cell transcriptomic data
The scRNA-seq raw data can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/) or [ENA](https://www.ebi.ac.uk/ena/browser/home) with an accession number.  
Then, [Cellranger](https://github.com/10XGenomics/cellranger) or [kb-python](https://github.com/pachterlab/kb_python) can be used to process the raw data to get the count data, the [Scanpy](https://github.com/scverse/scanpy) or [Seurat](https://github.com/satijalab/seurat) package help to preprocess the count data, such as Filtering, normalization, dimensionality reduction, clustering, etc.

## 2. Estimating RNA velocity from scRNA-seq data
For the conventional scRNA-seq data, which contains information on both spliced (mature mRNA, exon reads) and unspliced (pre-mRNA, intron reads) transcripts. One can utilize [Samtools](https://github.com/samtools/samtools), [Velocyto](https://github.com/velocyto-team/velocyto.py) or [kb-python](https://github.com/pachterlab/kb_python) to quantify spliced counts and unspliced counts for each gene. Then, RNA velocity could be estimated by various toolkits following the protocol under the default parameters, for example, [Velocyto](https://github.com/velocyto-team/velocyto.py), [Scvelo](https://github.com/theislab/scvelo), [CellDancer](https://github.com/GuangyuWangLab2021/cellDancer), [Dynamo](https://github.com/aristoteleo/dynamo-release), etc. It is necessary to adopt the appropriate toolkit according to different biological backgrounds.

However, the conventional RNA velocity relies on the mispriming of intron reads for current single-cell platforms, and thus the intron measures were biased and inaccurate, and it was scaled by the splicing rate and lacked real physical meanings (i.e., molecules/hour). The metabolic labeling-based method which measures both the historical or old and the new and nascent RNA of cells in a controllable way will be a better measurement for RNA velocity and transcriptomic dynamics. [Dynast](https://github.com/aristoteleo/dynast-release) toolkit can be used for analyzing labeling datasets and quantifying new RNA counts and old RNA counts. Then [Dynamo](https://github.com/aristoteleo/dynamo-release) package can estimate the absolute RNA velocity following the protocol under the default parameters.

## 3. Reconstructing the cell fate vector field from RNA velocity
One can reconstruct the cell fate vector field from RNA velocity with `dyn.vf.VectorField` by using [Dynamo](https://github.com/aristoteleo/dynamo-release). Then, the differential geometry (divergence, curl, acceleration, curvature, jacobian, etc) could be analyzed based on the reconstructed vector field. One can also learn an analytical function of vector field from sparse single-cell samples on the entire space robustly by `vector_field_function` function in [Dynamo](https://github.com/aristoteleo/dynamo-release).

## 4. Quantifying landscape-flux of global dynamics and thermodynamics
One could simulate stochastic dynamics by solving the Langevin equation based on the analytical function to get the steady-state probability distribution and quantify the non-equilibrium landscape and flux. 
For example, one can run `mouse_landscape_multi.py` in `./landscape-flux/mouse_retina` to generate the steady-state probability distribution of mouse retina development dynamics. The step can output grid data (`Xgrid.csv`, `Ygrid.csv`), probability distribution data (`p_tra.csv`), and stochastic force distribution data (`mean_Fx.csv`, `mean_Fy.csv`). Then, `plot_landscape_path.m` in `./landscape-flux` can be run to plot the landscape and least action path.

On the other hand, if one model the system with stochastic differential equations, one can run C++ following `toggle.c` in `./landscape-flux/toggle` to generate the steady-state probability distribution data `landscape.txt`, then run `toggle_landscape_path.m` to plot the global dynamics landscape and paths.

## 5. Loop flux decomposition
One can use [Dynamo](https://github.com/aristoteleo/dynamo-release) to calculate the least action path time $t_{ij}$ with `dyn.pd.least_action` between different cell-type transitions ($i \neq j$), the transition rate $P_{ij}$ between any two different cell-type transitions ($i \neq j$) is the reciprocal of time, $P_{ij}=1/t_{ij}$. And the transition matrix $M_{ij}$ = $P_{ij}$ ($i \neq j$), $$M_{ii}=-\sum_{j=1}^n P_{ij},  j \neq i.$$ Hence, one could calculate the transitions matrix, which satisfies that the sum of each row is equal to 0. Then, the transitions matrix can be used to run `loop_flux.m` to decompose the loop fluxes, which connect different cell states. Of course, one could adopt other appropriate toolkits to calculate the transition matrix according to different biological backgrounds, for example, [Scvelo](https://github.com/theislab/scvelo), [Cellrank](https://github.com/dpeerlab/cellrank), etc.


### Note: 
One can reproduce the results by running `.ipynb` files in `./notebook`, please note that the name of the notebook file may be inconsistent with the officially published article. Of course, users can also load user-defined data for analysis following the tutorials.

## Cite:
Zhu L, Yang S, Zhang K, Wang H, Fang X, Wang J. Uncovering underlying physical principles and driving forces of cell differentiation and reprogramming from single-cell transcriptomics. Proc Natl Acad Sci U S A. 2024 Aug 20;121(34):e2401540121. doi: [10.1073/pnas.2401540121](https://doi.org/10.1073/pnas.2401540121).

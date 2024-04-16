## Uncovering underlying physical principles and driving forces of cell differentiation & reprogramming using single-cell transcriptomics

![image](https://github.com/Zhu-1998/celldevelopment/blob/main/Graph.jpg)

# Analysis tutorials
## 1. Downloading and processing the single-cell transcriptomic data
The scRNA-seq raw data can be downloaded from [GEO](https://www.ncbi.nlm.nih.gov/geo/) or [ENA](https://www.ebi.ac.uk/ena/browser/home) with an accession number.  
Then, [Cellranger](https://github.com/10XGenomics/cellranger) or [kb-python](https://github.com/pachterlab/kb_python) can be used to process the raw data to get the count data, the [Scanpy](https://github.com/scverse/scanpy) or [Seurat](https://github.com/satijalab/seurat) package help to preprocess the count data, such as Filtering, normalization, dimensionality reduction, clustering, etc.

## 2. Estimating RNA velocity from scRNA-seq data
For the conventional scRNA-seq data, which contains information on both spliced (mature mRNA, exon reads) and unspliced (pre-mRNA, intron reads) transcripts. you can utilize [Samtools](https://github.com/samtools/samtools), [Velocyto](https://github.com/velocyto-team/velocyto.py) or [kb-python](https://github.com/pachterlab/kb_python) to quantify spliced counts and unspliced counts for each gene. Then, RNA velocity could be estimated by various toolkits following the protocol under the default parameters, for example, [Velocyto](https://github.com/velocyto-team/velocyto.py), [Scvelo](https://github.com/theislab/scvelo), [CellDancer](https://github.com/GuangyuWangLab2021/cellDancer), [Dynamo](https://github.com/aristoteleo/dynamo-release), etc. It is necessary to adopt the appropriate toolkit according to different biological backgrounds.

However, the conventional RNA velocity relies on the mispriming of intron reads for current single-cell platforms, and thus the intron measures were biased and inaccurate, and it was scaled by the splicing rate and lacked real physical meanings (i.e., molecules/hour). The metabolic labeling-based method which measures both the historical or old and the new and nascent RNA
of cells in a controllable way will be a better measurement for RNA velocity and transcriptomic dynamics. [Dynast](https://github.com/aristoteleo/dynast-release) toolkit can be used for analyzing labeling datasets and quantifying new RNA counts and old RNA counts. Then [Dynamo](https://github.com/aristoteleo/dynamo-release) package can estimate the absolute RNA velocity
following the protocol under the default parameters.

## 3. Reconstructing the dynamical vector field from RNA velocity
You can reconstruct the vector field from RNA velocity with `dyn.vf.VectorField()` by using [Dynamo](https://github.com/aristoteleo/dynamo-release). Then, the differential geometry (divergence, curl, acceleration, curvature, jacobian) could be analyzed based on the reconstructed vector field. You can also learn an analytical function of vector field from sparse single-cell samples on the entire space robustly by `vector_field_function` function in [Dynamo](https://github.com/aristoteleo/dynamo-release).

## 4. Quantifying landscape-flux of global dynamics and thermodynamics
You could simulate stochastic dynamics by solving the Langevin equation based on the analytical function to get the steay-state probability disturbation and quantify the non-equilibrium landscape and flux. 
For example, one can run `mouse_landscape_multi.py` in `./landscape-flux/mouse_retina` to generate the steay-state probability disturbation of mouse retina development dynamics. The step can output grid data (`Xgrid.csv`, `Ygrid.csv`), probability distribution data (`p_tra.csv`) and stochastic force distribution data (`mean_Fx.csv`, `mean_Fy.csv`). Then, `plot_landscape_path.m` in `./landscape-flux` can be runed to plot the landscape and least action path.
On the other hand, if you modeling the system with differential equations

## 5. Loop flux decomposition
You can use Dynamo to calculate the transition rate between any two cell-type transitions, and then 

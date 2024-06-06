# microhomology
The goal of this project is to use model-based statistical inference to identify the extent of microhomology involvement for nucleotide trimming and ligation during V(D)J recomination of adaptive immune receptor loci.

# Install
Everything R 4.1.3 and/or Python 3.11.5 based. Python and R packages that are required can be installed via [miniconda](https://docs.conda.io/en/latest/miniconda.html): 

```bash 
conda env create -f environment.yml
conda activate microhomology_jax
```

You will also need to install IGoR if you wish to annotate sequences using IGoR (Marcou et.al Nature Communications 2018)


# Requirements: 
Most of these analyses can be run on any machine.
However, some of the data preparation steps, such as sequence annotation using IGoR (Marcou et.al Nature Communications 2018), are computationally intensive and require a cluster to run efficiently.
This sequence annotation script is written specifically for a cluster set up to use the Slurm job scheduler. 
(Some minor modifications to the [sequence annotation script](scripts/annotate_with_igor.sh) could allow this step to be run locally or using a different cluster workload manager. 


# About the analysis

With this analysis, we want to quantify the extent of microhomology involvement for nucleotide trimming and ligation during V(D)J recombination.
See the manuscript for specific model and methods details: 

TODO: add manuscript details

The following packages were especially helpful in our analyses:

Python:

* `jax`
* `jaxopt`
* `pandas`

R:

- `data.table` (Dowle and Srinivasan, 2021)
- `tidyverse` (Wickham et. al, 2019) 
- `doParallel` (Corporation and Steve Weston, 2020)
- `cowplot` (Wilke, 2020)
- `Biostrings` (Pagès et. al, 2021)

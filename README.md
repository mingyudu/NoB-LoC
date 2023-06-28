# Nonparametric Bayesian Local Clustering (NoB-LoC) Algorithm for Proteomics Dataset

This repository contains the implementation of the Nonparametric Bayesian Local Clustering (NoB-LoC) algorithm for proteomics datasets. The NoB-LoC algorithm is a Bayesian approach for clustering proteins and samples based on their expression profiles.

## Usage

To compile the C++ file before implementation, use the following command in the terminal (replace `file_name.cpp` with the name of the corresponding C++ file):

```bash
R CMD SHLIB file_name.cpp
```

The repository includes three major code files:

1. [**`wrapper_biclustering_mcmc.R`**](https://github.com/mingyudu/NoB-LoC/blob/main/wrapper_biclustering_mcmc.R): 
   - Function: `biclustering_mcmc(matrix, pi_0, pi_1, beta, M, n.MCMC.iter, burn.in, seed)`
   - Description: Use this function to implement Markov Chain Monte Carlo (MCMC) for posterior inference.
   - Parameters:
     - `matrix`: The input proteomics matrix with a dimension of sample number by protein number.
     - `pi_0`: Probability of a protein being in any active protein set.
     - `pi_1`: Probability of a sample being active given the protein of interest being in an active protein set.
     - `beta`: Total mass parameter of protein clustering. The number of protein sets increases as `beta` increases.
     - `M`: Total mass parameter of sample sub-clustering. The number of sample sets increases as `M` increases.
     - `n.MCMC.iter`: Number of MCMC iterations.
     - `burn.in`: Burn-in period for discarding initial MCMC samples.
     - `seed`: Random seed number for reproducibility.

2. [**`wrapper_biclustering_summarize.R`**](https://github.com/mingyudu/NoB-LoC/blob/main/wrapper_biclustering_summarize.R):
   - Function: `biclustering_summarize(save.result, burn.in, protein_id, sample_id, prot.clust.output, sample.clust.output)`
   - Description: Use this function to summarize the clustering results for proteins and samples.
   - Parameters:
     - `save.result`: Output from the function `biclustering_mcmc` including the MCMC inference result.
     - `burn.in`: Number of burn-in samples to discard.
     - `protein_id`: Identifier for proteins.
     - `sample_id`: Identifier for samples.
     - `prot.clust.output`: Output file name for protein clustering summary.
     - `sample.clust.output`: Output file name for sample clustering summary.

3. [**`wrapper_biclustering_heatmap.R`**](https://github.com/mingyudu/NoB-LoC/blob/main/wrapper_biclustering_heatmap.R):
   - Function: `plot_bicluster_heatmap(matrix, infer.result, prot_df, sample_df, scale='column', outfile)`
   - Description: Use this function to visualize the clustering results using heatmaps.
   - Parameters:
     - `matrix`: The input proteomics matrix with a dimension of sample number by protein number.
     - `infer.result`: Output from the function `biclustering_summarize` containing the estimated clustering membership for proteins and samples.
     - `prot_df`: Data frame containing protein annotation information.
     - `sample_df`: Data frame containing sample annotation information.
     - `scale`: The scaling method for the heatmap. The default is 'column', which scales the values based on the columns. Alternatively, setting it to 'none' means no scaling is applied.
     - `outfile`: Output file name for the heatmap visualization.

Please refer to the individual code files for detailed usage instructions and additional information about the input and output formats.

## Acknowledgement

We would like to express our gratitude to Dr. Juhee Lee at UCSC for providing the code for the NoB-LoC algorithm. If you are interested in the NoB-LoC algorithm, please refer to Dr. Juhee Lee's paper: [Link to the paper](https://www.tandfonline.com/doi/full/10.1080/01621459.2013.784705)

## License

This project is licensed under the [MIT License](LICENSE).
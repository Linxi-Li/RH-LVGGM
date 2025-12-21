Readme

R code for“Robust Heterogeneity Adjustment for Gaussian Graphical Model with Latent Variables”
**Author:** Linxi Li, Rong Li, Qingzhao Zhang, and Shuangge Ma


### Main files
1. main.R
2. sim.R
3. comparison.R
4. demo.R

For simulations, check the file `sim.R', which uses functions from `Generate.data()' and  `Metric()'  to setup experiments. `main.R' provides core model implementations (e.g., OT.RH.lvggm), while `comparison.R' includes competing methods. 
For a quick start, execute the  `demo.R ` file, which includes example data and model fitting workflows.

### Workflow Summary:

Data Generation: Use Generate.data() to create synthetic datasets with  sparse/low-rank structures.

Model Fitting: Run OT.RH.lvggm() to estimate cluster-specific parameters and detect noise points.

Performance Evaluation: Call Metric() to quantify estimation accuracy and clustering performance.

Comparison Methods: Use comparison.R to run baselines (e.g., LVGGM, Glasso) and compare results.

-----------------------------------------------------------------------------------------------------------------
### R Package Dependencies
install.packages(c("MASS", "glasso", "flexclust"))


### Key Function

Generate.data() Input:
@ n: sample size
@ p: number of variables
@ r: number of latent variables
@ pi0: proportion of noise points
@ piK: mixing proportions for each cluster
@ K: number of clusters
@ sparsity: sparsity level for the precision matrix
@ mue: magnitude of the mean differences between clusters
@ sparsity.latent: sparsity level for the latent variable connections
@ seed: random seed for reproducibility
Generate.data() Output:
A list including:
@ data: generated data matrix
@ S: true sparse component of the precision matrix
@ L: true low-rank component of the precision matrix
@ Theta: true precision matrix (S - L)
@ Mu: true mean vectors for each cluster
@ member: true cluster memberships
-----------------------------------------------------------------------------------------------------------------
OT.RH.lvggm() Input:
@ data: n * p data matrix
@ G: number of clusters
@ tau.ini: initial cluster membership probabilities
@ lambda1.seq: sequence of lambda1 values (sparsity penalty)
@ lambda2.seq: sequence of lambda2 values (low-rank penalty)
@ npr.max: maximum proportion of noise points
@ erc: eigen ratio constraint
@ iter.max: maximum number of iterations
@ tol: convergence tolerance
OT.RH.lvggm() Output:
A list including:
@ lambda.1: optimal lambda1 value
@ lambda.2: optimal lambda2 value
@ S: estimated sparse component
@ L: estimated low-rank component
@ Theta: estimated precision matrix
@ cluster: final cluster assignments
@ pi: estimated mixing proportions
@ mean: estimated cluster means
@ logicd: log-likelihood of the data
-----------------------------------------------------------------------------------------------------------------
Metric() Input :
@ data.test: Test dataset (n × p matrix)
@ memb.test: True cluster labels for test data (vector of length n)
@ lambda1: Optimal hyperparameter λ₁ from model fitting
@ lambda2: Optimal hyperparameter λ₂ from model fitting
@ fit: Fitted model object (output of OT.RH.lvggm())
@ S.true: True sparse matrix (p × p matrix)
@ Theta.true: True precision matrix (p × p matrix)
@ S.hat: Estimated sparse matrix (p × p matrix)
@ Theta.hat: Estimated precision matrix (p × p matrix)
@ L.hat: Estimated low-rank matrix (p × p matrix)
@ L.true: True low-rank matrix (p × p matrix)
@ memb.true: True cluster labels for training data (vector of length n)
@ memb.hat: Estimated cluster labels from the model (vector of length n)
Metric() Output:
Returns a data frame with the following metrics:
Structural Estimation Errors:
@ FL.S: Frobenius norm error for sparse matrix S (normalized by S.true)
@ FL.L: Frobenius norm error for low-rank matrix L (normalized by L.true)
@ FL.Theta: Frobenius norm error for precision matrix Theta (normalized by Theta.true)
@ rank: Rank of the estimated low-rank matrix L.hat
@ Sparsity Detection Metrics:
@ TPR: True Positive Rate (correct detection of non-zero entries)
@ FPR: False Positive Rate (incorrect detection of non-zero entries)
Clustering Performance:
@ ARI.c: Adjusted Rand Index for cluster labels (excluding noise points)
@ ARI.o: Adjusted Rand Index for outlier detection (noise points vs. clusters)
@ ARI.pre: Adjusted Rand Index for test data cluster predictions
-----------------------------------------------------------------------------------------------------------------
`comparison.R' implements several alternative methods, including:
LVGGM: Latent Variable Graphical Model (Homogeneous)
H-LVGGM: Heterogeneous LVGGM
OR-LVGGM: Oracle LVGGM (using true labels)
Glasso: Graphical Lasso  (Homogeneous)
OR-Glasso: Oracle Glasso
RH-GGM: Robust Heterogeneous GGM
H-GGM: Heterogeneous GGM
For further details, refer to the comments within each function and the associated publication. When evaluating the performance metrics of comparative methods, minor code adjustments may be required to avoid potential errors. For example, homogeneous models do not require  cluster prediction.
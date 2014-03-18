Topic Modeling using Latent Dirichlet Allocation
================================================

A Parallel Stochastic Collapsed Variational Bayesian Inference for LDA (SCVB0) implementation in C++ using <b>OpenMP</b>.

<p>We have implemented a parallel implementation of SCVB0 algorithm proposed by James Foulds et al. We refer <a href="http://dl.acm.org/citation.cfm?id=2487575.2487697">Stochastic collapsed variational Bayesian inference for latent Dirichlet allocation</a>. All the notations used are same as mentioned in this paper.</p>
<p>We have used <b>New York Times</b> dataset available at <a href="http://archive.ics.uci.edu/ml/datasets/Bag+of+Words">UCI Machine Learning Repository</a>. </p>
<p>The dataset is divided into minibatches of size 100 documents each.</p>
<p>We have used OpenMp to parallelize the execution of algorithm. All the minibatches are divided among the available number of processors and then algorithm is executed parallelly. Results are updated in global matrices nPhi, nTheta and nZ</p>
<p>We have analyzed the perplexity convergence on KOS and NIPS datasets available on the same webpage as that of NYT dataset.</p>
<p>Use following commands to execute the code</p>
<p>Compile:
<code>make</code></p>
<p>Execute:
<code>$ ./fastLDA docword.txt iterations NumOfTopics</code></p>

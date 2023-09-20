---
layout: post
title:  "What is the evolutionary transcriptome?"
date:   2023-09-18 11:18:49 +0200
categories: jekyll update
---
The transcriptome is the totality of mRNA expressed from the genome of an organism. The transcriptome, as do other phenotypes (if considered as the first phenotype of the genome) or biological processes, changes over successive generations (`evolves`). The `evolutionary transcriptome` emphasises and contextualises the transcriptome as an evolving property. Through the comparative method, one can dissect this layer of an organism's biology.

### Example analysis

Below, I give two examples (out of several) to study the evolutionary transcriptome. Tutorials for studying this property of life can be accessed through open-source bioinformatic software applications such as [`myTAI`](https://github.com/drostlab/myTAI) and [`scTEI`](https://github.com/kullrich/scTEI). 

`myTAI` seeks to quantify how enriched evolutionary ancient (i.e. older) or novel (i.e. younger) genes are between developmental stages, based on a key metric used in evolutionary transcriptomics referred to as the **transcriptome age index** (`TAI`). Stages that are more conserved have a lower TAI and stages that are less conserved have a higher TAI. Only a single organism's transcriptome is needed and species-specific genes are accounted for. 

#### Transcriptome age index

Running a TAI analysis on an example _Arabidopsis thaliana_ developmental transcriptome, we get

```r
devtools::install_github("drostlab/myTAI")
library(myTAI)
# example dataset covering 7 stages of A thaliana embryo development
data("PhyloExpressionSetExample")
PlotSignature(
  PhyloExpressionSetExample, 
  measure = "TAI", 
  ylab = "TAI", 
  permutations = 20000)
```

![TAI_example](https://github.com/LotharukpongJS/LotharukpongJS.github.io/assets/80110649/4ba1917c-0fc9-4729-8f9e-30d89f37d405)

Do try this! 

We see here that the torpedo stage has the lowest TAI (thus inferred to be the most conserved stage) and according to the [flat-line test](https://drostlab.github.io/myTAI/reference/FlatLineTest.html) (`p_flt`) run on 20,000 permutations, we see that the pattern here deviates significantly from a flat line.

While powerful detectors for evolutionary traces in the global (and single-cell) transcriptome, we have not exhausted the list of approaches to study the evolutionary transcriptome.

#### Transcriptome distance

For example, one can infer the correlation (e.g. Pearson) or distance (e.g. Manhattan distance) between the developmental transcriptomes of two organisms. Some other studies also employ more sophisticated approaches using `information theory`. These include distance metrics such as the **JSD metric** (which is the square root of the Jensen-Shannon Divergence since it satisfies the triangle inequality - see the [summary by Lior Patcher](https://liorpachter.wordpress.com/tag/jensen-shannon-divergence/)).

In this example, we use the R package [`philentropy`](https://drostlab.github.io/philentropy/index.html) and the example _Arabidopsis thaliana_ developmental transcriptome dataset from [`myTAI`](https://github.com/drostlab/myTAI) as above.

> Step 1: normalise the gene abundance matrix (e.g. raw TPM `TPM.mat`) to have a sum total of 1 for each sample. This is automatically done under scipy's [`scipy.spatial.distance.jensenshannon`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.jensenshannon.html) and given as an option in `philentropy`.
> 
> Step 2: using the package `philentropy`, apply the [`distance()`](https://drostlab.github.io/philentropy/reference/distance.html) function using the method `jensen-shannon`.
> 
> Step 3: square root the Jensen-Shannon divergence to get JSD metric. See [Österreicher & Vajda (2003)](https://link.springer.com/article/10.1007/BF02517812), [Endres & Schindelin (2003)](http://www.yaroslavvb.com/papers/endres-new.pdf) and [Fuglede & Topsoe (2004)](https://ieeexplore.ieee.org/abstract/document/1365067).

```r
# obtaining the required dataset
library(myTAI)
# example dataset covering 7 stages of A thaliana embryo development (same as above)
data("PhyloExpressionSetExample")
TPM.mat <- PhyloExpressionSetExample[,-c(1,2)]
rownames(TPM.mat) <- PhyloExpressionSetExample[,2]
TPM.mat <- log2(TPM.mat) # optionally transform the data to stabilise the variance in highly expressed genes.

# Step 1 & 2
# Step 1 is covered using est.prob = "empirical" in philentropy::distance
# Otherwise, one must do something like apply(TPM.mat, 2, function(x) x/sum(x))
TPM.JSDiv <- 
  philentropy::distance(
    t(TPM.mat), 
    use.row.names = TRUE, 
    method = "jensen-shannon", 
    est.prob = "empirical")

# Step 3
TPM.JSD <- sqrt(TPM.JSDiv)

# Plot
heatmap(TPM.JSD, Rowv = NA, Colv = NA, symm = T)
```

<img width="500" alt="Heatmap of Arabidopsis embryo JSD in R" src="https://github.com/LotharukpongJS/LotharukpongJS.github.io/assets/80110649/f3895ea5-4219-4713-b670-9e995018aacd">


and if we look at the raw output `round(TPM.JSD, 4)`,
```r
         Zygote Quadrant Globular  Heart Torpedo   Bent Mature
Zygote   0.0000   0.0103   0.0177 0.0177  0.0203 0.0225 0.0266
Quadrant 0.0103   0.0000   0.0143 0.0155  0.0180 0.0226 0.0303
Globular 0.0177   0.0143   0.0000 0.0080  0.0133 0.0178 0.0297
Heart    0.0177   0.0155   0.0080 0.0000  0.0092 0.0159 0.0289
Torpedo  0.0203   0.0180   0.0133 0.0092  0.0000 0.0126 0.0287
Bent     0.0225   0.0226   0.0178 0.0159  0.0126 0.0000 0.0246
Mature   0.0266   0.0303   0.0297 0.0289  0.0287 0.0246 0.0000
```

As you can see, the mature embryo is very different to other stages. As you can also notice, this example using samples from the _same species_ is not really interesting and only gives some clues about the evolutionary transcriptome between developmental stages. The cool stuff starts when comparing developmental transcriptomes _between species_.

You can see an example [here](https://www.nature.com/articles/s41586-018-0734-6/figures/3). But do you notice something funny about how they calculated the JSD?

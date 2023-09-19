---
layout: post
title:  "What is the evolutionary transcriptome?"
date:   2023-09-18 11:18:49 +0200
categories: jekyll update
---
The transcriptome is the totality of mRNA expressed from the genome of an organism. The transcriptome, as do other phenotypes (if considered as the first phenotype of the genome) or biological processes, changes over successive generations (`evolves`). The `evolutionary transcriptome` emphasises and contextualises the transcriptome as an evolving property. Through the comparative method, one can dissect how it evolves.

Tutorials for studying this property of life can be accessed through open-source bioinformatic software applications such as [`myTAI`](https://github.com/drostlab/myTAI) and [`scTEI`](https://github.com/kullrich/scTEI). Running a TAI analysis on an example _Arabidopsis thaliana_ developmental transcriptome, we get

{% highlight R %}
devtools::install_github("drostlab/myTAI")
library(myTAI)
# example dataset covering 7 stages of A thaliana embryo development
data("PhyloExpressionSetExample")
PlotSignature(PhyloExpressionSetExample, measure = "TAI", ylab = "TAI", permutations = 20000)
{% endhighlight %}

![TAI_example](https://github.com/LotharukpongJS/LotharukpongJS.github.io/assets/80110649/4ba1917c-0fc9-4729-8f9e-30d89f37d405)

Do try this! 

These approaches seek to quantify how enriched evolutionary ancient (i.e. older) or novel (i.e. younger) genes are between developmental stages, based on a key metric used in evolutionary transcriptomics referred to as the transcriptome age index (`TAI`). Stages that are more conserved have a lower TAI and stages that are less conserved have a higher TAI. Only a single organism's transcriptome is needed and species-specific genes are accounted for. 

We see here that the torpedo stage has the lowest TAI and according to the [flat-line test](https://drostlab.github.io/myTAI/reference/FlatLineTest.html) (`p_flt`) run on 20,000 permutations, we see that the pattern here deviates significantly from a flat line.

While powerful detectors for evolutionary traces in the global (and single-cell) transcriptome, we have not exhausted the list of approaches to study the evolutionary transcriptome.

For example, one can infer the correlation (e.g. Pearson) or distance (e.g. Manhattan distance) between the developmental transcriptomes of two organisms. Some other studies also employ more sophisticated approaches using `information theory`. For example, one can use the JSD metric (which is the square root of the Jensen-Shannon Divergence since it satisfies the triangle inequality - see the [summary by Lior Patcher](https://liorpachter.wordpress.com/tag/jensen-shannon-divergence/)).

To do this

Step 1: normalise the gene abundance matrix (e.g. raw TPM) to have a total probability of 1 for each sample.

{% highlight R %}
apply(merged_mat, 2, function(x) x/sum(x))
{% endhighlight %}

Step 2: apply philentropy’s distance function using “jensen-shannon”. Here, we use the package [`philentropy`](https://drostlab.github.io/philentropy/index.html).

{% highlight R %}
philentropy::distance(t(mat_normalized), use.row.names = TRUE, method = "jensen-shannon")
{% endhighlight %}

Step 3: square root the Jensen-Shannon divergence to get JSD metric. See [Österreicher & Vajda (2003)](https://link.springer.com/article/10.1007/BF02517812), [Endres & Schindelin (2003)](chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/http://www.yaroslavvb.com/papers/endres-new.pdf) and [Fuglede & Topsoe (2004)](https://ieeexplore.ieee.org/abstract/document/1365067).

You can see an example [here](https://www.nature.com/articles/s41586-018-0734-6/figures/3). But do you notice something wrong with how they calculated this? The distance values are rather big...

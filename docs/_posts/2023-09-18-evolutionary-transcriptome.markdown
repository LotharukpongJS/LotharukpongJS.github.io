---
layout: post
title:  "What is the evolutionary transcriptome?"
date:   2023-09-18 11:18:49 +0200
categories: jekyll update
---
The transcriptome is the totality of mRNA expressed from the genome of an organism. The transcriptome, as do other phenotypes (if considered as the first phenotype of the genome) or biological processes, changes over successive generations (`evolves`). The `evolutionary transcriptome` emphasises and contextualises the transcriptome as an evolving property. Through the comparative method, one can dissect how it evolves.

Tutorials for studying this property of life can be accessed through open-source bioinformatic software applications such as [`myTAI`](https://github.com/drostlab/myTAI) and [scTEI](https://github.com/kullrich/scTEI), i.e.

{% highlight R %}
devtools::install_github("drostlab/myTAI")
library(myTAI)
# example dataset covering 7 stages of A thaliana embryo development
data("PhyloExpressionSetExample")
PlotSignature(PhyloExpressionSetExample)
{% endhighlight %}

Try this! These approaches seek to quantify how enriched evolutionary ancient (i.e. older) or novel (i.e. younger) genes are between developmental stages, based on a key metric used in evolutionary transcriptomics referred to as the transcriptome age index (`TAI`). Stages that are more conserved has a TAI and stages that are less conserved has a higher TAI. Here, only a single organisms' transcriptome is needed and species-specific genes are accounted for. 

While powerful detectors for evolutionary traces in the global (and single-cell) transcriptome, we have not exhausted the list of approaches to study the evolutionary transcriptome.

For example, one can infer the correlation (e.g. pearson) between the developmental transcriptomes of two organism.

I will explain more some other day.

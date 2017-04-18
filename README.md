spacemap
================

Description
-----------

The spaceMap R package constructs *de novo* networks from heterogenous data types in a high-dimensional context by applying a novel conditional graphical model. The underlying statistical model is motivated by applications in integrative genomics, where two (or more) -omic data profiles are modeled jointly to discover their interactions. spaceMap is particularly effective in learning networks in diverse applications that exhibit hub topology. In an application to breast cancer tumor data, spaceMap learned regulatory networks that generated novel hypotheses of cancer drivers, while also identifying known signatures of various molecular sub-types. Additionally, an accompanying network analysis toolkit is provided, which helps interpret the complex biological interactions from the inferred networks.

Installation
------------

The R package `spacemap` is available from GitHub. Please see the [most recent beta release](https://github.com/topherconley/spacemap/releases/tag/v0.45.0-beta) for installation.

Details
-------

### Learning robust networks

We built model selection and model aggregation tools into spaceMap for learning networks robust to the challenges of integrative genomics data due to high-dimensionality and complex noise structure. High throughput -omic experiments often have small sample sizes and the signal-to-noise can be low. Models applied to the resulting data can be prone to overfitting and high variability. By a procedure called *CV.Vote* based on K-fold cross validation, one can learn model tuning parameters that will help balance the tradeoff of power and false discovery rate (FDR). *CV.Vote* reports a network where edges must be present in a majority of networks fit under the training sets. Please see [Model Tuning](https://topherconley.github.io/spacemap/articles/tune_sim1.html) under the Vignettes tab for an example of how to tune networks with cross validation. Specific usage details of *CV.Vote* are documented in the function [cvVote](https://topherconley.github.io/spacemap/reference/cvVote.html).

Making your network robust to overfitting can be further enhanced through bootstrap aggregation---called *Boot.Vote*---especially when sample size relative to dimensionality is a real concern. This procedure learns networks on *B* bootstrap replicates of data under CV-selected tuning parameters. Only those edges with majority representation among these networks are reported in the final network. Specific usage details and an example of *Boot.Vote* are documented in the functions [bootEnsemble](https://topherconley.github.io/spacemap/reference/bootEnsemble.html) and [bootVote](https://topherconley.github.io/spacemap/reference/bootVote.html).

#### Simulations

Several simulation studies have been used to evaluate spaceMap in learning networks with prominent hub topology---a commonly encountered feature of biological networks. The GitHub repository [sim-spacemap](https://github.com/topherconley/sim-spacemap) contains code affiliated with these simulations featuring:

-   simulation of hub network topology
-   data generation according to the network topology
-   fitting of spaceMap and other graphical models to the data
-   evaluation of graphical models results relative to true network topology

Documentation of these simulation codes will be forthcoming.

### Network analysis toolkit

Once a network has been learned, interpretation can prove quite challenging without a proper set of tools. We included a network analysis toolkit (see here for documentation) as part of the spaceMap package to facilitate network interpretation with special focus on integrative genomic applications. Our toolkit enables:

-   annotation of nodes (e.g. gene coordinates, functional description)
-   identification of cis/trans regulatory information
-   prioritization of hub nodes
-   module analysis (community detection)
-   pathway enrichment analysis (GO/KEGG)

All these features are reported through structured tables and are easily incorporated into technical reports in a variety of formats. Moreover, the network analysis tool kit integrates the results into a network file (e.g. `.graphml` format) that can be exported to existing tools such as the [Cytoscape ecosystem](http://www.cytoscape.org/what_is_cytoscape.html). In this sense, spaceMap is not just a model; rather it is a tool for deriving meaning from integrative genomics data.

This toolkit is illustrated through the analysis of The Breast Cancer Proteogomics Landscape Study data from TCGA, which is hosted on the GitHub repository [neta-bcpls](https://topherconley.github.io/neta-bcpls/).

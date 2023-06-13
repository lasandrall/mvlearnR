# mvlearnR
The package **mvlearnR** and accompanying Shiny App (https://multi-viewlearn.shinyapps.io/MultiView_Modeling/) is intended for integrating data from multiple sources
(e.g. genomics, proteomics, metabolomics). Most existing software packages for multiview learning are decentralized,
making it difficult for users to perform comprehensive integrative analysis. The new package wraps statistical and machine
learning methods and graphical tools, providing a convenient and easy data integration workflow. For users with limited
programming language, we provide a Shiny Application to facilitate data integration. The methods have potential to offer
deeper insights into complex disease mechanisms. Currently, **mvlearnR** can be used for the following:


**Prefiltering of each omics data via differential analysis (DA)**. We provide both supervised and unsupervised options for DA or for filtering out noise variables prior to performing data integration.

**Integrating data from two sources using a variant of the popular unsupervised method for associating data from two views, i.e. canonical correlation analysis (CCA)**.

**Predicting a clinical outcome using results from CCA**. We provide four outcome data distribution type (i.e. gaussian, binomial, Poisson, and time-to-event data.)

**Supervised integrative analysis (one-step) for jointly associating data from two or more sources and classifying an outcome. This method allows to include covariates**.

**Supervised integrative analysis (one-step) for jointly associating structured data from two or more sources and classifying an outcome. This method allows to include covariates** 

**Visualizatiing results from the DA or our integrative analysis methods**. These plots include: volcano plots, UMAP plots, variable importance plots, discriminant plots, correlation plots, relevance network plots, loadings plots, and within- and between- view biplots. These visualization tools will help unravel complex relationships in the data. 

**Demonstrating our integration workflow via already uploaded synthetic and real molecular and clinical data pertaining to COVID-19**.

Currently, linear multivariate methods for integrative analysis and biomarker identification are provided in **mvlearnR**. However, we have developed integrative analysis methods for disease subtyping 
(https://arxiv.org/abs/2111.06209) and nonlinear integrative analysis methods for biomarker identification (https://arxiv.org/abs/2302.07930; https://arxiv.org/abs/2111.09964; https://arxiv.org/abs/2304.04692) 
that will eventually be added to **mvlearnR** and  the accompanying web application. Other methods we develop in the future will be added. 
Thus, we envision **mvlearnR** and our web application to be a one-stop place for comprehensive data integration, for both users of R (or Python) and non-users of these software. 

If you have any questions or suggestions for improvements, please email: ssafo@umn.edu

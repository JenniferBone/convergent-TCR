# Convergent-TCR
A respository of examples of computational techniques for analyzing T cell receptors (TCRs) in the immune system

## TCR Informatics in Metastatic Renal Cancer
Evaluation of patients’ T and B cell repertoires is a potentially powerful avenue for diagnosis and therapeutic-response monitoring of immunotherapies. New methodologies for identifying and tracking cancer-specific CDR3 regions from the estimated 10<sup>15</sup> - 10<sup>25</sup> total diversity of TCR CDR3 sequences could allow for early prediction of cancer states and patient responses. However, development of these diagnostic tools is hampered by an unmet need for data-analytics tools capable of modeling small sample sizes. While LLMs gain a greater ubiquitous hold in the small molecule informatics realms, there is still a strong need to leverage statistical methods capable of handling high-feature data with low sample size. This work uses a combination of motif analysis, Markov modeling, and new word similarity metrics to featurize TCRs, and cluster them by similarity in order to reveal underlying disease-specific signal in TCR data.

**Figure 1** Tree-maps showing the large TCR sequence diversity of two RCC patients treated with high-dose IL-2 cancer immunotherapy. 

![Tree Maps showing sequence diversity two patient T cell repertoires](https://github.com/JenniferBone/convergent-TCR/blob/2236649463b61055bab2e2c29bf9512c168558d8/figures/Tree_maps.png)


## Improved Featurization


## This Example
Most of the code associated with this work is in a manuscript under review. However, the current repo seeks to demonstrate a small portion of code dedicated to improving TCR similarity analysis. It calculates a distance matrix between the top most frequent CDR3s  

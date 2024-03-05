# Convergent-TCR
A respository of examples of computational techniques for analyzing T cell receptors (TCRs) in the immune system

## TCR Informatics in Metastatic Renal Cancer
Evaluation of patients’ T and B cell repertoires is a potentially powerful avenue for diagnosis and therapeutic-response monitoring of immunotherapies. New methodologies for identifying and tracking cancer-specific CDR3 regions from the estimated 10<sup>15</sup> - 10<sup>25</sup> total diversity of TCR CDR3 sequences could allow for early prediction of cancer states and patient responses. However, development of these diagnostic tools is hampered by an unmet need for data-analytics tools capable of modeling small sample sizes. While LLMs gain a greater ubiquitous hold in the small molecule informatics realms, there is still a strong need to leverage statistical methods capable of handling high-feature data with low sample size. This work uses a combination of motif analysis, Markov modeling, and new word similarity metrics to featurize TCRs, and cluster them by similarity in order to reveal underlying disease-specific signal in TCR data.

**Figure 1** Tree-maps showing the large BCR (IGH) and TCR (delta-chain) sequence diversity of two RCC patients treated with high-dose IL-2 cancer immunotherapy. The patient who responded to immunotherapy is colored in blue. However the repertoire is generally undistinguishable for the non-responding patient in red

![Tree Maps showing sequence diversity two patient T cell repertoires](https://github.com/JenniferBone/convergent-TCR/blob/2236649463b61055bab2e2c29bf9512c168558d8/figures/Tree_maps.png)


## Improved Featurization
Transition probability maps based on a first-order Markov model of TCR pair-wise motifs was developed to reduce the feature size of TCR tree-maps to a 400-feature space. Jack-knifing was then used alongside thresholding to reveal underlying patterns in responder versus non-responder TCR data.

**Figure 2** NLP featurization of TCR data shows the reduction of large TCR sequence diversity revealing signature motifs (colored in shades of yellow) in the non-responder patient.

![TPM Featurizer of two patient T cell repertoires](https://github.com/JenniferBone/convergent-TCR/blob/29e2fd620f0db973a0218df3be0a28a4b5d97b12/figures/TCR_NLP_featurizer_small.png)

## This Example
Most of the code associated with this work is in a manuscript under review. However, the current repo seeks to demonstrate a small portion of code dedicated to improving TCR similarity analysis (namely improving on Levenstein distance for TCR similarity clusters). It calculates a distance matrix between the top most frequent CDR3s observed in indivual patient repertoires. 

Please contact: bone.jennifer@gmail.com for questions on code described in the repo.

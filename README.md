# Metabolic Engineering for DPAHelix project

This repository corresponds to one of the five computational biology analyses of the DPAHelix project in the GOGEC Competition of 2026.

> **Objective**: Optimize the production of DPA through metabolic engineering.

We worked our model using notebooks from Google Colab (https://colab.research.google.com/). Important updates to our code were later added to this repository.

_(The following activities may not be ordered exactly as described in the final written report delivered in the competition)_

## Key resources

TODO: Add WSL and other stuff info

Databases:

1. KEGG: https://www.kegg.jp/
2. BioCyc: https://biocyc.org/
3. EMBL-EBI BioModels: https://www.ebi.ac.uk/biomodels/
4. BiGG: http://bigg.ucsd.edu/

Python libraries:

1. COBRApy: https://cobrapy.readthedocs.io/en/latest/index.html

We selected the following model:

> Henry2009 - Genome-scale metabolic network of Bacillus subtilis (iBsu1103): https://www.ebi.ac.uk/biomodels/MODEL1507180015#Files

The model provides:

1. 1,437 reactions associated with 1,103 genes.
2. Includes Gibbs free energy values for 1,403 reactions.
3. 653 irreversible reactions.
4. validation against an experimental dataset of 1,500 distinct conditions (Model accuracy: 89.7 - 93.1%).

We downloaded the model using:

```shell
wget "https://www.ebi.ac.uk/biomodels/services/download/get-files/MODEL1507180015/2/MODEL1507180015_url.xml" -O "model/MODEL1507180015_url.xml"
"
```

## Activity 1: Identification of key reactions

Prior any modeling analysis, we selecting our key reactions using KEGG.

## Activity 1: Simulation of precursor production under normal conditions

We analyzed the production of x (most important precursor of DPA biosynthesis) under normal conditions in _B. subtilis_.

We downloaded the model...

We set our key metabolic path to produce PHA...

We set conditions...

(TODO: Finish writing)

## Activity 2

_(Under development...)_
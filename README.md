# Project description 

In RNA sequencing experiments, the expression of mRNA is measured in cells or tissues using an experimental workflow that involves a number of enzymatic conversion steps. 

First, RNA molecules are caputured and converted into complementary DNA (cDNA). The pool of cDNA molecules is referred to as a "cDNA library" and typically includes millions of molecules. Conversion to cDNA is essential as DNA is much easier to process in the lab. 

Then, to generate sufficient concentrations of the "cDNA library" to make it measurable with sequencing instruments, one need to amplifiy the "cDNA library" exponentially by a reaction called PCR. This causes the cDNA concentration to increase by several orders of magnitude.

However, because of this exponential amplification step, the data generated contain a high fraction of duplicated sequences. The PCR reaction is known to be biased towards certain types of molecules, so it becomes important to be able to remove "PCR duplicates". 

To do this, a so called Unique Molecular Identifier (UMI) can be coupled to each RNA molecule when it's captured. The UMI becomes part of the cDNA molecule and is therefore also amplified. This making it possible to trace PCR copies to one single RNA molecule, remove any redundant copies and quantify the original number of RNA molecules. 

In this project, I will present a small simulation framework to generate cDNA libraries in silico and benchmark a set of  "duplicate removal algorithms" to see how well they perform in removing PCR duplicates.

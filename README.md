# HBNDmodel

## Introduction  
Hydrogen bond network dimensionality (HBND) describes how hydrogen-bond intermolecular interactions extend in a three-dimensional structure. The network arrangement is thought to influence properties where hydrogen bonding plays critical  role, such as crystal stability, mechanical  behaviour  and compound tabletability performance.  

The hydrogen bond network dimensionality problem can be formulated as a four-class classification task. Intermolecular hydrogen interactions between molecules may lead to enclosed motifs (0D), chains (1D), sheets (2D) or expansion of the network in all directions (3D).  

In this page we provide a hbnd dataset and a hbnd predictive model.   

## The dataset  

**Data retrieval**  
The dataset was obtained from the Cambridge Structural Database (CSD). The CSD was searched for all **organic** crystal structures of a **single chemical component**. Entries containing:  
- metals, salts, or ions were discarded, as these present additional challenges that were not addressed in the study.   
- disorder, errors or incomplete information about crystal atomic coordinates or hydrogen bonds were discarded, as they would not provide enough information for accurate hbnd label calculation.   

Of these, molecules with more than one crystal structure submitted to the database were removed. This step removes conflicting data where the compound is polymorphic and its different crystal arrangements are reported to have different network dimensionalities. 

**Label Assignment**   
Dimensionality was calculated through the computation of the square roots of the eigenvalues of the covariance matrix of the atomic coordinates of the supramolecular structures that resulted from two different expansions of the network, through hydrogen bond intermolecular interactions, using methods from the CSD Python API. Ratios for each dimension before and after the expansion were calculated, to deduce the number of directions in which the network grew.  

**Final Dataset**   
The final dataset includes 63 839 labelled examples across 4 classes: 0D (22971), 1D (30718), 2D (7013), 3D (3137). This can be found in the hbnd_dataset directory, under the name ***hbnd_dataset***. You may access the content with:  

	$ import pandas as pd
	$ hbnd = pd.read_pickle('hbnd_dataset')
	
The dataset contains 2 columns:   
- *id*, containing the molecule SMILES string  
- *label*, containing the correspondent hbnd tag  

## The model

**Feature calculation**  
The hbnd model was produced from the dataset above. One hundred and fifteen descriptors of two
and lower dimensions were calculated for each molecule using the RDkit package. 

**Data processing / Model fitting**  
Data processing included cleansing, scaling and class size balancing.   
The dataset was divided into a training (80%) and test set (20%), and different methods were considered for model fitting and hyperparameter optimisation.

The best performing model was produced by an SVM RBF routine, which we provide here:  
- **Input:** list of SMILES of the molecules to be predicted  
- **Output:** a list of (hbnd label, SMILES) pairs  
- **Performance:** 59% accuracy (compared to a 25% random threshold)  

This can be found in the hbnd_model directory, under the name ***hbnd_model***. For model usage, check the basic_tour walkthrough.  

See the [reference](https://pubs.rsc.org/en/content/articlelanding/2020/ce/d0ce00111b#!divAbstract) and basic_tour for more details. 

## Dependencies   
The code should be run using Python 3. It also requires the RDkit, Pandas, NumPy, and JobLib modules. Dependency installation via conda:  

      $ conda install pandas numpy pickle joblib matplotlib scikit-learn
      $ conda install -c conda-forge rdkit

## Usage 

**The model may be used as shown in the *walkthrough.ipynb* in the hbnd_model directory.**  

**Note**  
The model is only suitable to predict the HBND of **single chemical component organic crystals**.   

Using the model as is may only be suitable for initial initial screening of large libraries. The trust of output predictions can be increased by subjecting the model to a [confidence threshold](https://github.com/apfrade/ConfidenceMeasure/blob/master/README.md). This approach may enable further single compound HBND profilling.  

See the [reference](https://pubs.rsc.org/en/content/articlelanding/2020/ce/d0ce00111b#!divAbstract) for more details. 
 
	  
# References   

**Please cite us:**  

*A. P. Frade, P. McCabe and R. I. Cooper. “Increasing the performance, trustworthiness and practical value of machine learning models: a case study predicting hydrogen bond network dimensionalities from molecular diagrams”. 2020. CrystEngComm. DOI: 10.1039/D0CE00111B* 

Access the paper [here](https://pubs.rsc.org/en/content/articlelanding/2020/ce/d0ce00111b#!divAbstract).

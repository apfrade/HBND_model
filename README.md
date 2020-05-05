# HBNDmodel

## Introduction  
Hydrogen bond network dimensionality (HBND) describes how hydrogen-bond intermolecular interactions extend in a three-dimensional structure. The network arrangement is thought to influence properties where hydrogen bonding plays critical  role, such as crystal stability, mechanical  behaviour  and compound tabletability performance.  

The hydrogen bond network dimensionality problem can be formulated as a four-class classification task, and the four possible network dimensionality outcomes are schematically represented below. We provide a ready to use model for Hydrogen Bond Network Dimensionality prediction.  

- **Input:** list of SMILES of the molecules to be predicted  
- **Output:** a list of (hbnd label, SMILES) pairs  
- **Performance:** 59% accuracy (compared to a 25% random threshold)  

See the [reference](https://pubs.rsc.org/en/content/articlelanding/2020/ce/d0ce00111b#!divAbstract) for more details.  

## Disclamer   
1) The model is only prepared to predict the HBND of single chemical component organic crystals.  
2) Using the model as is may only be suitable for initial initial screening of large libraries. The trust of output predictions can be increased by subjecting the model to a [confidence threshold](https://github.com/apfrade/ConfidenceMeasure/blob/master/README.md). This approach may enable further single compound HBND profilling.  

See the [reference](https://pubs.rsc.org/en/content/articlelanding/2020/ce/d0ce00111b#!divAbstract) for more details. 

## Dependencies   
The code should be run using Python 3.  

     RDkit
     Pandas
     NumPy
     JobLib
   
Dependency installation via conda:  

      $ conda install pandas numpy pickle joblib matplotlib scikit-learn
      $ conda install -c conda-forge rdkit

## Installation  
Via pip:

      $ pip install -i https://test.pypi.org/simple/ HBNDmodel  
	  
# References   

**Please cite us:**  

*A. P. Frade, P. McCabe and R. I. Cooper. “Increasing the performance, trustworthiness and practical value of machine learning models: a case study predicting hydrogen bond network dimensionalities from molecular diagrams”. 2020. CrystEngComm. DOI: 10.1039/D0CE00111B* 

Access the paper [here](https://pubs.rsc.org/en/content/articlelanding/2020/ce/d0ce00111b#!divAbstract).

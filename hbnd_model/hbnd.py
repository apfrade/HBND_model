try:
    import os
    import numpy as np
    import pandas as pd
    from joblib import load
    from rdkit import Chem
    from rdkit.ML.Descriptors import MoleculeDescriptors 
except ImportError:
    raise ImportError(
        "In order to run the model you must have the libraries: " +
        "`numpy`, `pandas`, 'joblib' and 'rdkit' installed.") 

class PredictiveModel:

    def _get_descriptors(self, smiles_list):
        
        '''Calculates the descriptor list of a molecule. 
           It uses the rdkit package.
        '''
        
        descriptor_names = load('descriptor_names')
        desc_object = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
        
        invalid_molecules = []
        ids, descriptors = [], []
        
        for smiles in smiles_list:
        
            try:
                mol = Chem.MolFromSmiles(smiles)
                desc_list = list(desc_object.CalcDescriptors(mol))
                descriptors.append(desc_list)
                ids.append(smiles)
            
            except:
                invalid_molecules.append(smiles)
                
        if len(invalid_molecules) > 0:
            print('Some molecules could not be processed.')
            print(invalid_molecules)        
                
        if len(descriptors)==0:
            print('No molecules could be processed.')
            
        else: 
            df = pd.DataFrame(descriptors, columns=descriptor_names)
            df.insert(0, 'smiles', ids)
            return df
            

    def _predict_hbnd(self, df):
        
        ''' A function that predicts HBND.'''
        
        # partition data
        ids = df['smiles'].values
        descs = df.drop(['smiles'], axis=1)
        
        # load hbnd model
        hbnd_model =  load('hbnd_model')
        
        # predict
        hbnd = hbnd_model.predict(descs)
        
        # report
        res = []
        for molecule, dimensionality in zip(ids, hbnd):
            res.append((dimensionality, molecule))
        
        return res

    def hbnd_model(self, smiles_list,):
        
        '''Predicts the HBND of molecules given their SMILES.
        
        Parameters
        ----------
        smiles: 1d array
            Array of SMILES of the molecules to be predicted.
       
        Returns
        -------
        prediction_list: 1d list
            The list of tuples (SMILES, predicted HBND) containing the results of the prediction. 
        '''
        
        descriptors = self._get_descriptors(smiles_list)
        prediction_list = self._predict_hbnd(descriptors)
        
        return prediction_list

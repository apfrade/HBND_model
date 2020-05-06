try:
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
        
        meta = load('meta')
        desc_object = MoleculeDescriptors.MolecularDescriptorCalculator(meta['descriptor_names'].values)

        invalid_molecules = []
        ids, descriptors = [], []
        for smiles in smiles_list:

            try:
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    desc_list = list(desc_object.CalcDescriptors(mol))
                    descriptors.append(desc_list)
                    ids.append(smiles)
            except:
                invalid_molecules.append(smiles)
                
        if len(invalid_molecules) > 0 and len(descriptors)!=0:
            print('Some molecules could not be processed.')
            print(invalid_molecules)        
                
        if len(descriptors)==0:
            print('No molecules could be processed.')
            
        else: 
            df = pd.DataFrame(descriptors, columns=meta['descriptor_names'].values)
            df.insert(0, 'smiles', ids)
            return df 
            
            
    def _preprocess(self, df):

        meta = load('meta')
        scaled = meta['scaler'].transform(df.drop(['smiles'], axis=1))
        scaled = pd.DataFrame(scaled, columns=df.columns.values[1:])
        scaled.insert(0, 'smiles', df.smiles.values)
        
        scaled.fillna(value=np.nan, inplace=True)
        scaled.replace([np.inf, -np.inf], np.nan, inplace=True)
        scaled.dropna(inplace=True)
        scaled.reset_index(inplace=True, drop=True)
        
        return scaled
    
    
    def _predict_hbnd(self, df):

        ids = df['smiles'].values
        descs = df.drop(['smiles'], axis=1)
        hbnd_model =  load('hbnd_model')
        
        hbnd = hbnd_model.predict(descs)
        res = [(label, molecule) for molecule, label in zip(ids, hbnd)]
        
        return res

    def hbnd_model(self, smiles_list,):
        
        '''Predicts the HBND of molecules given their SMILES.
        
        :param 
        smiles_list: 1d array containing the SMILES of the molecules to be predicted.
       
        :return 
        predictions: 1d list of (predicted label, SMILES) tuples containing the results. 
        '''
        
        descriptors = self._get_descriptors(smiles_list)
        descriptors = self._preprocess(descriptors)
        predictions = self._predict_hbnd(descriptors)
        return predictions

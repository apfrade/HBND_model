try:
    import pickle
    import numpy as np
    import pandas as pd
    from joblib import load
    from rdkit import Chem
    from rdkit.ML.Descriptors import MoleculeDescriptors 
except ImportError:
    raise ImportError(
        "In order to run the model you must have the libraries: " +
        "`numpy`, `pandas`, 'joblib' and 'rdkit' installed.")    


def clean_df(df):
    
    ''' A function to clean rows of a dataframe
    
    It eliminates any rows containing None, inf or nan values. 
    
    Parameters
    ----------
    df: pandas dataframe
        The original dataframe.
    
    Returns
    -------
    df: pandas dataframe
        The cleaned version of the original df.
    '''
    
    df.fillna(value=np.nan, inplace=True)
    df.replace([np.inf, -np.inf], np.nan, inplace=True)
    df.dropna(inplace=True)
    df.reset_index(inplace=True, drop=True)
    
    return df

def get_mols(molecules):
    
    ''' A function that retrieves the mol representation of a molecule
    
    It uses the rdkit package to convert a list of molecules from their
    SMILES to their mol representation.
    
    Parameters
    ----------
    molecules: 1d array
        Array of SMILES of the molecules to be predicted.
    
    Returns
    -------
    mols: pandas dataframe
        The clean dataframe containing the SMILES and mol representation of each molecule.
    '''
    
    print('')
    
    mols = []
    for smiles in molecules:
        
        mol = Chem.MolFromSmiles(smiles)
        mols.append(mol)
        
        if mol is None: print ('Mol could not be calculated for %s' %smiles)
            
    mols = pd.DataFrame({'smiles': molecules, 'mols': mols})
    mols = clean_df(mols)
    
    return mols

def get_descriptors(mols):
    
    '''A function that retrieves the descriptor list of a molecule 
    
    It uses the rdkit package to calculate a list of numerical descriptors from the mol
    representation of each molecule. The list of descriptors to be calculated is given 
    by the file descriptor_names.
    
    Parameters
    ----------
    mols: pandas dataframe
        The clean dataframe containing the SMILES and mol representation of each molecule.
        
    Returns
    -------
    df: pandas dataframe
        The clean dataframe containing the SMILES representation and calculated descriptors for each molecule. 
    '''
    
    print('')
    
    # get list of descriptors
    descriptor_names = pd.read_pickle('data/descriptor_names')
    
    # retrive descriptors
    desc_object = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
    
    descriptors = []
    for mol in mols.mols.values:
        
        desc_list = list(desc_object.CalcDescriptors(mol))
        descriptors.append(desc_list)
        
        if None in desc_list: print ('Descriptors could not be calculated for %s' %df.smiles[df.mols == mol])
            
    df = pd.DataFrame(descriptors, columns=descriptor_names)
    df.insert(0, 'smiles', mols.smiles.values)
    df = clean_df(df)
    
    return df

def scale_df(df):

    '''A function that scales the desciptor values. 
    
    It scales the values of descriptors according to the scaler model stored in the file
    scaler_model.bin.
    
    Parameters
    ----------
    df: pandas dataframe
        The clean dataframe containing the SMILES representation and calculated descriptors for each molecule. 
       
    Returns
    -------
    scaled_df: pandas dataframe
        The scaled version of df containing SMILES and descriptors.
    '''
    
    # load scaling model
    scaler_model = load('data/scaler_model.bin')
    
    # scale descs
    X = df.drop(['smiles'], axis=1)
    X_test = scaler_model.transform(X)
    
    # rebuild df
    scaled_df = pd.DataFrame(X_test, columns=X.columns.values)
    scaled_df.insert(0, 'smiles', df.smiles.values)
    
    return scaled_df

def predict(scaled_df, report=False):
    
    ''' A function that predicts HBND.

    It predicts the HBND of a molecule from its correspondent set of
    scaled descriptors. The predictive model is stored in the file
    hbnd_model.sav.
    
    Parameters
    ----------
    scaled_df: pandas dataframe
        The scaled version of df containing SMILES and descriptors.
    
    report: bool, default=False
        Whether to print to the console the molecule SMILES and predicted HBND.
    
    Returns
    -------
    res: 1d list
        The list of tuples (SMILES, predicted HBND) containing the results of the prediction. 
    '''
    
    print('')
    
    # partition data
    ids = scaled_df['smiles'].values
    descs = scaled_df.drop(['smiles'], axis=1)
    
    # load hbnd model
    hbnd_model =  pickle.load(open('data/hbnd_model.sav', "rb"))
    
    # predict
    hbnd = hbnd_model.predict(descs)
    
    # report
    res = []
    if report: print('PREDICTIONS:\n\nHBND    SMILES' )
    for molecule, dimensionality in zip(ids, hbnd):
        
        res.append((dimensionality, molecule))
        if report: print(' %s  |  %s' %(dimensionality, molecule))
    
    return res

def hbnd_predictive_model(smiles, report=False):
    
    '''Predicts the HBND of molecules given their SMILES.
    
    Parameters
    ----------
    smiles: 1d array
        Array of SMILES of the molecules to be predicted.
    
    report: bool, default=False
        Whether to print to the console the molecule SMILES and predicted HBND.
    
    Returns
    -------
    prediction_list: 1d list
        The list of tuples (SMILES, predicted HBND) containing the results of the prediction. 
    '''
    
    mols = get_mols(smiles)
    descriptors = get_descriptors(mols)
    scaled_descriptors = scale_df(descriptors)
    prediction_list = predict(scaled_descriptors, report=report)
    
    return prediction_list

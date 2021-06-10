import os
# import sys
# import re
import json
import warnings
import pandas as pd
import numpy as np
from glob import glob
from tqdm.auto import tqdm
from multiprocessing import Pool, process
from time import time
from tqdm.contrib.concurrent import process_map 
from Bio.PDB import PDBParser, parse_pdb_header
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from IPython.display import display

PARAMS_FP='../config/config.json'

def find_pdb_files(src):
    '''
    Find al .pdb files in directory
    ----------
    Parameters
        src - str
            Path to directory containing pdb files. 
            Path expected to be of form  <src>/**/*.pdb
    -------
    Returns
    -------
        list of all .pdb filepaths in <src>
    '''
    filepath=os.path.join(src, '**/*.pdb')
    files=glob(filepath)
    return files

def get_protein_data_multiprocess(fn):
    '''
    Get protein structure's data for each atom_id, atom_name, model_id, chain_id, residue_name, 
    and atomic (x,y,z) coordinates. Currently runs slower than single processed version.
    ----------
    Parameters
    ----------  
        fn - str
            filepath to .pdb file to get protein structure object from
    -------
    Returns
    -------
        2d numpy array of shape (K, 9) array, where K is the number of atoms in the entire structure.
        Columns are (protein_id, model_id, chain_id, residue_id, residue_name, atom_name, atom_coord_x, atom_coord_y, atom_coord_z)
    '''
    def parse_pdb(fn):
        '''
        Get structure object from file
        ----------
        Parameters
        ----------
            fn - str
                filepath to .pdb file to get protein structure object from
        -------
        Returns
        -------
            Bio.PDB.Entity.Entity structure object from PDB file <fn>
        '''
        structure_name=os.path.basename(fn).split(".")[0]
        # print('structure_name',structure_name)

        structure=parser.get_structure(
            id=structure_name,
            file=fn
        )

        return structure

    parser=PDBParser()
    structure=parse_pdb(fn)
    structure_arr=np.array([])

    for model in structure.get_models():
        for chain in model.get_chains():
            for residue in chain.get_residues():
                for atom in residue.get_atoms():
                    structure_id, model_id, chain_id, residue_id, atom_name=atom.get_full_id()
                    data=[structure_id, model_id, chain_id, residue_id[1], atom_name[0]]
                    data=data+list(atom.get_coord())
                    structure_arr=np.append(structure_arr, data)

    num_cols=len(data)
    num_rows=structure_arr.shape[0]//num_cols                
    structure_arr=structure_arr.reshape(num_rows,num_cols)
    return structure_arr.tolist()

def get_protein_data(structure):
    '''
    Get protein structure's data for each atom_id, atom_name, model_id, chain_id, residue_name, 
    and atomic (x,y,z) coordinates.
    ----------
    Parameters
    ----------  
        structure - Bio.PDB.Structure object
            PDB structure object to get molecular data from
    -------
    Returns
    -------
        2d numpy array of shape (K, 9) array, where K is the number of atoms in the entire structure.
        Columns are (protein_id, model_id, chain_id, residue_id, residue_name, atom_name, atom_coord_x, atom_coord_y, atom_coord_z)
    '''
    num_cols=8
    num_atoms=pd.Series(structure.get_atoms()).shape[0]
    out_shape=(num_atoms, num_cols)

    structure_arr=list(np.empty(out_shape))
    # structure_arr=np.array([])
    
    row=0
    for atom in structure.get_atoms():
        structure_id, model_id, chain_id, residue_id, atom_name=atom.get_full_id()
        data=[structure_id, model_id, chain_id, residue_id[1], atom_name[0]]
        data=data+list(atom.get_coord())
        # structure_arr=np.append(structure_arr, data)
        structure_arr[row]=np.array(data)
        row+=1

    # num_cols=len(data)
    # num_rows=structure_arr.shape[0]//num_cols                
    # structure_arr=structure_arr.reshape(num_rows,num_cols)

    return np.array(structure_arr)

def extract_pdb_files_multiprocess(src, dst=None):
    '''
    Extract .pdb files into dataframe. Currently runs slower than single processed version
    ----------
    Parameters
    ----------
        src - str
            Path to directory containing pdb files. 
            Path expected to be of form  <src>/**/*.pdb
        dst - str, Default: None
            Path to save extracted .pdb files DataFrame to.
            If none do not save DataFrame
    -------
    Returns
    -------
        Dataframe of extracted .pdb file data
    ''' 
    files=find_pdb_files(src)[:5]

    structures=np.array(process_map(get_protein_data_multiprocess, files))


    structures_arr=np.array([])
    for structure in structures:
        structures_arr=np.append(structures_arr, structure)

    num_cols=8
    num_rows=structures_arr.shape[0]//num_cols                
    structures_arr=structures_arr.reshape(num_rows,num_cols)

    structure_cols=[
        'protein_id',
        'model_id',
        'chain_id',
        'residue_id',
        # 'residue_name',
        'atom_name',
        'atom_coord_x',
        'atom_coord_y',
        'atom_coord_z'
    ]

    structures=pd.DataFrame(
        data=structures_arr,
        columns=structure_cols
    )

    return structures

def extract_pdb_files(src, dst=None):
    '''
    Extract .pdb files into dataframe
    ----------
    Parameters
    ----------
        src - str
            Path to directory containing pdb files. 
            Path expected to be of form  <src>/**/*.pdb
        dst - str, Default: None
            Path to save extracted .pdb files to.
            If none do not save DataFrames
    -------
    Returns
    -------
        Dataframe of extracted .pdb file data if <dst> is None. Otherwise return None.
    ''' 
    start_time=time()
    def parse_pdb(fn):
        '''
        Get structure object from file
        ----------
        Parameters
        ----------
            fn - str
                filepath to .pdb file to get protein structure object from
        -------
        Returns
        -------
            Bio.PDB.Entity.Entity structure object from PDB file <fn>
        '''
        structure_name=os.path.basename(fn).split(".")[0]
        structure=parser.get_structure(
            id=structure_name,
            file=fn
        )
        return structure

    files=find_pdb_files(src)#[:10]

    #intilize parser once
    parser=PDBParser()
    iteration=0

    structure_cols=[
        'protein_id',
        'model_id',
        'chain_id',
        'residue_id',
        'atom_name',
        'atom_coord_x',
        'atom_coord_y',
        'atom_coord_z'
    ]

    if dst:
        for file in tqdm(files, total=len(files)):
            structure=parse_pdb(file)
            structure=get_protein_data(structure)
            df=pd.DataFrame(
                data=structure,
                columns=structure_cols
            )
            save_fn=os.path.basename(file).replace('.pdb','.csv')
            save_path=os.path.join(dst, os.path.basename(os.path.dirname(file)))
            os.makedirs(save_path, exist_ok=True)
            df.to_csv(os.path.join(save_path, save_fn), index=False)
    else:
        for file in tqdm(files, total=len(files)):
            structure=parse_pdb(file)
            if iteration==0:
                structures=get_protein_data(structure)
            else:
                structure_arr=get_protein_data(structure)
                structures=np.append(structures, structure_arr, axis=0)
            iteration+=1

        structures=pd.DataFrame(
            data=structures,
            columns=structure_cols
        )
        
        print('Completed data cleaning in %s seconds'%int(time()-start_time))

        return structures

if __name__=='__main__':
    with open(PARAMS_FP, "r") as file:
        params=json.load(file)

    params=params['cleaning-params']
    src=params['src']
    dst=params['dst']

    #hide discontinuous chain warnings
    warnings.simplefilter('ignore', PDBConstructionWarning)
    extract_pdb_files(src, dst)
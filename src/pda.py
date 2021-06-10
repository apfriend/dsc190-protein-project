import os
import sys
import json
import resource
import warnings
import gudhi as gd
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import argparse
from mpl_toolkits import mplot3d
from tqdm.auto import tqdm
from glob import glob
from scipy.spatial.distance import pdist, squareform
from Bio.PDB.PDBExceptions import PDBConstructionWarning

from Bio.PDB import PDBParser

def get_size_of(fn):
    '''
    Get size of file 
    ----------
    Parameters
        fn - str
            filepath to file
    -------
    Returns
    -------
        size of file in bytes
    '''
    st=os.stat(fn)
    return st.st_size

def limit_memory():
    '''limit memory to 18 gb used to catch out of memory error'''
    maxsize=18*(2**30)
    soft, hard=resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (maxsize, hard))
    print('Memory Limited to %i bytes'%maxsize)

def rips_complex(data, r_max=None):
    '''
    Compute Rips Complex for <data>
    ----------
    Parameters
    ----------
        data - list-like object
            List of molecular relative coordinates. 
            Expected dimension of (K, 3), where each column represents 
            a position in 3d space for K atoms.
        radius - float (optional)
            Maximum radius to compute Rips Complex from.
            If not passed then 
    -------
    Returns
    -------
        gudhi rips complex object
    '''
    if r_max is None:
        r_max=np.max(pdist(data))
    
    rips_complex=gd.RipsComplex(
        points=data,
        max_edge_length=r_max
    )

    return rips_complex

def persistence(csv_src, pdb_src, dst, r_max_list, dimensions=2):
    '''
    Compute persistence profile from rips comlex. 
    If persistence profile takes too much memory then a MemoryError
    will be raised and a value of None will be returned for both `birth_death`
    and `simplex_tree` 
    ----------
    Parameters
    ----------

        dimension - int (optional)
            Number of dimenions to compute persistence diagram of.
            If not passed thend defaults to 2.
    -------
    Returns
    -------
        (birth_death, simplex_tree)
            birth_death - features from rips complext simplex tree with their birth and death 
                          as a tuple in the second element
            simplex_tree - simplex tree from which features were created
    '''
    csv_files=glob(os.path.join(csv_src, '*', '*.csv'))
    pdb_files=glob(os.path.join(pdb_src, '*', '*.pdb'))
    # print(pdb_files)
    struture_keys=[os.path.basename(pdb).replace('.pdb','') for pdb in pdb_files]
    pdb_file_lookup=dict(zip(struture_keys, pdb_files))
    # print(pdb_file_lookup)
    parser=PDBParser()

    num_iter=len(r_max_list)*len(csv_files)
    # df_list=[None]*num_iter
    data=[None]*num_iter

    idx=0
    structure_tag=os.path.basename(csv_files[0]).replace('.csv','')
    for i, csv_file in tqdm(enumerate(csv_files), total=len(csv_files)):
        structure_tag=os.path.basename(csv_file).replace('.csv','')
        pdb_file=pdb_file_lookup[structure_tag]
        structure=parser.get_structure(structure_tag, pdb_file)
        num_chains=len(list(structure.get_chains()))
        df=pd.read_csv(csv_file)
        point_data=df.loc[:,['atom_coord_x','atom_coord_y','atom_coord_z']].values

        for r_max in tqdm(r_max_list, desc=structure_tag, total=len(r_max_list), leave=False): 
            simplicial_complex=rips_complex(point_data, r_max)
            simplex=simplicial_complex.create_simplex_tree(max_dimension=dimensions)
            simplex.compute_persistence()

            # fn=os.path.join(dst, 'radius-%i'%r_max,'%s.pyx'%structure_tag)
            # print(fn)
            # simplex.write_persistence_diagram(fn)
            betti=simplex.betti_numbers()
            betti_0=betti[0]
            if len(betti)>1:
                betti_1=betti[1]
            else:
                betti_1=0
           
            atoms=df.atom_name.unique()
            num_atoms=df.shape[0]

            structure_data=[
                structure_tag,
                num_atoms,
                num_chains, 
                betti_0, 
                betti_1, 
                r_max, 
                simplex.num_vertices(), 
                simplex.num_simplices()
            ]
            data[idx]=structure_data           
            idx+=1

    os.makedirs(dst, exist_ok=True)
    df2_dst=os.path.join(dst, 'birth_death_data.csv')
    col_names=[
        'protein_tag',
        'num_atoms',
        'num_chains',
        'betti_0',
        'betti_1',
        'max_radius',
        'num_vertices',
        'num_simplices'
    ]

    df2=pd.DataFrame(
        data=data,
        columns=col_names
    )
    df2.to_csv(df2_dst, index=False)
    print("done")

def get_protein_tag_table():
    '''
    Returns a dictionary of protein tags as keys, with a human readable name of
    what they are. Ex {'1trz':'Designer-Insulin}
    '''
    files=glob('../data/test/csv-data/*/*')
    file_group=[os.path.basename(os.path.dirname(file)) for file in files]
    tag=[os.path.basename(file).replace('.csv','') for file in files]
    tag_lookup=dict(zip(tag, file_group))
    return tag_lookup

def get_point_cloud_df(src):
    '''
    Get point cloud dataframe and add names of proteins as well as tags
    '''
    point_cloud_df=pd.read_csv(src)
    point_cloud_df=point_cloud_df.loc[point_cloud_df.max_radius==4]
    lookup_tag=get_protein_tag_table()
    point_cloud_df['protein_name']=point_cloud_df.protein_tag.apply(lambda x : lookup_tag[x])

    column_order=[
        'protein_tag',
        'protein_name',
        'num_atoms',
        'num_chains',
        'betti_0',
        'betti_1',
        'max_radius',
        'num_vertices',
        'num_simplices'
    ]
    return point_cloud_df.loc[:,column_order]



PARAMS_FP='../config/config.json'

if __name__=='__main__':
    with open(PARAMS_FP, "r") as file:
        params=json.load(file)

    params=params["analysis-params"]
    csv_src=params["csv-src"]
    pdb_src=params["pdb-src"]
    dst=params["dst"]
    test_radii=params["test-radii"]

    arg_parser=argparse.ArgumentParser()
    arg_parser.add_argument(
        '--test', 
        help='Test persistence diagrams from different max radii', 
        action='store_true'
    )
    args=arg_parser.parse_args()

    if args.test:
        warnings.simplefilter('ignore', PDBConstructionWarning)
        persistence(
            csv_src=csv_src,
            pdb_src=pdb_src,
            dst=dst,
            r_max_list=test_radii
        )

import os

from numpy.ma import default_fill_value
import pda
import sys
import gudhi as gd
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib._color_data as mcd
import seaborn as sns
from mpl_toolkits import mplot3d
from scipy.spatial.distance import pdist, squareform
from Bio.PDB import PDBParser
import nglview as nv
import warnings
import pda

def visualize_molecule(data, labels):
    '''
    Visualize molecular 3d coordinates
    ----------
    Parameters
    ----------
        data - list-like object
            List of molecular relative coordinates. 
            Expected dimension of (K, 3), where each column represents a position in 3d space for K atoms.
    '''
    # colors=[
    #     'tab:blue', 
    #     'tab:orange',
    #     'tab:green',
    #     'tab:red',
    #     'tab:purple',
    #     'tab:brown',
    #     'tab:pink',
    #     'tab:gray',
    #     'tab:olive',
    #     'tab:cyan'
    # ]

    colors = list({name for name in mcd.CSS4_COLORS if "xkcd:" + name in mcd.XKCD_COLORS})

    color_map=dict(zip(list(np.unique(labels)), colors[:len(np.unique(labels))]))

    fig=plt.figure()
    ax=plt.axes(projection='3d')
    ax.scatter(
        xs=data[:,0], 
        ys=data[:,1], 
        zs=data[:,2],
        c=[color_map[label] for label in labels]
    )
    plt.show();

def plot_rips_persistence_diagram(
    data=None, simplex_tree=None, radius=None, dimensions=2, dst=None, figure_size=(12,4)):
    '''
    Plot persistence diagram of point cloud
    ----------
    Parameters
    ----------
        data - list-like object
            List of molecular relative coordinates. 
            Expected dimension of (K, 3), where each column represents a position in 3d space for K atoms.
        radius - float (optional)
            Maximum radius to compute Rips Complex from.
            If not passed then 
        dimension - int (optional)
            Number of dimenions to compute persistence diagram of
            If not passed thend defaults to 2.
    '''
    if not data is None:            
        rips_complex=pda.rips_complex(data, radius)
        simplex_tree=rips_complex.create_simplex_tree(max_dimension=dimensions)
        birth_death=simplex_tree.persistence()
    elif not simplex_tree is None:
        birth_death=simplex_tree.persistence()
    else:
        raise ValueError("Null Value Error. Must pass either `data` or `simplex_tree`")



    fig, (ax1, ax2,ax3)=plt.subplots(
        nrows=1,
        ncols=3,
        figsize=figure_size
    )
    gd.plot_persistence_barcode(
        persistence=birth_death, 
        legend=True,
        axes=ax1
    )
    
    gd.plot_persistence_diagram(
        persistence=birth_death, 
        legend=True,
        axes=ax2
    )

    # # axs[1,0]=gd.plot_persistence_density(
    # gd.plot_persistence_density(
    #     persistence=birth_death,
    #     dimension=0, 
    #     legend=True
    #     # axes=axs[1,0]
    # )
    # axs[1,0].set_title("Persistence Density of 0 Dimensional Features")

    # axs[1,1]=gd.plot_persistence_density(
    gd.plot_persistence_density(
        persistence=birth_death,
        dimension=1, 
        legend=True,
        axes=ax3
    )
    fig.savefig('../figures/%s'%dst)

    ax3.set_title("Persistence Density of 1 Dimensional Features")

    plt.subplots_adjust(wspace=0.3)

def load_display_protein(csv_src, pdb_src, r_max, dimensions=2):
    '''
    Display protein diagram, point cloud, and persistence profile
    '''
    warnings.simplefilter('ignore', RuntimeWarning)

    parser=PDBParser()
    structure_name=os.path.basename(csv_src).replace('.csv','')
    structure=parser.get_structure(structure_name, pdb_src)
    df=pd.read_csv(csv_src)
    point_data=df.loc[:,['atom_coord_x','atom_coord_y','atom_coord_z']].values
    labels=df.chain_id.values

    # nv.show_biopython(structure);

    visualize_molecule(point_data, labels)
    
    simplicial_complex=pda.rips_complex(
        data=point_data,
        r_max=r_max,
    )
    simplex_tree=simplicial_complex.create_simplex_tree(max_dimension=dimensions)
    plot_rips_persistence_diagram(
        simplex_tree=simplex_tree,
        radius=r_max,
        dst='%s.png'%structure_name
    )

    print('betti_numbers: ',simplex_tree.betti_numbers())
    print('Number of chains: ', len(list(structure.get_chains())))

    return simplex_tree

def plot_feature_accuracy_by_radius(df):
    '''
    Plot how many features are lost or added as a function of radius
    ----------
    Parameters
    ----------
        df - Pandas DataFrame
            Dataframe containing stats from persistence profiles
    '''

    fig, axes=plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(12,6),
        sharex=True
    )

    sns.lineplot(
        data=df,
        x='max_radius',
        y='extra_features',
        hue='protein_tag',
        ax=axes[0]
    )

    sns.lineplot(
        data=df,
        x='max_radius',
        y='extra_features',
        hue='protein_name',
        ax=axes[1]
    )

    fig.patch.set_facecolor('white')

def plot_feature_count_by_radius(df):
    '''
    Plot how many features are lost or added as a function of radius
    ----------
    Parameters
    ----------
        df - Pandas DataFrame
            Dataframe containing stats from persistence profiles
    '''

    fig, axes=plt.subplots(
        nrows=1,
        ncols=2,
        figsize=(12,6),
        sharex=True
    )

    sns.lineplot(
        data=df,
        x='num_simplices',
        y='extra_features',
        hue='protein_name',
        ax=axes[1]
    )

    fig.patch.set_facecolor('white')
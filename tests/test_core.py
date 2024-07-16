import pytest

from ElliptiCBn.core import _order_nodes
from ElliptiCBn.core import calc_ellipticity
from ElliptiCBn.core import get_macrocycles
from ElliptiCBn.core import get_ellipticity
from ElliptiCBn.core import plot_results


import numpy as np
import pandas as pd

import os

def test__order_nodes():

    cycle_xyz = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]])
    results = _order_nodes(cycle_xyz=cycle_xyz)

    # make sure sorting order is correct
    assert np.array_equal(results,[0,2,1,3,0])

    # make sure array features are right
    assert len(results) == 5
    assert np.array_equal(results[0],results[-1])
    assert len(set(results)) == 4

    # send in different sort order
    cycle_xyz = np.array([[1,0,0],[0,1,0],[-1,0,0]])
    results = _order_nodes(cycle_xyz=cycle_xyz)

    # make sure sorting order is correct
    assert np.array_equal(results,[0,1,2,0])

    # make sure array features are right
    assert len(results) == 4
    assert np.array_equal(results[0],results[-1])
    assert len(set(results)) == 3


def test_calc_ellipticity():
    
    cycle_xyz = np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0]])
    pca_ellipticity, pca_vectors, original_ellipticity = calc_ellipticity(cycle_xyz)
    
    assert np.isclose(pca_ellipticity,0)
    assert np.array_equal(pca_vectors,
                          np.array([[ 0, 0,0],
                                    [ 0, 1,0],
                                    [ 0,-1,0],
                                    [ 1, 0,0],
                                    [-1, 0,0]],dtype=float))
    assert np.isclose(original_ellipticity,0.2118763179059997)

    cycle_xyz = np.array([[3,0,0],[-3,0,0],[0,1,0],[0,-1,0]])
    pca_ellipticity, pca_vectors, original_ellipticity = calc_ellipticity(cycle_xyz)
    
    assert np.isclose(pca_ellipticity,0.888888888888889)
    assert np.array_equal(pca_vectors,
                          np.array([[ 0,    0,0],
                                    [ 3,    0,0],
                                    [-3,    0,0],
                                    [ 0,  1/3,0],
                                    [ 0, -1/3,0]],dtype=float))
    assert np.isclose(original_ellipticity,0.606272373059569)



def test_get_macrocycles(example_xyz):
    

    counter_expect = [np.array([0,1,np.nan]),
                      np.array([0,np.nan]),]

    for counter, xyz in enumerate(example_xyz["*.xyz"]):

        out_df = get_macrocycles(xyz)

        assert issubclass(type(out_df),pd.DataFrame)

        assert np.array_equal(out_df.columns,
                              ["atom_index",
                               "atom_type",
                               "x","y","z",
                               "neighbors",
                               "neigh_pattern",
                               "molecule",
                               "molec_size",
                               "cycle"])
        
        # make sure we're getting the expected number of cycles
        assert np.array_equal(np.unique(out_df["cycle"]),
                                        counter_expect[counter],
                                        equal_nan=True)
        
        
        with pytest.raises(ValueError):
            out_df = get_macrocycles(xyz,
                                     min_num_carbons=20,
                                     max_num_carbons=20)
            
        # Make sure we find no cycles (only na) if we require huge cycles. This 
        # makes sure values for cycles are being passed in
        out_df = get_macrocycles(xyz,
                                 min_num_carbons=2000,
                                 max_num_carbons=20000)
        assert len(np.unique(out_df["cycle"])) == 1

                
def test_get_ellipticity(example_xyz,tmpdir):

    cwd = os.getcwd()
    os.chdir(tmpdir)

    counter_expect = [2,1]

    for counter, xyz in enumerate(example_xyz["*.xyz"]):
        atom_df = get_macrocycles(xyz)
        results, pca_vectors = get_ellipticity(atom_df)
        assert issubclass(type(results),pd.DataFrame)
        assert np.array_equal(results.columns,
                              ["id","size","pca_ellip","orig_ellip","nearby_atoms","bad_protons"])

        assert len(results["id"] == counter_expect[counter])
        assert len(pca_vectors) == counter_expect[counter]

    os.chdir(cwd)
    

def test_plot_results(example_xyz,tmpdir):
    
    cwd = os.getcwd()
    os.chdir(tmpdir)

    for xyz in example_xyz["*.xyz"]:
        
        atom_df = get_macrocycles(xyz)
        _, pca_vectors = get_ellipticity(atom_df)

        # pca vectors plot
        plot_results(atom_df,
                     html_file="yo.html",
                     pca_vector_list=pca_vectors)
        
        assert os.path.isfile("yo.html")
        os.remove("yo.html")

        # no pca vectors plot
        plot_results(atom_df,
                     html_file="yo.html",
                     pca_vector_list=None)
        
        assert os.path.isfile("yo.html")
        os.remove("yo.html")

        # no pca vectors, no molecules or cycles
        plot_results(atom_df,   
                     html_file="yo.html",
                     plot_structures=False,
                     plot_cycles=False)
        
        assert not os.path.isfile("yo.html")

        
    os.chdir(cwd)
    

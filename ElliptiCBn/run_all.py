"""
Wrapper function that runs a complete ElliptiCbn calculation.
"""
from .core import plot_results
from .core import get_macrocycles
from .core import get_ellipticity

import pandas as pd

import os

def _file_check(some_file,
                output_dir=".",
                overwrite=False):

    # If the file does not have a path specified, put in output_dir. If there
    # is a specified path, do *not* put in output_dir. So, /test/this/file.txt
    # and ../file.txt and ./file.txt go to the specified location; file.txt 
    # goes to output_dir/file.txt. 
    if os.path.split(some_file)[0] == "":
        some_file = os.path.join(output_dir,some_file)

    # If the file exists....
    if os.path.exists(some_file):

        # Throw error if we cannot overwrite
        if not overwrite:
            err = f"output file '{some_file}' already exists\n"
            raise FileExistsError(err)
        
        # Or delete the file
        else:
            os.remove(some_file)

    # Return path to file
    return some_file


def run_all(filename,
            min_num_carbons=10,
            max_num_carbons=20,
            guest_search_radius=3,
            summary_file="summary.xlsx",
            output_dir=".",
            overwrite=False):
    """
    Run an ElliptiCbn calculation.
    
    Parameters
    ----------
    filename : str or list
        xyz file name (or list of xyz files) to read
    min_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has less than 
        min_num_carbons
    max_num_carbons : int, default=20
        reject any macrocycle with a central cycle that has more than 
        max_num_carbons
    guest_search_radius : float, default=4
        look for guest atoms within this radius (in angstroms) of the atom 
        centroid.
    summary_file : str, default="summary.xlsx"
        write all cycles to this single summary file if there is more than one 
        xyz file specified. 
    output_dir : str, default="."
        write output to output_dir. 
    overwrite : bool, default=False
        overwrite existing output files
    """

    # Decide whether we have a single file or set of files coming in
    filenames = []
    if hasattr(filename,"__iter__"):
        if issubclass(type(filename),str):
            filenames.append(filename)
        else:
            for f in filename:
                filenames.append(f)

    if len(filenames) == 0:
        err = "No file(s) specified\n"
        raise ValueError(err)

    # If output_dir already exists...
    if os.path.exists(output_dir):

        # If it exists but is not a directory
        if not os.path.isdir(output_dir):

            # If no overwrite, throw error
            if not overwrite:
                err = f"output directory {output_dir} already exists but is not a directory\n"
                raise FileExistsError(err)
            
            # Otherwise, remove and create directory
            else:
                os.remove(output_dir)
                os.mkdir(output_dir)
    
    # Make directory
    else:
        os.mkdir(output_dir)
            
    # If more than one file, make sure the summary_file is defined and appears
    # to be a spreadsheet
    if len(filenames) > 1:

        out_ext = summary_file.split(".")[-1]
        if out_ext not in ["xlsx","csv"]:
            err = f"summary_file {summary_file} must have a .xlsx or .csv extension\n"
            raise ValueError(err)

        # Check to see if file exists
        summary_file = _file_check(some_file=summary_file,
                                   output_dir=output_dir,
                                   overwrite=overwrite)

    # Go through each file
    dfs = []
    for filename in filenames:

        # Check output csv and html files
        file_root = os.path.basename(filename)

        excel_file = f"{file_root}.xlsx"
        excel_file = _file_check(some_file=excel_file,
                                 output_dir=output_dir,
                                 overwrite=overwrite)

        html_file = f"{file_root}.html"
        html_file = _file_check(some_file=html_file,
                                output_dir=output_dir,
                                overwrite=overwrite)
        # Get macrocycles
        atom_df = get_macrocycles(filename,
                                  min_num_carbons=min_num_carbons,
                                  max_num_carbons=max_num_carbons)
        
        # Get ellipticities and write to csv
        ellipticities, pca_vectors = get_ellipticity(atom_df,
                                                     guest_search_radius=guest_search_radius)
        ellipticities.to_excel(excel_file,index=False)

        # Plot results 
        fig = plot_results(atom_df,
                           html_file,
                           pca_vector_list=pca_vectors)
        
        # Record for summary file
        ellipticities["file"] = filename
        dfs.append(ellipticities)

    # If more than one filename,write to output
    if len(dfs) > 1:

        out_df = pd.concat(dfs,ignore_index=True)
        if out_ext == "xlsx":
            out_df.to_excel(summary_file,index=False)
        elif out_ext == "csv":
            out_df.to_csv(summary_file,index=False)
        else:
            err = "`out_ext` should be xlsx or csv\n"
            raise ValueError(err)
        
        return None
    
    return fig
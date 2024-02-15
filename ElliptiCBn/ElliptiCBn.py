"""
Core functions to run ElliptiCBn.
"""

# dependencies
import sklearn
from sklearn.decomposition import PCA 
import plotly.graph_objects as go
import networkx as nx 

# typical scientific computing dependencies
import numpy as np
import pandas as pd

# standard library
import os

def parse_xyz(filename):
    """
    Parse an xyz file and return a dataframe with all atoms. 
    
    Parameters
    ----------
    filename : str
        filename of xyz file to load and parse
    
    Returns
    -------
    atom_df : pandas.DataFrame
        dataframe with all atoms in xyz file. has 'atom_type' carbon with atom
        and x,y,z columns with coordinates. 
    """

    # dictionary to hold atom coordinates as lists of lists
    out_dict = {"atom_type":[],
                "x":[],
                "y":[],
                "z":[]}

    line_counter = 0
    with open(filename) as f:
        for line in f:

            # Skip first two lines
            line_counter += 1
            if line_counter <= 2:
                continue

            # If first entry is C or O, record coordinates as floats
            col = line.split()
            out_dict["atom_type"].append(col[0])
            out_dict["x"].append(float(col[1]))
            out_dict["y"].append(float(col[2]))
            out_dict["z"].append(float(col[3]))

    return pd.DataFrame(out_dict)


def get_network_components(xyz_data,
                           bond_dist=2.5):
    """
    Use networkx to identify individual molecules in the file. 
    Each atom is a node that is connected to other atoms if the 
    distance to a nearby carbon is less than the specified bond distance. 
    By default, this bond distance is set to 3 angstroms. Thus, strongly
    connected components (structures) are atoms that are all connected by less 
    than bond_dist atoms. 
    
    Parameters
    ----------
    xyz_data: numpy.ndarray 
        2-d numpy array containing the xyz coords of all carbon atoms in file
    bond_dist : float, default=2.5
        Maximum angstrom distance to define connected carbons
    
    Returns
    -------
    idx_list: list of ints
        Contains list of numbers for carbons in each component      
    points: list of tuples
        Contains pairs of carbon xyz coordinates that are directly connected
    components: list
        List of lists, with a list for each component identified. xyz points
        are stored as tuples in each component list
    """
    
    dist_matrix = sklearn.metrics.pairwise_distances(xyz_data)

    G = nx.DiGraph()
    for i in range(len(xyz_data)):
        G.add_node(i,coord=xyz_data[i])

    for i in range(len(G.nodes)):
        for j in range(1,len(G.nodes)):
            if dist_matrix[i,j] < bond_dist:
                G.add_edge(i,j)
                G.add_edge(j,i)
                
    components = nx.strongly_connected_components(G)
    component_list = [np.array(list(c),dtype=int) for c in components]

    return component_list


def filter_aspect_ratio(component,
                        atom_xyz,
                        pca_variance_ratio=3):
    """
    Use the aspect ratio of structures to identify plane-like versus rod-like
    structures. Use a principal component analysis to remove structures with
    high variance in only one axis. 
    
    Parameters
    ----------
    component : numpy.ndarray 
        array of indices corresponding to this component on atom_xyz
    atom_xyz : numpy.ndarray
        array of coordinates. first dimension can be indexed by component, 
        second dimension is length 3 and holds x, y, and z. 
    pca_variance_ratio : float, default=3
        keep components where the ratio of pca1/pca2 is < pca_variance_ratio.

    Returns
    -------
    passes : bool
        True if the component has an aspect ratio less than pca_variance_ratio; 
        False otherwise.
    """
    
    # Get component coordinates
    component_xyz = atom_xyz[component]

    # Do 2-dimensional PCA on the component
    pca = PCA(n_components=2)
    pca = pca.fit(np.array(component_xyz))
    
    # Calculate aspect ratio for the two components. 
    variance = pca.explained_variance_ratio_
    ratios = variance[0]/variance[1]
    
    if ratios < pca_variance_ratio:
        return True
    
    return False


def get_central_cycle(component,
                      carbon_xyz,
                      oxygen_xyz,
                      oxygen_dist_cutoff=2.9,
                      min_num_carbons=10,
                      max_num_carbons=20,
                      min_cycle_cc_bond_length=1.3,
                      max_cycle_cc_bond_length=1.7):
    """
    Extract the indices of carbons corresponding to the central cycle in a 
    cucurbituril macrocycle. 

    Parameters
    ----------
    component : numpy.ndarray 
        array of indices corresponding to this component on atom_xyz
    carbon_xyz : numpy.ndarray
        array of carbon coordinates. first dimension can be indexed by
        component, second dimension is length 3 and holds x, y, and z. 
    oxygen_xyz : numpy.ndarray
        array of oxygen coordinates. 
    oxygen_dist_cutoff : float, default=2.9
        when selecting the central cucurbituril macrocycle, identify carbons by
        removing any carbon closer than oxygen_dist_cutoff to an oxygen
    min_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has less than 
        min_num_carbons
    max_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has more than 
        max_num_carbons
    min_cycle_cc_bond_length: float, default=1.3
        minimum length to identify cc bonds in the macrocycle
    max_cycle_cc_bond_length: float, default=1.7
        maximum length to identify cc bonds in the macrocycle

    Returns
    -------
    cycle : numpy.ndarray or None
        returns array of indexes corresponding to central carbon cycle. if no
        cycle is found passing the filter parameters, return None. 
    """

    # Get distances between all component carbons and all oxygens
    component_xyz = carbon_xyz[component]
    co_dist = sklearn.metrics.pairwise_distances(component_xyz,
                                                 oxygen_xyz)
    
    # Create boolean mask holding carbons that are more than oxygen_dist_cutoff
    # away from the nearest oxygen
    to_keep = np.zeros(len(component),dtype=bool)
    for i in range(len(component)):
        min_co_dist = np.min(co_dist[i,:])
        if min_co_dist > oxygen_dist_cutoff:
            to_keep[i] = True

    # Central cycle consists of carbons relatively distant from oxygen
    central_cycle = component[to_keep]
    central_cycle_xyz = carbon_xyz[central_cycle]
    if len(central_cycle_xyz) == 0:
        return None

    # Get distances between all carbons in the central cycle
    cc_dist = sklearn.metrics.pairwise_distances(central_cycle_xyz)

    # Create boolean mask holding carbons that have a single bond to a neighbor
    to_keep = np.zeros(len(central_cycle),dtype=bool)
    for i in range(len(central_cycle)):

        d = cc_dist[i,:]
        min_cc_dist = np.ma.masked_equal(d,0,copy=False).min()
        if min_cycle_cc_bond_length < min_cc_dist and min_cc_dist < max_cycle_cc_bond_length:
            to_keep[i] = True

    # Filter central cycle based on atom-atom distances
    central_cycle = central_cycle[to_keep]

    # Final check of size. Should be an even number of carbons between min_
    # and max_num_carbons
    size = len(central_cycle)
    if size % 2 == 0:
        if size < min_num_carbons or size > max_num_carbons:
            return None
        
    return central_cycle

def _order_nodes(cycle_xyz):
    """
    Create an array with indexes for a cycle in continuous order based on 
    distance. Assumes cycle_xyz forms some kind of cycle. The final array 
    with have length len(cycle_xyz) + 1, with the first and last entries 
    being identical. If you traverse the array from beginning to end, you 
    traverse the cycle step-by-step. Both the starting node and direction 
    (clockwise vs. anticlockwise0 are arbitrary. 
    """

    # Get euclidean distances between all points, then set diagnol to 
    # np.nan
    D = sklearn.metrics.pairwise_distances(cycle_xyz)
    D[np.diag_indices_from(D)] = np.nan

    # Arbitrarily start with node 0.
    nodes = [0]
    
    # For all remaining points...
    for i in range(len(D)-1):
        
        # Get the index of the node that is closest to the last node in 
        # nodes.
        new_node = np.nanargmin(D[nodes[-1]])
        
        # Set all distances from the last node to np.nan (so we will never
        # pick up again). 
        D[nodes[-1],:] = np.nan
        D[:,nodes[-1]] = np.nan
        
        # append the nearest new node to the list
        nodes.append(new_node)
    
    # Append the first node again, forming a complete cycle
    nodes.append(nodes[0])
    
    return np.array(nodes,dtype=int)


def calc_ellipticity(cycle_xyz):
    """
    Calculate the ellipticity of a macrocycle. 
    
    Parameters
    ----------
    cycle_xyz : numpy.ndarray
        L x 3 array holding xyz coordinates of atoms in the macrocycle
    
    Returns
    -------
    pca_ellipticity : float
        ellipticity of the cycle calculated using a PCA analysis
    pca_vectors : numpy.ndarray
        5 x 3 array holding the centroid and then the coordinates of 
        four points away from the centroid defining the PCA vectors in 
        the cycle_xyz coordinate space
    original_ellipticity : float
        ellipticity calculated from the centroid, maximum centroid-
        carbon distance, and cycle perimeter.
    """
    
    # Fit a 2D pca to the macrocycle coordinates
    pca = PCA(n_components=2)
    pca = pca.fit(cycle_xyz)

    # Get variance along both axies
    variance = pca.explained_variance_ratio_

    # Get ellipticity from variance
    pca_ellipticity = (variance[0]-variance[1])/variance[0]

    # Get centroid and all distances to the centroid
    centroid =  np.sum(cycle_xyz,axis=0)/cycle_xyz.shape[0]
    all_dists = np.sqrt(np.sum((cycle_xyz - centroid)**2,axis=1))
        
    # Get fraction to strech each axis
    fx0 = variance[0]/(variance[0]+variance[1])
    fx1 = variance[1]/(variance[0]+variance[1])
    
    # Scalar gives us lenght scale for vectors to plot
    scalar = np.max(all_dists)/fx0
    
    # v0p and v0m are plus and minus vectors from the centroid along the
    # first PCA axis
    v0p = pca.inverse_transform([scalar*fx0,0])
    v0m = pca.inverse_transform([-scalar*fx0,0])
    
    # v1p and v1m are plus and minus vectors from the centroid along the
    # second PCA axis
    v1p = pca.inverse_transform([0,scalar*fx1])
    v1m = pca.inverse_transform([0,-scalar*fx1])
    
    # Build an array holding the centroid and PCA vectors
    pca_vectors = np.row_stack([centroid,v0p,v0m,v1p,v1m])
    
    # calculate the ellipticity using the original method (max centroid-atom 
    # distance and perimeter). 
    
    # Get perimeter
    node_order = _order_nodes(cycle_xyz)
    perim_dists = np.sqrt(np.sum((cycle_xyz[node_order[1:]]-cycle_xyz[node_order[:-1]])**2,axis=1))
    perimeter = np.sum(perim_dists)
    
    # Get a and b and use to calculate ellipticity
    a = np.max(all_dists)
    b = np.sqrt((perimeter**2)/(2*(np.pi**2)) - a**2)
    original_ellipticity = (a-b)/a
    
    return pca_ellipticity, pca_vectors, original_ellipticity
    
    
def get_macrocycles(filename,
                    bond_dist=2.5,
                    aspect_ratio_filter=3,
                    oxygen_dist_cutoff=2.9,
                    min_num_carbons=10,
                    max_num_carbons=20,
                    min_cycle_cc_bond_length=1.3,
                    max_cycle_cc_bond_length=1.7):
    """
    Identify the macrocycles present in an xyz coordinate file. 
    
    Parameters
    ----------
    filename : str
        xyz file name to read
    bond_dist : float, default=2.5
        any atoms closer than bond distance (in angstroms) are identified as 
        part of a single molecule
    aspect_ratio_filter : float, default=3
        reject any identified cycles that have a PCA aspect ratio greater than
        aspect_ratio_filter. An aspect ratio of 1 corresponds to a square; an
        aspect ratio of 10 would be long and skinny. 
    oxygen_dist_cutoff : float, default=2.9
        when selecting the central cucurbituril macrocycle, identify carbons by
        removing any carbon closer than oxygen_dist_cutoff to an oxygen
    min_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has less than 
        min_num_carbons
    max_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has more than 
        max_num_carbons
    min_cycle_cc_bond_length: float, default=1.3
        minimum length to identify cc bonds in the macrocycle
    max_cycle_cc_bond_length: float, default=1.7
        maximum length to identify cc bonds in the macrocycle
    
    Returns
    -------
    atom_df : pandas.DataFrame
        dataframe holding all atoms in the xyz file. the "cycle" column holds
        the macrocyle a particular atom was identifed as being a part of. The 
        "molecule" column indicates idnetified molecules; "molec_size" holds 
        the number of atoms that were found in the containing molecule. 
    """
    
    print(f'Analyzing {filename}.',flush=True)

    # Read all atoms
    atom_df = parse_xyz(filename)
    
    # Prep dataframe for new data
    atom_df["cycle"] = np.nan
    atom_df["molecule"] = np.nan
    atom_df["molec_size"] = np.nan

    # Make an array of heavy atom coordinates
    heavy_index = atom_df.index[atom_df["atom_type"] != "H"]
    heavy_atom_xyz = np.array(atom_df.loc[heavy_index,
                                          ["x","y","z"]])

    # Get strongly connected components for atoms closer than bond_dist
    components = get_network_components(heavy_atom_xyz,
                                        bond_dist=bond_dist)
    
    # Create a list of molecules (list of heavy atom dataframes)
    molecules = []
    for i, c in enumerate(components):
        molecules.append(atom_df.loc[heavy_index[c],:])
        atom_df.loc[heavy_index[c],"molecule"] = i
        atom_df.loc[heavy_index[c],"molec_size"] = len(c)
        

    # Go through each molecule
    cycle_counter = 0
    for m in molecules:

        # Grab carbons and oxygens
        carbon_df = m.loc[m["atom_type"] == "C",:]
        oxygen_df = m.loc[m["atom_type"] == "O",:]

        # Create arrays of carbon and oxygen coordinates 
        carbon_xyz = np.array(carbon_df.loc[:,["x","y","z"]])
        oxygen_xyz = np.array(oxygen_df.loc[:,["x","y","z"]])

        # Too few carbons
        if len(carbon_xyz) < min_num_carbons:
            continue

        # No oxygens. Can't be right molecule type. 
        if len(oxygen_xyz) == 0:
            continue

        # Get strongly connected carbon components within molecule. This may
        # be different than total strongly connected components because we 
        # dropepd some heavy atoms. 
        components = get_network_components(carbon_xyz,
                                            bond_dist=bond_dist)

        # For each component
        for comp in components:

            # component is too tiny to even analyse
            if len(comp) < min_num_carbons:
                continue

            # Filter on aspect ratio (find components that are plane-like rather
            # than rod-like).
            if not filter_aspect_ratio(comp,
                                       atom_xyz=carbon_xyz,
                                       pca_variance_ratio=aspect_ratio_filter):
                continue
        
            
            # Get the central cycle. 
            central_cycle = get_central_cycle(comp,
                                              carbon_xyz=carbon_xyz,
                                              oxygen_xyz=oxygen_xyz,
                                              oxygen_dist_cutoff=oxygen_dist_cutoff,
                                              min_num_carbons=min_num_carbons,
                                              max_num_carbons=max_num_carbons,
                                              min_cycle_cc_bond_length=min_cycle_cc_bond_length,
                                              max_cycle_cc_bond_length=max_cycle_cc_bond_length)

            # No cycle found that matches our search criteria
            if central_cycle is None:
                continue

            # Record that the atoms in question are part of this particular cycle
            idx = atom_df.index[carbon_df.index[central_cycle]]
            atom_df.loc[idx,"cycle"] = cycle_counter
            cycle_counter += 1
    
    print(f"{cycle_counter} macrocycles identified.")
    print("",flush=True)

    return atom_df

def get_ellipticity(atom_df):
    """
    Calculate the ellipticity for all atoms in atom_df.
    
    Parameters
    ----------
    atom_df : pandas.DataFrame
        pandas dataframe holding x,y,z coordiantes of atoms with 
        macrocycles identified in the 'cycle' column 
    
    Returns
    -------
    results : pandas.DataFrame
        dataframe with ellipticities calcualted for all cycles in 
        atom_df. columns are "id" (identifier for cycle), "size" 
        (number of carbon atoms in the central cycle), 
        "pca_ellip" (ellipticity calculated by PCA), and 
        "orig_ellip" (ellipticity calculated using original method
    pca_vectors : list
        list of 5x3 numpy arrays defining ellipse vectors corresponding
        to ellipse a and b. 
    """
    
    # Create lists to store results
    cycle_sizes = []
    pca_ellipticities = []
    pca_vectors = []
    orig_ellipticites = []
    
    # Get sorted array of non-na cycles to look for
    cycle_ids = np.unique(atom_df["cycle"])
    cycle_ids = cycle_ids[np.logical_not(np.isnan(cycle_ids))]
    cycle_ids.sort()
    
    print(f"Calculating ellipticities for {len(cycle_ids)} macrocycles.\n",
          flush=True)    

    # Go through each cycle
    for cycle in cycle_ids:

        # Extract coordinates
        cycle_xyz = np.array(atom_df.loc[atom_df["cycle"] == cycle,
                                        ["x","y","z"]])
        
        # Calcualte ellipticities
        ellip, vec, original_ellip = calc_ellipticity(cycle_xyz)
        
        # Record results
        cycle_sizes.append(len(cycle_xyz))
        pca_ellipticities.append(ellip)
        pca_vectors.append(vec)
        orig_ellipticites.append(original_ellip)
    
    # Create output 
    out_dict = {"id":cycle_ids,
                "size":cycle_sizes,
                "pca_ellip":pca_ellipticities,
                "orig_ellip":orig_ellipticites}
    results =  pd.DataFrame(out_dict)
    
    print("Results:")
    print(results)
    print("",flush=True)

    return results, pca_vectors
        

def plot_results(atom_df,
                 html_file=None,
                 bond_cutoff_dist=1.8,
                 plot_structures=True,
                 plot_cycles=True,
                 min_molecule_size=10,
                 pca_vector_list=None):

    python_file = os.path.abspath(__file__)
    package_dir = os.path.dirname(python_file)
    
    all_gos = []    
    if plot_structures:


        atom_csv = os.path.join(package_dir,"atom-colors.csv")
        color_df = pd.read_csv(atom_csv)
        hex_colors = dict(zip([e.upper() for e in color_df["elem"]],color_df["hex"]))

        atom_colors = []
        for atom in atom_df["atom_type"]:

            try:
                atom_colors.append(hex_colors[atom.upper()])
            except KeyError:
                atom_colors.append("#FF66FF")

        atom_df = atom_df.copy()
        atom_df["color"] = atom_colors
        
        atom_df = atom_df.loc[atom_df["molec_size"] >= min_molecule_size,:]
        
        atom_go = go.Scatter3d(name=None,
                               x=atom_df.x,
                               y=atom_df.y,
                               z=atom_df.z,
                               mode="markers",
                               showlegend=False,
                               marker={"size":5,
                                       "color":atom_df.color})
        all_gos.append(atom_go)


        # Get distances between all atoms
        xyz = np.array(atom_df.loc[:,["x","y","z"]])
        dists = sklearn.metrics.euclidean_distances(xyz)
        
        # Bonds are atoms closer than cutoff. This will be symmetrical about 
        # diagonal matrix. So take only bonds where the index for 0 is less than
        # the index for 1: above the diagonal. 
        bonds = np.argwhere(dists < bond_cutoff_dist)
        bonds = bonds[bonds[:,0] < bonds[:,1]]
        
        # Remove spurious HH bonds
        hh_mask = np.logical_and(np.array(atom_df.loc[atom_df.index[bonds[:,0]],"atom_type"] == "H"),
                                 np.array(atom_df.loc[atom_df.index[bonds[:,1]],"atom_type"] == "H"))
        bonds = bonds[np.logical_not(hh_mask)]

        # Create bond coordinates 
        bond_coord = [[],[],[]]
        for b in bonds:
            
            # For each atom...
            for i in range(2):
                coord = xyz[b[i],:]
                
                # append coordinates in x, y, and z
                for j in range(3):
                    bond_coord[j].append(coord[j])

            # Append "None" so plotly draws segments separately
            for j in range(3):
                bond_coord[j].append(None)

        # Create graphical object for bonds
        bond_go = go.Scatter3d(name=None,
                               x=bond_coord[0],
                               y=bond_coord[1],
                               z=bond_coord[2],
                               mode='lines',
                               showlegend=False,
                               connectgaps=False,
                               line={"color":"gray",
                                     "width":2})

        all_gos.append(bond_go)

    if plot_cycles and "cycle" in atom_df.columns:

        # Read color palette from file
        palette_txt = os.path.join(package_dir,"palette-colors.txt")
        with open(palette_txt) as f:
            color_palette = [c.strip() for c in f.readlines()]
        color_palette = [c for c in color_palette if c != ""]
        
        cycle_counter = 0
        cycles = np.unique(atom_df["cycle"])
        cycles = cycles[np.logical_not(np.isnan(cycles))]
        cycles.sort()
        for cycle in cycles:
    
            cycle_df = atom_df.loc[atom_df.cycle == cycle,:]    
            cycle_xyz = np.array(cycle_df.loc[:,["x","y","z"]])
            
            nodes = cycle_df.index[_order_nodes(cycle_xyz=cycle_xyz)]
            cycle_df = cycle_df.loc[nodes,:]

            cycle_go = go.Scatter3d(name=f"cycle {int(cycle)}",
                                    x=cycle_df.x,
                                    y=cycle_df.y,
                                    z=cycle_df.z,
                                    mode="lines",
                                    line={"width":20,
                                          "color":color_palette[cycle_counter]})
            all_gos.append(cycle_go)
            
            if pca_vector_list is not None:
                
                pca_vectors = pca_vector_list[cycle_counter]

                vector_coord = [[],[],[]]
                for i in range(4):
                    for j in range(3):
                        
                        vector_coord[j].append(pca_vectors[i+1,j])
                        vector_coord[j].append(pca_vectors[0,j])
                        vector_coord[j].append(None)
                
                        # Create graphical object for bonds
                        vector_go = go.Scatter3d(name=None,
                                                 x=vector_coord[0],
                                                 y=vector_coord[1],
                                                 z=vector_coord[2],
                                                 mode='lines',
                                                 showlegend=False,
                                                 connectgaps=False,
                                                 line={"color":color_palette[cycle_counter],
                                                       "width":10})
                
                all_gos.append(vector_go)
                
            cycle_counter += 1
                
                
    fig = None
    if len(all_gos) > 0:

        axes_kwargs = {}
        layout = {'xaxis': axes_kwargs,
                  'yaxis': axes_kwargs,
                  'showlegend':True,
                  'height':800,
                  'width':800}

        fig = go.Figure(data=all_gos,
                        layout=layout)

        fig.update_scenes(xaxis_visible=False,
                          yaxis_visible=False,
                          zaxis_visible=False)

    if html_file is not None:
        print(f"Saving plot to {html_file}",flush=True)
        fig.write_html(html_file)
    
    return fig



# xxxxx
def run_all(filename,
            bond_dist=2.5,
            aspect_ratio_filter=3,
            oxygen_dist_cutoff=2.9,
            min_num_carbons=10,
            max_num_carbons=20,
            min_cycle_cc_bond_length=1.3,
            max_cycle_cc_bond_length=1.7,
            summary_file="summary.csv",
            overwrite=False):
    """
    Identify the macrocycles present in an xyz coordinate file. 
    
    Parameters
    ----------
    filename : str or list
        xyz file name (or list of xyz files) to read
    bond_dist : float, default=2.5
        any atoms closer than bond distance (in angstroms) are identified as 
        part of a single molecule
    aspect_ratio_filter : float, default=3
        reject any identified cycles that have a PCA aspect ratio greater than
        aspect_ratio_filter. An aspect ratio of 1 corresponds to a square; an
        aspect ratio of 10 would be long and skinny. 
    oxygen_dist_cutoff : float, default=2.9
        when selecting the central cucurbituril macrocycle, identify carbons by
        removing any carbon closer than oxygen_dist_cutoff to an oxygen
    min_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has less than 
        min_num_carbons
    max_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has more than 
        max_num_carbons
    min_cycle_cc_bond_length: float, default=1.3
        minimum length to identify cc bonds in the macrocycle
    max_cycle_cc_bond_length: float, default=1.7
        maximum length to identify cc bonds in the macrocycle
    summary_file : str, default="summary.csv"
        write all cycles to this single summary file if there is more than one 
        xyz file specified. 
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

    # If more than one file, make sure the summary file is present
    if len(filenames) > 1:
        out_ext = summary_file.split(".")[-1]
        if out_ext not in ["xlsx","csv"]:
            err = f"summary_file {summary_file} must have a .xlsx or .csv extension\n"
            raise ValueError(err)

        if os.path.exists(summary_file):
            if not overwrite:
                err = f"summary file '{summary_file}' already exists\n"
                raise FileExistsError(err)
            else:
                os.remove(summary_file)

    # Go through each file
    dfs = []
    for filename in filenames:

        # Check output csv and html files
        file_root = os.path.basename(filename)

        csv_file = f"{file_root}.csv"
        if os.path.exists(csv_file):
            if not overwrite:
                err = f"output file '{csv_file}' already exists\n"
                raise FileExistsError(err)
            else:
                os.remove(csv_file)
    
        html_file = f"{file_root}.html"
        if os.path.exists(html_file):
            if not overwrite:
                err = f"output file '{html_file}' already exists\n"
                raise FileExistsError(err)
            else:
                os.remove(html_file)
        
        # Get macrocycles
        atom_df = get_macrocycles(filename,
                                  bond_dist=bond_dist,
                                  aspect_ratio_filter=aspect_ratio_filter,
                                  oxygen_dist_cutoff=oxygen_dist_cutoff,
                                  min_num_carbons=min_num_carbons,
                                  max_num_carbons=max_num_carbons,
                                  min_cycle_cc_bond_length=min_cycle_cc_bond_length,
                                  max_cycle_cc_bond_length=max_cycle_cc_bond_length)
        
        # Get ellipticities and write to csv
        ellipticities, pca_vectors = get_ellipticity(atom_df)
        ellipticities.to_csv(csv_file,index=False)

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

        
        

        

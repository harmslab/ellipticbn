"""
Core functions to run ElliptiCBn.
"""

# dependencies
import sklearn
from sklearn.decomposition import PCA 
import plotly.graph_objects as go
import networkx as nx 
import MDAnalysis as mda

# typical scientific computing dependencies
import numpy as np
import pandas as pd

# standard library
import os

def _order_nodes(cycle_xyz):
    """
    Create an array with indexes for a cycle in continuous order based on 
    distance. Assumes cycle_xyz forms some kind of cycle. The final array 
    with have length len(cycle_xyz) + 1, with the first and last entries 
    being identical. If you traverse the array from beginning to end, you 
    traverse the cycle step-by-step. Both the starting node and direction 
    (clockwise vs. anticlockwise) are arbitrary. 
    """

    # Get euclidean distances between all points, then set diagonal to 
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
    centroid =  np.mean(cycle_xyz,axis=0)
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
                    min_num_carbons=10,
                    max_num_carbons=20):
    """
    Identify the macrocycles present in an xyz coordinate file. 
    
    Parameters
    ----------
    filename : str
        xyz file name to read
    min_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has less than 
        min_num_carbons
    max_num_carbons : int, default=10
        reject any macrocycle with a central cycle that has more than 
        max_num_carbons
    
    Returns
    -------
    atom_df : pandas.DataFrame
        dataframe holding all atoms in the xyz file. the "cycle" column holds
        the macrocyle a particular atom was identifed as being a part of. The 
        "molecule" column indicates idnetified molecules; "molec_size" holds 
        the number of atoms that were found in the containing molecule. 
    """
    
    print(f'Analyzing {filename}.',flush=True)

    if min_num_carbons >= max_num_carbons:
        err = f"min_num_carbons ({min_num_carbons}) must be smaller than max_num_carbons ({max_num_carbons})\n"
        raise ValueError(err)

    allowed_atoms = ["C","N","O","H"]
    chno_lines = []
    non_chno_dict = {"atom_type":[],
                     "x":[],
                     "y":[],
                     "z":[]}
    with open(filename) as f:

        line_counter = -1
        for line in f:
            line_counter += 1
            if line_counter == 1:
                second_line = line

            if line_counter < 2:
                continue

            col = line.split()
            if col[0] in allowed_atoms:
                chno_lines.append(line)
                continue

            non_chno_dict["atom_type"].append(col[0])
            non_chno_dict["x"].append(float(col[1]))
            non_chno_dict["y"].append(float(col[2]))
            non_chno_dict["z"].append(float(col[3]))
    
    chno_lines.insert(0,f"{len(chno_lines)}\n")
    chno_lines.insert(1,second_line)

    with open("tmp-xyz.xyz","w") as f:
        f.write("".join(chno_lines))

    univ = mda.Universe("tmp-xyz.xyz",
                        guess_bonds=True)

    out = {"atom_index":[],
           "atom_type":[],
           "x":[],
           "y":[],
           "z":[],
           "neighbors":[],
           "neigh_pattern":[]}

    G = nx.DiGraph()

    for i, a in enumerate(univ.atoms):

        G.add_node(i)
        
        element = a.element
        x, y, z = a.position
        neighbors = a.bonded_atoms.indices
        neighbor_elem = univ.atoms[neighbors].elements.copy()
        neighbor_elem.sort()

        for n in neighbors:
            G.add_edge(i,n)

        out["atom_index"].append(i)
        out["atom_type"].append(element)
        out["x"].append(x)
        out["y"].append(y)
        out["z"].append(z)
        out["neighbors"].append(neighbors)
        out["neigh_pattern"].append("".join(neighbor_elem))

    # Try to delete temporary file. May fail on windows because mdanalysis 
    # keeps file pointer open and locks file
    try:
        os.remove("tmp-xyz.xyz")
    except PermissionError:
        pass

    atom_df = pd.DataFrame(out)
    molecules = nx.strongly_connected_components(G)
    molecules = [np.array(list(m),dtype=int) for m in molecules]

    atom_df["molecule"] = -1
    atom_df["molec_size"] = -1
    for i, m in enumerate(molecules):
        atom_df.loc[m,"molecule"] = i
        atom_df.loc[m,"molec_size"] = len(m)

    search_patterns = ["CHNN","CNN"]
    cnn_atom_df = atom_df.loc[atom_df["neigh_pattern"].isin(search_patterns),:]
    macrocycles = np.unique(cnn_atom_df["molecule"])

    cycle_counter = 0
    atom_df["cycle"] = np.nan
    for m in macrocycles:

        m_df = cnn_atom_df.loc[cnn_atom_df["molecule"] == m,:]

        size = len(m_df)

        if size % 2 != 0:
            continue
        if size < min_num_carbons or size > max_num_carbons:
            continue

        atom_df.loc[m_df.atom_index,"cycle"] = cycle_counter

        cycle_counter += 1
    
    # Dataframe with coordinates of non chno atoms            
    non_chno_df = pd.DataFrame(non_chno_dict)
    non_chno_df["neighbors"] = [np.array([],dtype=int) for _ in non_chno_df.index]
    non_chno_df["neigh_pattern"] = ""
    non_chno_df["molecule"] = -1
    non_chno_df["molec_size"] = -1
    non_chno_df["cycle"] = np.nan
    non_chno_df["atom_index"] = np.arange(len(non_chno_df.index),dtype=int) + np.max(atom_df["atom_index"])
    non_chno_df = non_chno_df.loc[:,atom_df.columns]

    # Record all atoms, both chno (atom_df) and the ones we ignored earlier
    if len(non_chno_df.index) > 0:
        atom_df = pd.concat([atom_df,non_chno_df],ignore_index=True)
    
    print(f"{cycle_counter} macrocycles identified.")
    print("",flush=True)

    return atom_df

def get_ellipticity(atom_df,guest_search_radius=3):
    """
    Calculate the ellipticity and some quality control for all atoms in atom_df.
    
    Parameters
    ----------
    atom_df : pandas.DataFrame
        pandas dataframe holding x,y,z coordiantes of atoms with 
        macrocycles identified in the 'cycle' column 
    guest_search_radius : float, default=4
        look for guest atoms within this radius (in angstroms) of the atom 
        centroid.
    
    Returns
    -------
    results : pandas.DataFrame
        dataframe with ellipticities calcualted for all cycles in 
        atom_df. columns are "id" (identifier for cycle), "size" 
        (number of carbon atoms in the central cycle), 
        "pca_ellip" (ellipticity calculated by PCA), and 
        "orig_ellip" (ellipticity calculated using original method
        "nearby_atoms" (number of atoms within guest_search_cutoff of the centroid)
        "bad_protons" (number of carbons with protons facing into the cycle rather than out)
    pca_vectors : list
        list of 5x3 numpy arrays defining ellipse vectors corresponding
        to ellipse a and b. 
    """
    
    # Create lists to store results
    cycle_sizes = []
    pca_ellipticities = []
    pca_vectors = []
    orig_ellipticites = []
    nearby_atoms = []
    bad_protons = []
    
    # Get sorted array of non-na cycles to look for
    cycle_ids = np.unique(atom_df["cycle"])
    cycle_ids = cycle_ids[np.logical_not(np.isnan(cycle_ids))]
    cycle_ids.sort()
    
    print(f"Calculating ellipticities for {len(cycle_ids)} macrocycles.\n",
          flush=True)    

    # Go through each cycle
    for cycle in cycle_ids:

        cycle_df = atom_df.loc[atom_df["cycle"] == cycle,:]

        # Extract coordinates
        cycle_xyz = np.array(cycle_df.loc[:,["x","y","z"]])
        
        # Calculate ellipticities
        ellip, vec, original_ellip = calc_ellipticity(cycle_xyz)

        # Get centroid
        centroid =  np.mean(cycle_xyz,axis=0)
        
        # Get atom that are not part of the molecule with the cycle but are near
        # the cycle centroid -- possible guests
        molec = np.unique(atom_df.loc[atom_df["cycle"] == cycle,"molecule"])[0]
        non_cycle_xyz = np.array(atom_df.loc[atom_df["molecule"] != molec,
                                             ["x","y","z"]])
        if len(non_cycle_xyz) == 0:
            nearby = 0
        else:
            D = sklearn.metrics.pairwise_distances([centroid],non_cycle_xyz)
            nearby = np.sum(D < guest_search_radius)

        # Get protons from cycle -- should be facing away from center (negative
        # dot product).
        dot_products = []
        for idx in cycle_df.index:
            
            # Grab a carbon from the cycle
            this_row = cycle_df.loc[idx,:]
            c_vec = np.array(this_row.loc[["x","y","z"]])

            # Get length-one vector going from carbon to centroid
            c_to_cent = centroid - c_vec
            vec_length = np.sqrt(np.sum(c_to_cent**2))
            c_to_cent = c_to_cent/vec_length

            # go through the carbon neighbors
            neighbors = this_row["neighbors"]
            for n in neighbors:

                # If a neighbor id a hydrogen
                neigh_row = atom_df.loc[n,:]
                if neigh_row["atom_type"] == "H":

                    # Get a length-one vector going from carbon to hydrogen
                    h_vec = np.array(neigh_row.loc[["x","y","z"]])
                    c_to_h = h_vec - c_vec
                    vec_length = np.sqrt(np.sum(c_to_h**2))
                    c_to_h = c_to_h/vec_length

                    # Dot product. Will be negative if facing opposite ways,
                    # positive if facing same way
                    dot_products.append(np.dot(c_to_cent,c_to_h))

                    continue
            
        # bad hydrogens have positive dot products
        dot_products = np.array(dot_products)
        bad_h = np.sum(dot_products > 0)

        # Record results
        cycle_sizes.append(len(cycle_xyz))
        pca_ellipticities.append(ellip)
        pca_vectors.append(vec)
        orig_ellipticites.append(original_ellip)
        nearby_atoms.append(nearby)
        bad_protons.append(bad_h)
    
    # Create output 
    out_dict = {"id":cycle_ids,
                "size":cycle_sizes,
                "pca_ellip":pca_ellipticities,
                "orig_ellip":orig_ellipticites,
                "nearby_atoms":nearby_atoms,
                "bad_protons":bad_protons}
    results =  pd.DataFrame(out_dict)
    
    print("Results:")
    print(results)
    print("",flush=True)

    return results, pca_vectors
        

def plot_results(atom_df,
                 html_file=None,
                 plot_structures=True,
                 plot_cycles=True,
                 pca_vector_list=None):
    """
    Plot structures resulting from an ElliptiCBn calculation. 
    
    Parameters
    ----------
    atom_df : pandas.DataFrame
        dataframe holding atoms with x,y,z coordinates. If ellipses are to be 
        drawn, it must have a 'cycle' column as well. 
    html_file : str, optional
        if specified, write the plotly result to an html file
    plot_structures : bool, default=True
        draw molecular structures (atoms/bonds)
    plot_cycles : bool, default=True
        draw calculated ellipses
    pca_vector_list : list, optional
        list of 5x3 PCA array returned by get_ellipticity (one array for each
        ellipse). Used to draw ellipse a and b vectors. If not specified, do 
        not draw vectors. 

    Returns
    -------
    fig : plotly.graph_objects.Figure 
        plotly figure holding plot
    """

    python_file = os.path.abspath(__file__)
    package_dir = os.path.abspath(os.path.join(os.path.dirname(python_file),
                                               "data"))
    
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
        
        atom_go = go.Scatter3d(name=None,
                               x=atom_df.x,
                               y=atom_df.y,
                               z=atom_df.z,
                               mode="markers",
                               showlegend=False,
                               marker={"size":5,
                                       "color":atom_df.color})
        all_gos.append(atom_go)

        bond_list = []
        for idx in atom_df.index:
            row = atom_df.loc[idx,:]
            atom_i = int(row["atom_index"])
            
            for atom_j in row["neighbors"]:
                bond = [atom_i,int(atom_j)]

                bond.sort()
                bond_list.append(tuple(bond))

        bonds = list(set(bond_list))

        # Get all coordinates
        xyz = np.array(atom_df.loc[:,["x","y","z"]])

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

        # Get cycles        
        cycle_counter = 0
        cycles = np.unique(atom_df["cycle"])
        cycles = cycles[np.logical_not(np.isnan(cycles))]
        cycles.sort()

        # Expand color palette unless there are too many cycles
        num_repeats = int(np.ceil(len(cycles)/len(color_palette)))
        if num_repeats < 1:
            num_repeats = 1
        color_palette = color_palette*num_repeats

        # Go through each cycle
        for cycle in cycles:
    
            # Get cycle dataframe and coord
            cycle_df = atom_df.loc[atom_df.cycle == cycle,:]    
            cycle_xyz = np.array(cycle_df.loc[:,["x","y","z"]])
            
            # Get cycle in order of nodes
            nodes = cycle_df.index[_order_nodes(cycle_xyz=cycle_xyz)]
            cycle_df = cycle_df.loc[nodes,:]

            # Create graph object corresponding to cycle 
            cycle_go = go.Scatter3d(name=f"cycle {int(cycle)}",
                                    x=cycle_df.x,
                                    y=cycle_df.y,
                                    z=cycle_df.z,
                                    mode="lines",
                                    line={"width":20,
                                          "color":color_palette[cycle_counter]})
            all_gos.append(cycle_go)
            
            # If a pca vector is passed in, draw them
            if pca_vector_list is not None:
                
                pca_vectors = pca_vector_list[cycle_counter]
                vector_coord = [[],[],[]]

                # For each line to draw...
                for i in range(4):

                    # for each cartesian coordinate...
                    for j in range(3):
                        
                        # Record coordinates
                        vector_coord[j].append(pca_vectors[i+1,j])
                        vector_coord[j].append(pca_vectors[0,j])
                        vector_coord[j].append(None)
                
                        # Create graphical object for pca
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
                
    # create the figure with the graph objects
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

        # Save html if requested
        if html_file is not None:
            print(f"Saving plot to {html_file}",flush=True)
            fig.write_html(html_file)
    
    return fig




        

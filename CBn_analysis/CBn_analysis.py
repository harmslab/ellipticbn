#!/usr/bin/env python3
# coding: utf-8

import numpy as np
import networkx as nx #dependency
from dataclasses import dataclass 
from itertools import combinations
from sklearn.decomposition import PCA #dependency
import plotly.express as px  #dependency
import pandas as pd
import sys
import os




@dataclass
class atom_data:
    name: str
    atoms: list
    x: list
    y: list
    z: list


def parse_XYZ(file):
    """
    Function designed to extract the x, y, and z coordinates and names 
    of all atoms in the file, as well as the name of the file. Coordinates
    in file are assumed to be in angstroms
    
    Parameters
    ----------
    file: xyz file 
          File containing list of str
    
    Returns
    -------
    class object: atom data
                  Dataclass with name as a str, atoms as list, x,y,and z coords as lists
                  
    """
    
    with open(file) as f:
        atoms = []
        xdata = []
        ydata = []
        zdata = []
        
        lines = f.read()
        
        
        split = lines.split('\n')
        
        i = 0
        for line in split[2:-1]:
            line = ','.join(line.split())  
            i +=1
            
            line = line.split(',')
            
            atoms.append(line[0]+str(i))
            xdata.append(float(line[1]))
            ydata.append(float(line[2]))
            zdata.append(float(line[3]))

                


    return atom_data(file,atoms,xdata,ydata,zdata)




def get_xyz_stack(xyz_results):
    """
    This function is used to extract the carbon and oxygen atom data
    to be used in further analysis and visualization of structures in 
    the file. It will use the lists of data contained in the dataclass
    returned from the parse_xyz function to return numpy arrays of 
    the carbon and oxygen coordinates.
    
    Parameters
    ----------
    pdb_results: class object
                 Dataclass with name as a str, atoms as list, x,y,and z coords as lists
    
    Returns
    -------
    carbon_xyz: numpy.ndarray
                2-dimensional numpy array with 3 columns as x, y, z coords of carbon atoms
    
    O_xyz: numpy.ndarray
           2-dimensional numpy array with 3 columns as x, y, z coords of oxygen atoms
           
    """
    
    carbons = []
    O_atoms = []
    
    atom_compiled = {}
    
    # get atom numbers
    for atom in xyz_results.atoms:
        nums = int(''.join([char for char in atom if char.isdigit()]))
        atom_compiled[nums] = []

    # get atom characters
    for key in atom_compiled.keys():
        for atom in xyz_results.atoms:
            num = int(''.join([char for char in atom if char.isdigit()]))
            if num == key:
                atom_compiled[key].append(''.join([char for char in atom if not char.isdigit()]))

                
    for key,atoms in atom_compiled.items():
        for atom in atoms:
            atom_for_index = atom+str(key)

            idx_for_coords = xyz_results.atoms.index(atom_for_index)

            x = xyz_results.x[idx_for_coords]
            y = xyz_results.y[idx_for_coords]
            z = xyz_results.z[idx_for_coords]

            if atom == 'C':
                carbons.append([x,y,z])
            elif atom == 'O':
                O_atoms.append([x,y,z])
            
                
    carbon_xyz = np.array(carbons)
    O_xyz = np.array(O_atoms)


    
    return carbon_xyz, O_xyz




def get_network_components(xyz_data, bond_dist):
    """
    This function uses networkx to identify individual structures in the file.
    Each carbon atom is a node and is connected to another carbon if the 
    distance to a nearby carbon is less than the specified bond distance. 
    By default, this bond distance is set to 3 angstroms. Thus, strongly
    connected components (structures) are clusters of carbon atoms that are
    less than 3 angstroms away from an adjacent carbon.
    
    Parameters
    ----------
    xyz_data: numpy ndarray 
              2-d numpy array containing the xyz coords of all carbon atoms in file
    
    bond_dist: int
               Maximum angstrom distance to define connected carbons
    
    Returns
    -------
    idx_list: list of ints
              Contains list of numbers for carbons in each component
              
    points: list of tuples
            Contains pairs of carbon xyz coordinates that are directly connected
            
    components: list of lists
                Contains a separate list for each component identified. xyz points
                are stored as tuples in each component list
    """
    
    nodes = [(i[0],i[1],i[2]) for i in xyz_data]

    G = nx.DiGraph()

    G.add_nodes_from(nodes)

    atom_atom_pairs = list(combinations(nodes,2))

    cc_bond_length = bond_dist

    distances = []
    points = []

    for pair in atom_atom_pairs:
        x1 = pair[0][0]
        x2 = pair[1][0]

        y1 = pair[0][1]
        y2 = pair[1][1]

        z1 = pair[0][2]
        z2 = pair[1][2]

        xcalc = (x2-x1)**2
        ycalc = (y2-y1)**2
        zcalc = (z2-z1)**2

        d = np.sqrt(xcalc+ycalc+zcalc)

        if d <= cc_bond_length:
            G.add_edge(pair[0],pair[1])
            G.add_edge(pair[1],pair[0])
            points.append([pair[0],pair[1]])
        distances.append(d)
        
    components = nx.strongly_connected_components(G)

    nx.draw(G)
    
    component_list = []
    idx_list = []
    
    for i in components:
        if len(i) > 1:
            component_list.append(list(i))
            indexes_per_component = []

            for xyz_point in i:
                indexes_per_component.append(nodes.index(xyz_point))

            idx_list.append(indexes_per_component)
            
    return idx_list ,points, component_list




def filter_linear(indices, points, components):
    """
    Some files have structures that are highly linear that are not of interest. 
    This function will filter out linear components that we do not need to use 
    in the analysis and visualization. This function uses a principal components
    analysis to remove structures with high variance in only one axis. 
    
    Parameters
    ----------
    indices: list of ints
             Numbers associated with each carbon found in each structure
             
    points: list of tuples
            Contains pairs of carbon xyz coordinates that are directly connected
    
    components: list of lists
                Contains a separate list for each component identified. xyz points
                are stored as tuples in each component list

    Returns
    -------
    trimmed_indices: list of ints
                     Numbers associated with each carbon found in each structure. 
                     Excluding highly linear structures
    
    trimmed_points: list of tuples
                    Contains pairs of carbon xyz coordinates that are directly connected.
                    Excluding highly linear structures
                    
    trimmed_components: list of lists
                        Contains a separate list for each component identified. xyz points
                        are stored as tuples in each component list. Excluding highly
                        linear structures
    """
    
    trimmed_components = []
    trimmed_points = []
    trimmed_indices = []
    
    for idx,comp in enumerate(components):
        new = np.array(comp)
        

        pca = PCA(n_components=2)
        pcomp = pca.fit(new)
        
        
        variance = pca.explained_variance_ratio_
        
        ratios = variance[0]/variance[1]

        if ratios < 3:
            trimmed_components.append(comp)
            trimmed_indices.append(indices[idx])
            
    for point_pair in points:
        for tcomp in trimmed_components:
            if point_pair[0] in tcomp and point_pair[1]  in tcomp:
                trimmed_points.append(point_pair)
        
    return trimmed_indices, trimmed_points, trimmed_components




def calc_centroids(connected_components):
    """
    This function calculates the centroid of each structure using
    xyz coordinates. These centroids are used in a distance calculation
    in another function.
    
    Parameters
    ----------
    connected_components: list of lists
                          Contains a separate list for each component identified. xyz points
                          are stored as tuples in each component list. Excluding highly
                          linear structures
                          
    Returns
    -------
    centroids: list of tuples
               Each tuple in this list contains the x, y, an z coordinates of a structure's
               centroid.
    """
                           
    
    centroids = []
    for comp in connected_components:
        xvals = []
        yvals = []
        zvals = []
        
        for xyz in comp:
            xvals.append(xyz[0])
            yvals.append(xyz[1])
            zvals.append(xyz[2])
        
        length = len(comp)
        
        centroid = (sum(xvals)/length,sum(yvals)/length,sum(zvals)/length)
        
        centroids.append(centroid)
        

    return centroids




def draw_interactive(filename, points, components, centroids,guests=False):
    """
    This function is used to show the structures as interactive plotly graphs 
    in 3 dimensional scatterplots. All carbon atoms and their centroids
    are shown in these graphs. 
    
    Parameters
    ----------
    filename: str
              The name of the file being run
              
    points: list of tuples
            Contains pairs of carbon xyz coordinates that are directly connected.
            Excluding highly linear structures
            
    components: list of lists
                Contains a separate list for each component identified. xyz points
                are stored as tuples in each component list. Excluding highly
                linear structures
                
    centroids: list of tuples
               Each tuple in this list contains the x, y, an z coordinates of a structure's
               centroid.
               
               
    Yields
    -------
    3d scatterplot and saves as HTML file
    
    """
    
    data_centroid = pd.DataFrame(centroids,columns=['x','y','z'])
    centroid_colors = []
    
    for i,centroid in enumerate(centroids):
        centroid_colors.append(f'Centroid {i+1}')

    data_centroid['color'] = centroid_colors    
        
    
    
    already_seen = []
    point_colors = []
    
    for pair in points:
        if pair[0] not in already_seen:
            for idx, c in enumerate(components):
                if pair[0] in c:
                    point_colors.append(f'Structure {idx+1}')
            already_seen.append(pair[0])
            
        if pair[1] not in already_seen:
            for idx, c in enumerate(components):
                if pair[1] in c:
                    point_colors.append(f'Structure {idx+1}')
            already_seen.append(pair[1])
    
    data_points = pd.DataFrame(already_seen,columns=['x','y','z'])
    data_points['color'] = point_colors

    d = pd.concat([data_points,data_centroid])
    
    xlim = [min(d['x'])-1,max(d['x'])+1]
    ylim = [min(d['y'])-1,max(d['y'])+1]
    zlim = [min(d['z'])-1,max(d['z'])+1]
    
    fig = px.scatter_3d(d,x='x',y='y',z='z',color='color')
    fig.update_layout(scene=dict(xaxis=dict(range=xlim),yaxis=dict(range=ylim),zaxis=dict(range=zlim)))
    
    #fig.show()


        
    
    nav = os.path.basename(filename)
    
    if guests:
        fig.write_html(f'./{filename}_data/{nav}_full_rings_with_guests.html')
    else:
        fig.write_html(f'./{filename}_data/{nav}_full_rings_without_guests.html')

    
    return (xlim,ylim,zlim)



def draw_interactive_single(filename, components, oxygens, centroids,limits):
    """
    This function is used to show the structures as interactive plotly graphs 
    in 3 dimensional scatterplots. Only carbon atoms in the center ring
    of each structure is drawn, along with its centroid. These center ring
    components are used for the ellipticity calculation in the calc_write_distances
    function.
    
    Parameters
    ----------
    filename: str
              The name of the file being run
            
    components: list of lists
                Contains a separate list for each component identified. xyz points
                are stored as tuples in each component list. Excluding highly
                linear structures
                
    oxygens: numpy ndarray
             2d numpy array containing xyz coordinates of all oxygen atoms
             
             
    centroids: list of tuples
               Each tuple in this list contains the x, y, an z coordinates of a structure's
               centroid.
               
               
    Returns
    -------
    center_ring_components: list of lists
                            Each list contains the carbons in the center ring of a structure.
                            These are used for ellipticity calculations
    
    
    Yields
    ------
    3d scatterplot and saves as HTML file
    """
    
    
   
    
    data_centroid = pd.DataFrame(centroids,columns=['x','y','z'])
    centroid_colors = []
    
    for i,centroid in enumerate(centroids):
        centroid_colors.append(f'Centroid {i+1}')

    data_centroid['color'] = centroid_colors    
    
    
    
    component_name = []
    high_component_name = [] 
    carbons_to_plot = []
    simple_carbons_to_plot = []
    
    center_ring_components = []
    
    unusual_structure = False
    
    
    for i,comp in enumerate(components):
        distances = []
        component_carbons = []

        for carbon in comp:
            c_dist_to_o = []

            cx = carbon[0]
            cy = carbon[1]
            cz = carbon[2]

            for oxygen in oxygens:
                ox = oxygen[0]
                oy = oxygen[1]
                oz = oxygen[2]

                xcalc = (cx-ox)**2
                ycalc = (cy-oy)**2
                zcalc = (cz-oz)**2

                d = np.sqrt(xcalc+ycalc+zcalc)
                c_dist_to_o.append(d)

            distances.append(min(c_dist_to_o))

            if min(c_dist_to_o) > 2.95:
                component_carbons.append(carbon)
                simple_carbons_to_plot.append(carbon)
                component_name.append(f'Structure {i+1}')


        center_ring_components.append(component_carbons)
 
    
    ### Add conditional here for only doing this step for structures beyond ~20? 24?
    # go through carbon, carbon
    final_component_carbons =[]
    for i, comp in enumerate(center_ring_components):
        c_distances = []
        fixed_comp = []
        
        for carbon in comp:
            cdist_to_c = []
            cx = carbon[0]
            cy = carbon[1]
            cz = carbon[2]
            
            for carbon2 in comp:
                if carbon != carbon2:
                    c2x = carbon2[0]
                    c2y = carbon2[1]
                    c2z = carbon2[2]
                    
                    cxcalc = (cx-c2x)**2
                    cycalc = (cy-c2y)**2
                    czcalc = (cz-c2z)**2
                    
                    d = np.sqrt(cxcalc+cycalc+czcalc)
                    cdist_to_c.append(d)
            c_distances.append(min(cdist_to_c))
            
            if 1.50 < min(cdist_to_c) < 1.57:
                carbons_to_plot.append(carbon)
                fixed_comp.append(carbon)
                high_component_name.append(f'Structure {i+1}')
            
        final_component_carbons.append(fixed_comp)


    
    # if too high carbons, dataframe is carbons to plot and high component name. Otherwise its component carbons and component name
    for comp in center_ring_components:
        if len(comp) > 20:
            unusual_structure = True
    
    if unusual_structure:
        component_dataframe = pd.DataFrame(carbons_to_plot,columns=['x','y','z'])
        component_dataframe['color'] = high_component_name

    else:
        component_dataframe = pd.DataFrame(simple_carbons_to_plot,columns=['x','y','z'])
        component_dataframe['color'] = component_name
        
    d = pd.concat([component_dataframe,data_centroid])
    
  
    fig = px.scatter_3d(d,x='x',y='y',z='z',color='color')
    fig.update_layout(scene=dict(xaxis=dict(range=limits[0]),yaxis=dict(range=limits[1]),zaxis=dict(range=limits[2])))
    #fig.show()
    nav = os.path.basename(filename)
    fig.write_html(f'./{filename}_data/{nav}_single_rings.html')

   
    if unusual_structure:
        return final_component_carbons
    else:
        return center_ring_components





def calc_write_distances(name,structure_data,single_rings,centroid_data):
    """
    This is a function that will calculate the distances of each carbon
    atom from the centroid in each structure, as well as the ellipticity
    of each structure. Once calculated, this function will write this 
    information out to an excel spreadsheet, with each sheet corresponding
    to a given structure and the sheet name as the ellipticity value.
    
    Parameters
    ----------
    name: str
          Name of the file 
          
    structure_data: list of lists
                    Contains a separate list for each component identified. xyz points
                    are stored as tuples in each component list. Excluding highly
                    linear structures
                    
    single_rings: list of lists
                  Each list contains the carbons in the center ring of a structure.
                  These are used for ellipticity calculations
                  
    centroid_data: list of tuples
                   Each tuple in this list contains the x, y, an z coordinates of a structure's
                   centroid.
                   
    Returns
    -------
    compiled_data: pandas DataFrame 
                   DataFrame object with x,y,z coordinates and distances from centroid in every structure
                   
    final_data_list: list of pandas DataFrames
                     List containing DataFrames with x,y,z coordinates and distances from centroid in
                     a given structure
                     
    Yields
    ------
    Excel spreadsheet with each sheet representing a different structure and the ellipticity of that 
    structure as the name of the sheet
    """

    final_data_list = []
    nav = os.path.basename(name)
    writer = pd.ExcelWriter(f'{name[:-4]}_data/{nav[:-4]}_distance_data.xlsx',engine='xlsxwriter')
    
    for idx, structure in enumerate(structure_data):
        distances = []
        xvals = []
        yvals = []
        zvals = []
        for i in range(len(structure)):
            
            d = np.sqrt(sum((np.array(structure[i])-np.array(centroid_data[idx]))**2))
            distances.append(np.round(d,4))
            
            xvals.append(structure[i][0])
            yvals.append(structure[i][1])
            zvals.append(structure[i][2])
            

        data_to_write = pd.DataFrame({'Carbons':['C']*len(xvals),'X':xvals,'Y':yvals,
                                      'Z':zvals,'Dist from centroid':distances})
        # Ellipticity calculation
        new = np.array(single_rings[idx])

        pca = PCA(n_components=2)
        pcomp = pca.fit(new)

        variance = pca.explained_variance_ratio_

        ellipticity = np.round((variance[0]-variance[1])/variance[0],2)
        # end ellipticity calculation
        
        data_to_write.to_excel(writer,sheet_name=f'Ellipticity {idx+1} - {ellipticity}',index=False)
        
        final_data_list.append(data_to_write)
        
    
    writer.save()   #### GOING TO DEPRECATE. NEED CHANGED
    compiled_data = pd.concat(final_data_list)
    
    return compiled_data, final_data_list




def run_all(f):
    """
    Main function that will run everything above
    
    Parameters
    ----------
    f: xyz file
       file on which to run all the above
       
    Returns
    -------
    None
    """

    xyz_data = parse_XYZ(f)

    if not os.path.exists(f'{xyz_data.name[:-4]}_data'):
        os.mkdir(f'{xyz_data.name[:-4]}_data')


    cxyz,oxyz = get_xyz_stack(xyz_data)


    carbon_indices, p, c = get_network_components(cxyz,3)
 
    prelim_cents = calc_centroids(c)


    ti, tp, tc = filter_linear(carbon_indices,p,c)


    cents = calc_centroids(tc)

    if len(tc) != len(c):
        d = draw_interactive(xyz_data.name[:-4],p,c,prelim_cents,guests=True)
        
    lims = draw_interactive(xyz_data.name[:-4],tp,tc,cents)
    f = draw_interactive_single(xyz_data.name[:-4],tc,oxyz,cents,lims)

    full, compiled_data = calc_write_distances(xyz_data.name,tc,f,cents)
    
    return


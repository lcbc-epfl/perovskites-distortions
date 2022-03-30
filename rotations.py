# Funtions to rotate the molecules in the A cation position, useful to avoid all molecules being parallel 
import random
from copy import deepcopy
from pymatgen.core.structure import Structure   
import numpy as np 

def get_center_of_mass(
    structure: Structure,
    sites_indexes: list
    ) -> np.array:
    """
    Calculates center of mass for sites specified as sites_indexes 
    """
    center = np.zeros(3)
    total_weight= 0
    for index in sites_indexes:
        site = structure[index]
        wt = site.species.weight
        center += site.coords * wt
        total_weight += wt
    return center/ total_weight

def get_a_cation_sites(
    struct: Structure,
    a_cation_center_species: str='C',
    radius_cutoff: float = 3.5, # in Angstroms
) -> list:
    """
    Get the indices of the A cation sites in the structure.
    Args:
        struct (Structure): structure
        a_cation_center_species (str, optional): Species of the A cation center. Defaults to 'C'.
        radius_cutoff (float, optional): Radius cutoff for finding the neighbours of the A cation center. \
            Defaults to 3.5 Angstroms.
    """
    # Get the central atom of all A cations, by default Carbon (e.g, A_cation = formamidinium)
    carbon_sites = [(index, site) for index, site in enumerate(struct) if site.species_string == a_cation_center_species]
    fa_molecules = []
    count = 0
    for index, carbon in carbon_sites:
        # Get the indices of all atoms bonded to the central atom
        neighbors = struct.get_neighbors(carbon, r= radius_cutoff)
        print([site.species_string for site in neighbors]) # Make sure we're getting the right atoms
        fa_molecules.append([index,])
        fa_molecules[count] += ([site.index for site in neighbors])
        count += 1
    # assert all elements in list have the same length
    assert len(set([len(molecule) for molecule in fa_molecules])) == 1
    return fa_molecules

def rotate_a_cation(
    struct: Structure,
    fa_molecules: list,
    maximum_angle: int = 180, # in degrees
    rotation_axis: list=[0,1,0],
):
    """ 
    Rotates A cation molecules by random angles chosen between -maximum_angle and +maximum_angle.

    Args:
        struct (Structure): structure
        fa_molecules (list): list of lists with the indices of the sites for each molecule (e.g. [[0,1,2,3], [4,5,6,7]])
        maximum_angle (int, optional): Random angles will be chosen between [-maximum_angle, +maximum_angle] Defaults to 180 degrees.
        rotation_axis (list, optional): Vector to rotate around. \
            For formamidinium, the vector defined between the two Nitrogens (for FA). Defaults to [0,1,0].
    Returns:
        Structure: structure with A cation molecules rotated 
    """
    struct_rotated = deepcopy(struct)
    # Need to perturb structure to break symmetries. \
    # Otherwise, it breaks some of the FA molecules when rotating (does weird shit with the H bonded to the C)
    struct_rotated.perturb(distance = 0.000001)
    for sites in fa_molecules:
        angle = random.randrange(-maximum_angle, maximum_angle)
        #print(f"Rotating by {angle} degrees")
        center_of_mass = get_center_of_mass(struct_rotated, sites)
        struct_rotated.rotate_sites(
            indices = sites, 
            theta = angle * np.pi / 180,
            axis = rotation_axis,
            anchor = center_of_mass, # struct_rotated[sites[0]].coords,
            to_unit_cell = False,
            )
    return struct_rotated
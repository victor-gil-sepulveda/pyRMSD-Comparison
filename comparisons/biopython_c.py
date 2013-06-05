'''
Created on 03/06/2013

@author: victor
'''
import Bio.PDB
import numpy
import time
from pyRMSD.RMSDCalculator import RMSDCalculator

if __name__ == '__main__':
    TRAJECTORY_FILE = "data/ubi_amber_100.pdb"
    
    structures = Bio.PDB.PDBParser().get_structure("UBIQ", TRAJECTORY_FILE)
    number_of_models =  len(structures)
    number_of_atoms = len(list(structures[0].get_atoms()))
    
    # Getting coordinates
    coordinates = []
    for structure in structures:
        structure_coords = []
        for atom in structure.get_atoms():
            coords = [coord for coord in atom.coord]
            structure_coords.append(numpy.array(coords))
        coordinates.append(numpy.array(structure_coords))
    coordinates = numpy.array(coordinates, dtype=numpy.float64)
    
    # Computation with Biopython
    start_t = time.time()
    biopython_rmsd_values = [] 
    super_imposer = Bio.PDB.Superimposer()
    for i in range(0, len(structures)-1):
        reference = structures[i]
        for j in range(i+1, len(structures)):
            mobile = structures[j]
            super_imposer.set_atoms(list(reference.get_atoms()), list(mobile.get_atoms()))
            biopython_rmsd_values.append(super_imposer.rms)
    stop_t = time.time()
    print "Biopython's computation time:", stop_t - start_t

    # Computation with pyRMSD, Biopython uses KABSCH's
    start_t = time.time()
    pyrmsd_rmsd_values = RMSDCalculator(coordinates, "KABSCH_SERIAL_CALCULATOR").pairwiseRMSDMatrix()
    stop_t = time.time()
    print "pyRMSD's computation time:", stop_t - start_t
    
    # Comparison
    numpy.testing.assert_array_almost_equal(biopython_rmsd_values, pyrmsd_rmsd_values, 12)
    print "Done"
    

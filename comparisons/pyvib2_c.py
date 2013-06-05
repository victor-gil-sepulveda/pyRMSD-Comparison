'''
Created on 04/06/2013

@author: victor
'''
from pyviblib.calc.qtrfit import fit
import prody 
import numpy
import time
from pyRMSD.RMSDCalculator import RMSDCalculator

if __name__ == '__main__':
    TRAJECTORY_FILE = "data/ubi_amber_100.pdb"
    
    # We use prody to load the coordinates
    pdb_data = prody.parsePDB(TRAJECTORY_FILE)
    num_atoms = pdb_data.numAtoms()
    total_frames = pdb_data.numCoordsets()
    trajectory = prody.PDBEnsemble("UBI")
    trajectory.setAtoms(pdb_data)
    trajectory.addCoordset(pdb_data.getCoordsets())
    coordinates = numpy.copy(pdb_data.getCoordsets())
    
    # From docs, it needs dim 4 coordinates. It looks like it is a straight FORTRAN conversion.
    pyvib_coords = numpy.concatenate(([[[0]]*num_atoms]*total_frames,coordinates), axis = 2)
    # We want to use all atoms ...
    all_atoms_numbers = numpy.array(zip(range(num_atoms),range(num_atoms)))
    # Preparing weights (from the code, it must be num_atoms +1, again because of the presumed FORTRAN origin...)
    weights = numpy.array([1.]*(num_atoms+1))
    # Matrix creation
    start_t = time.time()
    pyvib_rmsd_data = []
    for i in range(0, total_frames-1):
        for j in range(i+1, total_frames):
            pyvib_rmsd_data.append( fit(all_atoms_numbers,   # Atom pairings
                                          pyvib_coords[i],   # Reference conformation
                                          pyvib_coords[j],   # Fitting conformation
                                          weights            # Weights
                                   )["rms"])   
    stop_t = time.time()
    print "PyVib2's computation time:", stop_t - start_t
    
    # Computation with pyRMSD
    start_t = time.time()
    pyrmsd_rmsd_values = RMSDCalculator(coordinates, "QTRFIT_SERIAL_CALCULATOR").pairwiseRMSDMatrix()
    stop_t = time.time()
    print "pyRMSD's computation time:", stop_t - start_t

    # Comparison
    numpy.testing.assert_array_almost_equal(pyvib_rmsd_data, pyrmsd_rmsd_values, 12)
    print "Done"
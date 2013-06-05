'''
Created on 03/06/2013

@author: victor
'''

import prody
import time
import numpy
from pyRMSD.RMSDCalculator import RMSDCalculator

if __name__ == '__main__':
    TRAJECTORY_FILE = "data/ubi_amber_100.pdb"
    
    # Coordinates loading
    pdb_data = prody.parsePDB(TRAJECTORY_FILE)
    num_atoms = pdb_data.numAtoms()
    total_frames = pdb_data.numCoordsets()
    
    trajectory = prody.PDBEnsemble("UBI")
    trajectory.setAtoms(pdb_data)
    trajectory.addCoordset(pdb_data.getCoordsets())
    coordinates = numpy.copy(pdb_data.getCoordsets())
    
    # RMSD matrix computation
    start_t = time.time()
    prody_rmsd_values = []
    for i in range(0, total_frames):
        trajectory.setCoords(pdb_data.getCoordsets()[i])
        trajectory.superpose()
        prody_rmsd_values.extend(trajectory.getRMSDs()[1:])
        trajectory.delCoordset(0)
    stop_t = time.time()
    print "Prody's computation time:", stop_t - start_t
    
    # Computation with pyRMSD, Prody uses KABSCH's
    start_t = time.time()
    pyrmsd_rmsd_values = RMSDCalculator(coordinates, "KABSCH_SERIAL_CALCULATOR").pairwiseRMSDMatrix()
    stop_t = time.time()
    print "pyRMSD's computation time:", stop_t - start_t

    # Comparison
    numpy.testing.assert_array_almost_equal(prody_rmsd_values, pyrmsd_rmsd_values, 12)
    print "Done"
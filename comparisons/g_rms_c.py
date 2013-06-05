'''
Created on 03/06/2013

@author: victor
'''
import subprocess
import time
import os
from pyRMSD.utils.proteinReading import Reader
from pyRMSD.RMSDCalculator import RMSDCalculator
import numpy
import collections

class XPMConverter:
    """
    Reads an xpm file from g_rms and creates a complete matrix representation.
    
    Original script from: 
        Tsjerk Wassenaar (tsjerkw at gmail.com) 
    
    Modifications: 
        Victor Gil Sepulveda  (victor.gil.sepulveda at gmail.com)
    """
    def XPMConverter(self):
        pass
    
    def unquote(self,s):
        return s[1+s.find('"'):s.rfind('"')]
    
    def uncomment(self,s):
        return s[2+s.find('/*'):s.rfind('*/')]
    
    def col(self,c):
        color = c.split('/*')
        value = self.unquote(color[1])
        color = self.unquote(color[0]).split()
        return color[0], value
    
    def convert(self, xpm_file_handler_in):
        
        # Read in lines until we find the start of the array
        meta = [xpm_file_handler_in.readline()]
        while not meta[-1].startswith("static char *gromacs_xpm[]"):
            meta.append(xpm_file_handler_in.readline())
        
        # The next line will contain the dimensions of the array
        dim = xpm_file_handler_in.readline()
        
        # There are four integers surrounded by quotes
        nx, ny, nc, nb = [int(i) for i in self.unquote(dim).split()] #@UnusedVariable
        
        # The next dim[2] lines contain the color definitions
        # Each pixel is encoded by dim[3] bytes, and a comment
        # at the end of the line contains the corresponding value
        colors = dict([self.col(xpm_file_handler_in.readline()) for i in range(nc)])
        matrix = collections.deque()
        for i in xpm_file_handler_in:
            if i.startswith("/*"):
                continue
            j = self.unquote(i)
            z = [float(colors[j[k:k+nb]]) for k in range(0,nx,nb)]
            # z contains row's values, and this rows go from N to 0 (as we can see in the image,
            # it is reversed)
            matrix.appendleft(z)
        xpm_file_handler_in.close()
        return matrix

if __name__ == '__main__':
    TRAJECTORY_FILE = "data/ubi_amber_1000.pdb"
    
    # Computation with g_rmsd
    # Generate the expect script that will automate input control
    expect_script_str = ("".join(open("data/expect_script","r").readlines()))%(TRAJECTORY_FILE, TRAJECTORY_FILE)
    open("data/expect_script_tmp","w").write(expect_script_str)
    # Use 'expect' to spawn the program
    start_t = time.time()
    subprocess.call(["expect", "data/expect_script_tmp"])
    # Load the matrix
    g_rmsd_matrix = XPMConverter().convert(open("data/matrix.xpm","r"))
    stop_t = time.time()
    print "g_rms's computation time:", stop_t - start_t
    os.system("rm data/matrix.xpm data/rmsd.xvg data/expect_script_tmp")
    
    # Computation with pyRMSD
    start_t = time.time()
    coordinates = Reader().readThisFile(TRAJECTORY_FILE).read()
    RMSDCalculator(coordinates, "QTRFIT_SERIAL_CALCULATOR").oneVsFollowing(0) # This is to mimic g_rmsd pipeline 
    pyrmsd_rmsd_values = RMSDCalculator(coordinates, "QTRFIT_SERIAL_CALCULATOR").pairwiseRMSDMatrix()
    stop_t = time.time()
    print "pyRMSD's computation time:", stop_t - start_t
    
    # Convert g_rmsd matrix to 'condensed'
    dim = len(g_rmsd_matrix)
    c_m_values = []
    for i in range(dim-1):
        for j in range(i+1,dim):
            c_m_values.append(g_rmsd_matrix[i][j])
    
    rmsd = numpy.sqrt(((numpy.array(c_m_values) - numpy.array(pyrmsd_rmsd_values))**2).sum()/len(c_m_values))
    print "RMSD", rmsd
    numpy.testing.assert_almost_equal(c_m_values, pyrmsd_rmsd_values, 3)
    print "Done"
                     
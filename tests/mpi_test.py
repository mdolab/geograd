from mpi4py import MPI
from mpitest import mpitest as mpt 
import numpy as np 


if __name__ == "__main__":
    comm = MPI.COMM_WORLD

    n_long = 30
    A = np.arange(0, n_long*3, dtype=np.float_).reshape(3,n_long)
    print(A[:,0])
    checksum = np.sum(A)

    res = mpt.scat_gath_test(A)
    if (comm.rank == 0):
        print('Correct sum: '+str(checksum))
        print('MPI sum: '+ str(res[0]))

    print(np.sum(A*(1+2+3+4)*2-res[1]))
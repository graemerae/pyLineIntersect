import math
import numpy as np

def lineIntersect3D(PA,PB):
	# Find intersection point of lines in 3D space, in the least squares sense.
	# PA :          Nx3-matrix containing starting point of N lines
	# PB :          Nx3-matrix containing end point of N lines
	# P_Intersect : Best intersection point of the N lines, in least squares sense.
	# distances   : Distances from intersection point to the input lines
	# Graeme Rae 2017
	# Translated to python from Matlab Script originally written by:
	# Anders Eikenes, 2012

	Si = PB - PA #N lines described as vectors
	t1=np.sqrt(np.sum(Si*Si,axis=1))
	ni = Si /  np.transpose(t1*np.ones((3,1)))
	nx = ni[:,0]; ny = ni[:,1]; nz = ni[:,2]
	SXX = np.sum(nx*nx-1)
	SYY = np.sum(ny*ny-1)
	SZZ = np.sum(nz*nz-1)
	SXY = np.sum(nx*ny);
	SXZ = np.sum(nx*nz);
	SYZ = np.sum(ny*nz);
	S = [[SXX,SXY,SXZ],[SXY,SYY,SYZ],[SXZ,SYZ,SZZ]]
	S=np.array(S)
	CX  = np.sum(PA[:,0]*(nx*nx-1) +  PA[:,1]*(nx*ny)  +  PA[:,2]*(nx*nz))
	CY  = np.sum(PA[:,0]*(nx*ny)   +  PA[:,1]*(ny*ny-1) + PA[:,2]*(ny*nz))
	CZ  = np.sum(PA[:,0]*(nx*nz)   +  PA[:,1]*(ny*nz)  +  PA[:,2]*(nz*nz-1))
	C   = [CX,CY,CZ]
	C= np.array(C)
	P_intersect = np.linalg.lstsq(S,C); 
	N=len(PA)
	distances=np.zeros((N,1));
	for i in range(0,N): 
		Pdiff= P_intersect[0]-PA[i,:]
		ui=np.dot(Pdiff,Si[i,:].T)/ np.dot(Si[i,:],Si[i,:].T) ;
		distances[i]=np.linalg.norm(P_intersect[0]-PA[i,:]-ui*Si[i,:]);
	return P_intersect[0],distances.T[0]


#Example 
pa=[[3,3,0],[3,8,0],[3,3,0]]
pb=[[6,5,0],[6,8,0],[6,8,0]]
pa = np.array(pa)
pb = np.array(pb)
Points,Distances=lineIntersect3D(pa,pb)
print Points
print Distances



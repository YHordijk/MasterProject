import numpy as np
import scm.plams as plams



def align_molecule_to_plane(mol, i, j, k, plane='xy'):
	'''
	Aligns a plams.Molecule object to a plane
	i, j, k are indices of atoms forming the plane
	'''
	v = {'xy':np.array([0,0,1]),
		 'yx':np.array([0,0,1]),
		 'xz':np.array([0,1,0]),
		 'zx':np.array([0,1,0]),
		 'yz':np.array([1,0,0]),
		 'zy':np.array([1,0,0])}[plane]

	a, b, c = np.array(mol.atoms[i].coords), np.array(mol.atoms[j].coords), np.array(mol.atoms[k].coords)

	R = align_plane((a-b), (c-b), v)
	mol.rotate(R)
	mol.translate([-x for x in mol.get_center_of_mass()])


def align_plane(a, b, v=np.array([0,0,1])):
	'''
	Construct rotation matrix to align the normal of the plane a-b to v
	'''

	n = np.cross(a, b)
	n = n/np.linalg.norm(n)
	k = np.cross(n, v)
	s = np.linalg.norm(k)
	c = n@v

	I = np.eye(3)
	skew = np.array([
		[0, -k[2], k[1]],
		[k[2], 0, -k[0]],
		[-k[1], k[0], 0]])

	R = I + skew + (skew@skew) * (1-c)/(s*s)

	return R

mol = plams.Molecule(r"C:\Users\Yuman Hordijk\Desktop\Scripts\MasterProject\calculations\achiral_catalyst.H_Ph_AlF3.vacuum\TS.002\output.xyz")
print(mol)
align_molecule_to_plane(mol, 1, 2, 3, 'xy')
print(mol)
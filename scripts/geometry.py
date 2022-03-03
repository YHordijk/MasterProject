import numpy as np
import scm.plams as plams


normal_to_plane = 	{'xy':np.array([0,0,1]),
					 'yx':np.array([0,0,1]),
					 'xz':np.array([0,1,0]),
					 'zx':np.array([0,1,0]),
					 'yz':np.array([1,0,0]),
					 'zy':np.array([1,0,0])}


def align_molecule_to_plane(mol, plane_idx=[0,1,2], plane='xy'):
	'''
	Aligns a plams.Molecule object to a plane
	i, j, k are indices of atoms forming the plane
	'''

	p = normal_to_plane[plane]

	if type(mol) is plams.Molecule:
		a, b, c = [np.array(mol.atoms[i].coords) for i in plane_idx]
		Rplane = align_plane((a-b), (c-b), p)
		mol.rotate(Rplane)
	else:
		a, b, c = [np.array(mol.positions[i]) for i in plane_idx]
		Rplane = align_plane((a-b), (c-b), p)
		mol.apply_rotmat(Rplane)


def rotate_molecule_in_plane(mol, align_idx=[0,1], v=np.array([1,0,0]), plane='xy'):
	'''
	Aligns a molecule to vector v in plane
	By default aligns to the x-axis in the xy-plane
	''' 

	k = normal_to_plane[plane]

	v = v/np.linalg.norm(v)
	if type(mol) is plams.Molecule:
		a, b = [np.array(mol.atoms[i].coords) for i in align_idx]
	else:
		a, b = [np.array(mol.positions[i]) for i in align_idx]
	c = (a-b)
	c = c/np.linalg.norm(c)
	angle = -np.arccos(c@v) 
	if angle < -np.pi/2: angle *= -1
		
	R = np.array([
		[np.cos(angle), -np.sin(angle), 0],
		[np.sin(angle),  np.cos(angle), 0],
		[0, 0, 1]])
	if type(mol) is plams.Molecule:
		mol.rotate(R)
	else:
		mol.apply_rotmat(R)

def center_molecule(mol, center_idx):
	if type(mol) is plams.Molecule:
		d = np.array(mol.atoms[center_idx].coords)
		mol.translate(-d)
	else:
		d = np.array(mol.positions[center_idx])
		mol.center(d)


def align_plane(a, b, v=np.array([0,0,1])):
	'''
	Construct rotation matrix to align the normal of the plane a-b to v
	'''

	n = np.cross(a, b)
	n = n/np.linalg.norm(n)

	return align_axis(n, v)



def align_axis(a, v=np.array([1,0,0])):
	a = a/np.linalg.norm(a)
	k = np.cross(a, v)
	s = np.linalg.norm(k)
	c = a@v

	I = np.eye(3)
	w = skew(k)

	R = I + w + (w@w) * (1-c)/(s*s)

	return R

def skew(v):
	return np.array([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])


# mol = plams.Molecule(r"C:\Users\Yuman Hordijk\Desktop\Scripts\MasterProject\calculations\achiral_catalyst.H_Ph_AlF3.vacuum\TS.002\output.xyz")
# print(mol)
# align_molecule(mol, plane_idx=[3, 2, 1], plane='xy')
# print(mol)
import os
import numpy as np
import pubchempy as pcp
import molviewer2.data as data
import periodictable as pt
import networkx as nx
import matplotlib.pyplot as plt
from math import cos, sin, pi
import molviewer2.screen
import itertools as it


#####	====================	molecule loading	====================	#####

xyz_path = './molviewer2/data/resources/xyz//'

def save_to_xyz(mol, path, comment=''):
	s = f'{mol.natoms}\n'
	s += comment + '\n'
	for e, p in zip(mol.elements, mol.positions):
		s += f'{e:2}\t{p[0]:>10.6f}\t{p[1]:>10.6f}\t{p[2]:>10.6f}\n'

	with open(path, 'w+') as f:
		f.write(s)


def _load_molecule_from_path(path, name=None):
	elements = []
	positions = []
	with open(path, 'r') as xyz:
		for i, l in enumerate(xyz.readlines()[2:]):
			elements.append(l.split()[0])
			positions.append([float(x) for x in l.split()[1:]])

	return Molecule(name=name, elements=elements, positions=np.asarray(positions))


def _load_molecule_from_pubchem(name):
	try:
		name = int(name)
	except:
		pass

	mol = pcp.get_compounds(name, ('name', 'cid')[type(name) is int], record_type='3d')
	if len(mol) == 0:
		print(f'No compound or 3d structure with name {name} found on Pubchem. Please try again with CID.')
		return 

	# print(f'{len(mol)} compounds found on Pubchem with name {name}.')
	mol = mol[0]
	positions = np.asarray([[a.x, a.y, a.z] for a in mol.atoms])
	positions = np.where(positions == None, 0, positions).astype(float)
	elements = np.asarray([a.element for a in mol.atoms])

	mol = Molecule(name=name, elements=elements, positions=positions)
	path = xyz_path + name + '.xyz'
	save_to_xyz(mol, path)

	return mol


def load_molecule(arg):
	#test if its a path
	if os.path.exists(arg):
		return _load_molecule_from_path(arg, name=arg)
	elif os.path.exists(xyz_path + arg + '.xyz'):
		return _load_molecule_from_path(xyz_path + arg + '.xyz', name=arg)
	else:
		return _load_molecule_from_pubchem(arg)


def load_plams_molecule(mol):
	if type(mol) is not list:
		mol = [mol]

	molobjs = []
	for m in mol:
		name = m.name
		elements = [a.symbol for a in m.atoms]
		positions = np.array([a.coords for a in m.atoms])
		molobjs.append(Molecule(name=name, elements=elements, positions=positions))
	return molobjs


#####	====================	helper functions	====================	#####

l2_norm = lambda u, v: np.linalg.norm(u-v)


#####	====================	main molecule class	====================	#####

class Molecule:
	def __init__(self,
				name=None,
				elements=None, 
				positions=None):

		self.name = name
		self.elements = elements
		self.positions = positions

		self.update_molecule()


	def __repr__(self):
		s = f'{self.name}\n'
		for e, p in zip(self.elements, self.positions):
			s += f'{e:2}\t{p[0]:>10.6f}\t{p[1]:>10.6f}\t{p[2]:>10.6f}\n'
		return s


	def update_molecule(self):
		#this should be run anytime changes are made to the molecules positions and elements
		self.atom_numbers = np.asarray([pt.elements.symbol(e).number for e in self.elements])
		self.masses = np.asarray([pt.elements[n].mass for n in self.atom_numbers])
		self.radii = np.asarray([pt.elements[n].covalent_radius for n in self.atom_numbers])
		self.colours = np.asarray([data.ATOM_COLOURS[n] for n in self.atom_numbers])
		self.guess_bond_matrix()		
		self.get_graph_representation()
		self.get_distance_matrix()
		self.get_unique_bonds()
		self.get_unique_bond_angles()
		self.get_pairs_4_apart()
		self.original_pos = self.positions


	####	GRAPH functions	####
	def get_graph_representation(self):
		G = nx.Graph() #graph with atom index as nodes and bonds as edges
		for i, el in enumerate(self.elements):
			G.add_node(i, element=el)
		for i, p in enumerate(self.bond_tuples):
			bo = self.B[p[0],p[1]]
			c = (0,0,0)
			if bo == 1:
				c = (1,0,0)
			elif bo == 2:
				c = (0,1,0)
			elif bo == 3:
				c = (0,0,1)
			G.add_edge(*p, weight=bo, color=c)

		self.graph_rep = G
		return G


	def draw_graph_representation(self, use_colors=True, draw_element_label=True):
		#get atom_colors:
		C = data.ATOM_COLOURS #reference to atomic colours
		#grab atom colors and normalize ([0, 255] -> [0, 1])
		if use_colors: atom_colors = [(float(C[i][0])/255, float(C[i][1])/255, float(C[i][2])/255) for i in self.atom_numbers]
		else: atom_colors = None

		#grab labels 
		if draw_element_label: labels = {i: e for i, e in enumerate(self.elements)}
		else: labels = None

		#get graph_positions
		#initial pos:
		init_pos = {i: (p[0], p[1]) for i, p in enumerate(self.positions)}
		#generate spring layout
		graph_pos = nx.spring_layout(self.graph_rep, iterations=500, dim=2, pos=init_pos)

		#get edge colors:
		edge_color = nx.get_edge_attributes(self.graph_rep, 'color')

		nx.draw(self.graph_rep, pos=init_pos, node_color=atom_colors, labels=labels, edge_color=edge_color.values())
		plt.show()


	def get_distance_matrix(self):
		A = np.zeros((self.natoms, self.natoms))
		G = self.graph_rep
		n = G.nodes()
		for i in range(self.natoms):
			for j in range(i+1, self.natoms):
				try:
					A[i,j] = A[j,i] = len(nx.shortest_path(G, i, j)) - 1
				except:
					A[i,j] = A[j,i] = -1
		self.A = A
		return A


	####	HELPER functions	####
	@property
	def natoms(self):
		return len(self.elements)


	def get_unique_bonds(self):
		B = self.B
		self.unique_bonds = []
		for i in range(self.natoms):
			for o in np.where(B[i,i+1:]>=1)[0]:
				self.unique_bonds.append((i,o+i+1))

		self.unique_bonds = np.asarray(self.unique_bonds)
	

	def get_unique_bond_angles(self):
		B = self.B
		self.unique_bond_angles = []
		for i in range(self.natoms):
			for a, c in list(it.combinations(np.where(B[i,:]>=1)[0], 2)):
				self.unique_bond_angles.append((a, i, c))

		self.unique_bond_angles = np.asarray(self.unique_bond_angles)


	# def get_unique_bond_torsions(self):
	# 	B = self.B
	# 	self.unique_bond_torsions = []
	# 	for i in range(self.natoms):
	# 		for a, c in list(it.combinations(np.where(B[i,i+1:]>=1)[0], 2)):
	# 			self.unique_bond_torsions.append((a+i+1, i, c+i+1))
				
	# 	self.unique_bond_torsions = np.asarray(self.unique_bond_torsions)


	def get_pairs_4_apart(self):
		pairs = []
		for i in range(self.natoms):
			for n in np.where(self.A[i,i+1:] >= 4)[0]:
				pairs.append((i,n+i+1))
			
		self.unique_pairs_4 = np.asarray(pairs)


	def get_rotation_matrix(self, x=None, y=None, z=None):
		R = np.eye(3)

		if not x is None:
			c = cos(x)
			s = sin(x)
			R = R @ np.array(([ 1, 0, 0],
						      [ 0, c,-s],
						      [ 0, s, c]))

		if not y is None:
			c = cos(y)
			s = sin(y)
			R = R @ np.array(([ c, 0, s],
						      [ 0, 1, 0],
						      [-s, 0, c]))

		if not z is None:
			c = cos(z)
			s = sin(z)
			R = R @ np.array(([ c,-s, 0],
						   	  [ s, c, 0],
						   	  [ 0, 0, 1]))

		return R


	@property
	def center_of_mass(self):
		return np.sum(self.masses.reshape(-1,1) * self.positions, 0)/np.sum(self.masses)


	def idx_atoms_by_distance_to_point(self, p):
		dist = np.sum((self.positions-p)**2, axis=1)
		idx = np.argsort(dist)
		return idx


	def unique_atom_types(self):
		return set(self.elements)


	def unique_atom_numbers(self):
		return set(self.atom_numbers)


	####	OPERATIONS on molecule	####
	def shake(self, magnitude=.1):
		self.positions = self.positions + magnitude * np.random.rand(*self.positions.shape)


	def center(self, p=None):
		if p is None:
			p = self.center_of_mass

		self.positions = self.positions - p


	def rotate(self, x=None, y=None, z=None):
		R = self.get_rotation_matrix(x=x,y=y,z=z)
		self.positions = (R @ self.positions.T).T


	def guess_bond_matrix(self, **params):
		#prepare data:
		N = self.natoms 
		B = np.zeros((N,N)).astype(int)
		self.B = B
		max_iter = params.get('max_iter', 5*N)

		#max valences for each element in mol
		max_valences = np.asarray([int(data.MAX_VALENCE[e][0]) for e in self.atom_numbers])

		#covalent radii
		atom_radii = np.asarray([pt.elements[n].covalent_radius for n in self.atom_numbers])

		#calculate cartesian distances
		cart_dist = np.zeros((N,N))
		for i, a1 in enumerate(self.positions):
			for j, a2 in enumerate(self.positions[i:]):
				cart_dist[i, i+j] = cart_dist[i+j, i] = l2_norm(a1, a2)

		#helper functions
		nbonds = lambda i: np.sum(B[i], axis=-1) #number of bonds to atom i
		overbonded = lambda i: nbonds(i) >= max_valences[i]
		saturated = lambda i: nbonds(i) == max_valences[i]
		all_saturated = lambda: all(nbonds(i) == max_valences[i] for i in range(N))
		sort_by_distance_to_atom = lambda a, lst: lst[np.argsort(cart_dist[a][lst])]

		def connected_atoms(i):
			a = np.where(B[i, np.arange(N)] > 0, np.arange(N), None)
			return a[a != None].astype(int)

		def unsaturated_connected_atoms(i):
			a = connected_atoms(i)
			return np.asarray([b for b in a if not saturated(b)])

		#make initial guess for bonds
		#make single bonds based on distance first:
		for i, a1 in enumerate(atom_radii):
			for j, a2 in enumerate(atom_radii[i+1:]):
				if cart_dist[i, i+j+1] < a1 + a2 + 0.4:
					B[i, i+j+1] = B[i+j+1, i] = 1
					if overbonded(i): break

		#generate bond_tuples for use in graph construction
		self.bond_tuples = []
		for i, a1 in enumerate(B):
			for j, a2 in enumerate(a1[i+1:]):
				if a2 > 0: self.bond_tuples.append((i,i+j+1))

		#get hybridisations:
		hybridisations = np.zeros(N).astype(int)
		for a in range(N):
			c = nbonds(a)
			if max_valences[a] == 1:
				hybridisations[a] = 0

			elif max_valences[a] == 2:
				if c == 2:
					hybridisations[a] = 3
				elif c == 1:
					hybridisations[a] = 2

			elif max_valences[a] == 3:
				if c == 3:
					hybridisations[a] = 3
				elif c == 2:
					hybridisations[a] = 2
				elif c == 1:
					hybridisations[a] = 1

			elif max_valences[a] == 4:
				if c == 4:
					hybridisations[a] = 3
				elif c == 3:
					hybridisations[a] = 2
				elif c == 2:
					hybridisations[a] = 1

		sorted_atoms = np.argsort(max_valences) #sort atoms based on max_valences
		
		for iter in range(max_iter):
			if all_saturated():
				print('✔️ Bond order guessing succesful!')
				break #final condition

			B = np.zeros((N,N)).astype(int)
			self.B = B

			#make initial guess for bonds
			#make single bonds based on distance first:
			for i, a1 in enumerate(atom_radii):
				for j, a2 in enumerate(atom_radii[i+1:]):
					if cart_dist[i, i+j+1] < a1 + a2 + 0.4:
						B[i, i+j+1] = B[i+j+1, i] = 1
						if overbonded(i): break

			for a in sorted_atoms:
				if saturated(a): continue #if already satisfied, skip
				#else get unsaturated connected atoms
				unsatconn = unsaturated_connected_atoms(a)
				np.random.shuffle(unsatconn)
				neighbour_hybrids = [hybridisations[x] for x in unsatconn]
				neighbour_hybrids = sorted(neighbour_hybrids, key=lambda x: cart_dist[a,x])

				if hybridisations[a] == 1 and 1 in neighbour_hybrids:
					x = neighbour_hybrids.index(1)
					B[unsatconn[x],a] = 3
					B[a,unsatconn[x]] = 3

				if saturated(a): continue
				#sort them by distance and loop
				for b in unsatconn:
					B[a,b] += 1
					B[b,a] += 1
					if saturated(a): break

		if not all_saturated(): print('❌ Bond order guessing failed 😥')




if __name__ == '__main__':
	mol = load_molecule('benzene')
	print(mol)
	scr = screen.Screen(size=(1600,900))
	# scr.draw_mode = 'simple'
	mol.center()
	# mol.positions = mol.positions + np.random.rand(*mol.positions.shape)*0.5*2-0.25*2
	# mm.optimise(mol)
	scr.draw_molecule(mol)
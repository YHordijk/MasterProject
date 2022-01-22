import numpy as np
import molecule 
import data


def get_uff_atom_types(mol):
	B = mol.B
	types = []
	for i in range(mol.natoms):
		e = mol.elements[i]
		if e == 'H':
			if B[i].sum() > 1:
				types.append('H_b') #briding hydrogen
			else:
				types.append('H_') #generic hydrogen
		elif e == 'He':
			types.append('He4+4') #helium
		elif e == 'Li': 
			types.append('Li')
		elif e == 'Be':
			types.append('Be3+2')
		elif e == 'C':
			if all(B[i] <= 1):
				types.append('C_3')
			if sum(B[i] == 2) >= 1:
				types.append('C_2')
			if sum(B[i] == 3) == 1:
				types.append('C_1')
		elif e == 'O':
			if all(B[i] <= 1):
				types.append('O_3')
			if sum(B[i] == 2) >= 1:
				types.append('O_2')
			if sum(B[i] == 3) == 1:
				types.append('O_1')
		elif e == 'N':
			if all(B[i] <= 1):
				types.append('N_3')
			if sum(B[i] == 2) >= 1:
				types.append('N_2')
			if sum(B[i] == 3) == 1:
				types.append('N_1')
		elif e == 'S':
			if all(B[i] <= 1):
				types.append('S_3+2')
			if sum(B[i] == 2) >= 1:
				types.append('S_2')

	return types



def uff_atom_type_to_idx(types):
	atom_types = FF_UFF_PARAMETERS['Atom']
	idx = []
	for t in types:
		idx.append(atom_types.index(t))
	return np.asarray(idx)



def load_uff_parameters(path):
	keys = ['Atom', 'r1', 'theta0', 'x1', 'D1', 'zeta', 'Z1', 'Vi', 'Uj', 'Xi', 'Hard', 'Radius']
	params = {k:[] for k in keys}
	with open(path, 'r') as prm:
		for l in prm.readlines():
			if l.startswith('param'):
				s = l.strip().split()[1:]
				[params[k].append(s[i]) for i, k in enumerate(keys)]

	return params


class UFFOptimiser:
	def __init__(self, mol):
		self.mol = mol
		self.uff_types = get_uff_atom_types(mol)
		self.uff_types_idx = uff_atom_type_to_idx(uff_types)
		en = np.asarray(list(data.ELECTRONEGATIVITIES.values())).astype(float)
		self.bonds = mol.unique_bonds

	def optimise(self, mol, yield_every=1):
		uff_types = get_uff_atom_types(mol)
		uff_types_idx = uff_atom_type_to_idx(uff_types)
		mol.frames = []
		en = np.asarray(list(data.ELECTRONEGATIVITIES.values())).astype(float)
		
		for i in range(5000):
			#bond length part
			F = np.zeros_like(mol.positions)
			bonds = mol.unique_bonds
			ia = bonds[:,0]
			iauff = uff_types_idx[ia]
			ib = bonds[:,1]
			ibuff = uff_types_idx[ib]
			pos = mol.positions[bonds]
			bond_orders = mol.B[ia,ib]
			a = pos[:,0,:]
			b = pos[:,1,:]
			ab = a - b
			r = np.linalg.norm(ab, axis=1, keepdims=True)
			r1 = np.asarray(np.asarray(FF_UFF_PARAMETERS['r1']), dtype=float).reshape(-1,1)
			ra = r1[iauff]
			rb = r1[ibuff]
			rBO = -0.1332*(ra+rb)*np.log(bond_orders).reshape(-1,1)
			Xi = en[ia]
			Xj = en[ib]
			rEN = ra*rb*(np.sqrt(Xi)-np.sqrt(Xj))**2/(Xi*ra+Xj*rb)
			r0 = ra + rb + rBO + rEN
			Fr = -2*(r-r0)*ab/r
			np.add.at(F, ia, Fr)
			np.add.at(F, ib, -Fr)
			
			#bond angle part
			angles = mol.unique_bond_angles
			pos = mol.positions[angles]
			r1 = pos[:,0,:]
			r2 = pos[:,1,:]
			r3 = pos[:,2,:]
			a = r1 - r2
			b = r3 - r2
			na = 1/np.linalg.norm(a, axis=1, keepdims=True)
			nb = 1/np.linalg.norm(b, axis=1, keepdims=True)
			ahat = a*na
			bhat = b*nb
			ab = np.sum(a*b, axis=1, keepdims=True)
			z = na*nb*ab
			theta = np.arccos(z)
			theta0 = np.asarray(np.asarray(FF_UFF_PARAMETERS['theta0'])[uff_types_idx[angles[:,1]]], dtype=float).reshape(-1,1) * np.pi/180
			Funiv = 2*(theta - theta0)/np.sqrt(1-z*z)
			Fa = Funiv * (na*(bhat-z*ahat))
			Fb = Funiv * (nb*(ahat-z*bhat))

			np.add.at(F, angles[:,0], Fa)
			np.add.at(F, angles[:,1], (-Fa-Fb)/2)
			np.add.at(F, angles[:,2], Fb)


			#vdw part
			pairs = mol.unique_pairs_4
			if pairs.size > 0:
				pos = mol.positions[pairs]
				a = pos[:,0,:]
				b = pos[:,1,:]
				r = np.linalg.norm(a-b, axis=1, keepdims=True)
				chi = np.asarray(np.asarray(FF_UFF_PARAMETERS['x1']), dtype=float).reshape(-1,1)
				D = np.asarray(np.asarray(FF_UFF_PARAMETERS['D1']), dtype=float).reshape(-1,1)
				chii = chi[pairs[:,0]]
				chij = chi[pairs[:,1]]

				chiij = (chii*chij)**.5
				Di = D[pairs[:,0]]
				Dj = D[pairs[:,1]]
				Dij = (Di*Dj)**.5
				r2 = r*r
				zeta = (chiij/r2)**6
				Fvdw = 12*Dij*(zeta-1)*zeta/r2*(a-b)
				np.add.at(F, pairs[:,0], Fvdw)
				np.add.at(F, pairs[:,1], -Fvdw)




			mol.positions = mol.positions + F/100

			mol.center()

			if i%yield_every == 0: yield


FF_UFF_PARAMETERS_PATH = './data/resources/forcefields/UFF/UFF.prm'
FF_UFF_PARAMETERS = load_uff_parameters(FF_UFF_PARAMETERS_PATH)


# if __name__ == '__main__':
# 	get_uff_atom_types(molecule.load_molecule('benzene'))
# 	get_uff_atom_types(molecule.load_molecule('ethene'))
# 	types = get_uff_atom_types(molecule.load_molecule('ethyne'))
import scm.plams as plams
import os, paths 


def write_mol(mol, path):
	with open(path, 'w+') as file:
		file.write(f'{len(mol.atoms)}\n')
		file.write(f'{", ".join(R+"="+n for R, n in mol.substituents.items())}\n')
		for a in mol.atoms:
			file.write(f'{a.symbol}\t{a.coords[0]: .8f}\t{a.coords[1]: .8f}\t{a.coords[2]: .8f}\n')

def sub_path(name):
	return os.path.join(paths.SGT_substituents, name + '.sub')

def reaction_path(name):
	return os.path.join(paths.SGT, name + '.tmplt')


def generate_stationary_points(template_name, substituents={}):
	'''
	This function generates stationary points according to a given template
	It will split the file up into the different molecules
	It will then find substituents at each molecule and try to substitute them

	template_name: name of template file
	substituent: dictionary with substituent group as key and substituent name as value
	'''

	substituents['default'] = 'H'

	def parse_contents(lines, issub=False):
		lines = [l.strip() for l in lines]

		#read from file
		name = lines[0]
		struct = [l.split() for l in lines[1:]]
		elements = [s[0] for s in struct]
		coords = [[float(x) for x in s[1:4]] for s in struct]
		tags = [s[4:] for s in struct]
		
		#create molecule
		mol = plams.Molecule()
		mol.name = name
		atoms = [plams.Atom(symbol=s, coords=c) for s, c in zip(elements, coords)]
		
		if issub: mol.connector = [None, None]
		else: mol.connector = {}

		for i, a, t in zip(range(len(atoms)), atoms, tags):
			a.tag = t
			if not issub:
				for f in t:
					if f.startswith('R'):
						name = f[0:-1]
						if not name in mol.connector: mol.connector[name] = [None, None]
						if f[-1] == 'a': mol.connector[name][0] = a
						elif f[-1] == 'b': mol.connector[name][1] = a

			if issub:
				if t[0] == 'a':   mol.connector[0] = i
				elif t[0] == 'b': mol.connector[1] = i
					
		if issub: mol.connector = tuple(mol.connector)
		else: mol.connector = {R: tuple(c) for R, c in mol.connector.items()}

		[mol.add_atom(a) for a in atoms]
		return mol

	#First parse substituents
	sub_mols = {}
	for R, p in substituents.items():
		with open(sub_path(p), 'r') as file:
			sub_mols[R] = parse_contents(file.readlines(), issub=True)



	with open(reaction_path(template_name), 'r') as reaction:
		content = reaction.read().split('\n\n')
		mols = [parse_contents(c.split('\n')) for c in content]
		for m in mols:
			m.substituents = {R: sub_mols[R].name for R in m.connector.keys()}
			for R in m.connector:
				if not R in sub_mols: 	s = sub_mols['default'].copy()
				else: 					s = sub_mols[R].copy()

				la, lb = s.atoms[s.connector[0]], s.atoms[s.connector[1]]
				m.substitute(m.connector[R], s, (la, lb))
			mpath = os.path.join(paths.xyz, template_name, '_'.join(list(sorted(m.substituents.values()))), f'{m.name}.xyz')
			m.path = mpath
			os.makedirs(os.path.dirname(mpath), exist_ok=True)
			write_mol(m, mpath)

	return mols

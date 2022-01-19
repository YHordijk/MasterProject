import scm.plams as plams
import mol_viewer
import molviewer.molecule as mol


def load_reaction(path):
	with open(path, 'r') as reaction:
		stationary_points = []
		structures = []
		since_last_empty = 0
		for i, l in enumerate(reaction.readlines()):
			if l.strip() == '':
				since_last_empty = 0
				if len(struct) > 0:
					structures.append(struct)
					struct = []
			else:
				since_last_empty += 1

			if since_last_empty == 1:
				stationary_points.append(l.strip())
				struct = []
			if since_last_empty > 1:
				struct.append(l.strip())

		mols = {}
		for sp, struct in zip(stationary_points, structures):
			elements = [s.split()[0] for s in struct]
			coords = [[float(x) for x in s.split()[1:4]] for s in struct]
			mols[sp] = mol.Molecule(name=sp, elements=elements, positions=coords)

		# print(mols[''])
		mol_viewer.MoleculeDrawer(mols['TS']).mainloop()


load_reaction(r"D:\Users\Yuman\Desktop\MasterProject\resources\struct_generator_templates\test.tmplt")
# mol_viewer.MoleculeDrawer('penicillin').mainloop()
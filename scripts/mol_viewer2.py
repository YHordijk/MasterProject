import molviewer2.molecule as molecule
import molviewer2.screen as screen
import struct_generator, sys

def show(mols, simple=False):
	if not type(mols) is list:
		mols = [mols]
	try:
		mols = molecule.load_plams_molecule(mols)
	except: pass
	scr = screen.Screen(size=(1600,900))
	if simple:
		scr.draw_mode = 'simple'

	[mol.center() for mol in mols]
	scr.draw_molecules(mols)


if __name__ == '__main__':
	if len(sys.argv) > 1:
		mol = molecule.load_molecule(sys.argv[1])
		mol.substituents = {}
		for a in sys.argv[2:]:
			f, arg = a.split('=')
			if f == 'n':
				mol.name = arg
			if f.startswith('R'):
				mol.substituents[f] = arg
			if f == 'reaction':
				mol.reaction = arg
				
		show(mol, simple=False)
	else:
		# xyz = molecule.load_molecule(r"D:\Users\Yuman\Desktop\MasterProject\calculations\achiral_catalyst_Cl_Br_ZnCl2\3.achiral_catalyst.Cl_Br_ZnCl2.Rcat\output.xyz")
		xyz = struct_generator.generate_stationary_points('no_catalyst', {'R1':'Cl', 'R2':'Br', 'Rcat':'SnCl4'})
		show(xyz, simple=False)
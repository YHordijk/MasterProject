import molviewer2.molecule as molecule
import molviewer2.screen as screen
import struct_generator


def show(mols, simple=False):
	mols = molecule.load_plams_molecule(mols)
	scr = screen.Screen(size=(1600,900))
	if simple:
		scr.draw_mode = 'simple'
	scr.draw_molecules(mols)


if __name__ == '__main__':
	xyz = struct_generator.generate_stationary_points('achiral_catalyst', {'R1':'Cl', 'R2':'Br', 'Rcat':'SnCl4'})
	show(xyz, simple=False)
import os
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"

import molviewer2.molecule as molecule
import molviewer2.screen as screen
import sys, json
import numpy as np


def show(mols, simple=False):
	if not type(mols) is list:
		mols = [mols]
	try:
		mols = molecule.load_plams_molecule(mols)
	except: raise
	scr = screen.Screen(size=(1600,900))
	if simple:
		scr.draw_mode = 'simple'
	[mol.center() for mol in mols]
	scr.draw_molecules(mols)


def show_results(results, simple=False):
	mol_files = [r.data['opt'].get('output xyz', r.data['opt']['input xyz']) for r in results]
	mols = [molecule.load_molecule(f) for f in mol_files]
	scr = screen.Screen(size=(1600, 900))
	if simple:
		scr.draw_mode = 'simple'
	[mol.center() for mol in mols]
	scr.draw_molecules(mols)


if __name__ == '__main__':
	if len(sys.argv) > 1:
		mol = molecule.load_molecule(sys.argv[1])
		mol.substituents = {}
		stage = 'int'
		for a in sys.argv[2:]:
			f, arg = a.split('=')
			if f == 'n':
				mol.name = arg
			if f.startswith('R'):
				mol.substituents[f] = arg
			if f == 'reaction':
				mol.reaction = arg
			if f == 'normalmode':
				d = np.array([float(x) for x in arg.split(',')])
				mol.normalmode = d.reshape(d.size//3,3)
			if f == 'stage':
				stage = arg
			if f == 'res':
				with open(arg, 'r') as j:
					res = json.load(j)
					mol.name = res['stationary_point']
					mol.substituents = res['substituents']
					mol.reaction = res['reaction']
					if 'vibrations' in res:
						if res['vibrations']['freq_nimag'] > 0 and stage == 'out':
							d = np.array(res['vibrations']['freq_imag_displ'])
							mol.normalmode = d.reshape(d.size//3,3)

				
		show(mol, simple=False)
	else:
		xyz = molecule.load_molecule(r"D:\Users\Yuman\Desktop\MasterProject\resources\struct_generator_templates\urea_tBu_Ph\cat.xyz")
		show(xyz, simple=True)
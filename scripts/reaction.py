import paths
import results_database2 as results_database
from utility import hartree2kcalmol as h2k
import struct_generator, mol_viewer2


def get_reaction_calculations(template, substituents):
	mols = struct_generator.generate_stationary_points(template, substituents)
	hashes = struct_generator.get_hashes_for_calc(template, substituents)
	ids = {}
	results = {}
	with results_database.DatabaseManager() as dbm:
		for st, hash in hashes.items():
			ids[st] = dbm.get_id_by_hash(hash)
			results[st] = dbm.get_by_id(ids[st])
			if results[st] is None: print('not found: ', st)

	return results


def print_energies(results):
	for SP, res in results.items():
		print(f'{SP:5}: ΔG = {h2k(res.energies["gibbs_energy"]):.2f} kcal/mol')


def show_reaction(template, substituents=None, simple=False):
	if substituents is None:
		substituents = {}
	mols = struct_generator.generate_stationary_points(template, substituents, keep_dummy=True)
	[print(m.name, m.flags) for m in mols]
	mol_viewer2.show(mols, simple=simple)

# results = get_reaction_calculations('no_catalyst', {'R1':'H', 'R2':'H'})
# print_energies(results)

# results = get_reaction_calculations('achiral_catalyst', {'R1':'H', 'R2':'H', 'Rcat':'I2'})
# print_energies(results)

show_reaction('urea_tBu_Ph', {'Rch':'S', 'R2':'Ph'}, simple=False)
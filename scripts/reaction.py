import paths, os
import results_database3 as results_database
from utility import hartree2kcalmol as h2k
import struct_generator2, mol_viewer2, job_results3

join = os.path.join

def get_reaction_calculations(template, substituents):
	mols = struct_generator2.generate_stationary_points(template, substituents)
	hashes = struct_generator2.get_hashes_for_calc(template, substituents)
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
	mols = struct_generator2.generate_stationary_points(template, substituents, keep_dummy=True)
	mol_viewer2.show(list(mols.values()), simple=simple)


def get_reaction_results(template, substituents):
	print(f'Searching for reactions matching: reaction={template}, ' + ', '.join(f'{R}={sub}' for R, sub in substituents.items()))
	results = job_results3.get_all_results(calc_dir=join(paths.master, 'calculations_test'))

	#filter the results
	results = [r for r in results if r.reaction == template] #first filter by template
	results = [r for r in results if all(sub == substituents[R] for R, sub in r.substituents)]
	
	print(f'Found {len(results)} results')
	for r in results:
		print(f'\t{r.stationary_point}, {r.status}')
	return results


# results = get_reaction_calculations('no_catalyst', {'R1':'H', 'R2':'H'})
# print_energies(results)

# results = get_reaction_calculations('achiral_catalyst', {'R1':'H', 'R2':'H', 'Rcat':'I2'})
# print_energies(results)

# show_reaction('urea_tBu_Ph', {'Rch':'S', 'R2':'Ph'}, simple=False)
get_reaction_results('urea_tBu_Ph', {'R1':'H', 'R2':'Ph', 'Rch':'O'})
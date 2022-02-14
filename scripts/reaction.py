import paths, os
import results_database3 as results_database
from utility import hartree2kcalmol as h2k
import struct_generator2, mol_viewer2, job_results3
import matplotlib.pyplot as plt


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


def get_reaction_profile(res):
	def plot_path(energies, label=None):
		plt.plot(range(len(energies)), energies, label=label)


	reaction = res[0].reaction
	res = {r.stationary_point: r for r in res}

	if reaction == 'urea_tBu_Ph':
		Renergies = [res['sub'].gibbs + res['cat'].gibbs + res['H'].gibbs + res['rad'].gibbs,
					 res['sub_cat_complex'].gibbs + res['H'].gibbs + res['rad'].gibbs,
					 res['TSR'].gibbs + res['H'].gibbs,
					 res['P1R_cat_complex'].gibbs + res['H'].gibbs,
					 res['P2R_cat_complex'].gibbs,
					 res['P2'].gibbs + res['cat'].gibbs]
		Renergies = [h2k(R - Renergies[0]) for R in Renergies]
		Senergies = [res['sub'].gibbs + res['cat'].gibbs + res['H'].gibbs + res['rad'].gibbs,
					 res['sub_cat_complex'].gibbs + res['H'].gibbs + res['rad'].gibbs,
					 res['TSS'].gibbs + res['H'].gibbs,
					 res['P1S_cat_complex'].gibbs + res['H'].gibbs,
					 res['P2S_cat_complex'].gibbs,
					 res['P2'].gibbs + res['cat'].gibbs]
		Senergies = [h2k(R - Senergies[0]) for R in Senergies]
		labels = ['R', 'RC', 'TS', 'P1C', 'P2C', 'P']


		fig, ax = plt.subplots(1,1) 

		plot_path(Renergies, '(R)')
		plot_path(Senergies, '(S)')

		ax.set_xticks(range(len(labels)))
		ax.set_xticklabels(labels)
		plt.legend()
		plt.show()


def get_all_reactions(res_dir=paths.results):
	...


class Reaction:
	def __init__(self, reaction, substituents):
		self.reaction = reaction 
		self.substituents = substituents
		self.results = get_reaction_results(self.reaction, self.substituents)

	def print_reaction(self):
		print(f'Reaction: {self.reaction}')
		print(f'')


	@property
	def asymmetric(self):
		return len(self.enantiomers) > 0

	@property
	def enantiomers(self):
		return {r.enantiomer for r in self.results if r.enantiomer != 'N/A'}

	@property
	def complete(self):
		return all(r.status in ['Success', 'Warning', 'Failed'] for r in self.results)

	



# show_reaction('squaramide', {'Rch':'S', 'Rch2':'O', 'R1':'H', 'R2':'Ph', 'Rc1':'Ph', 'Rc2':'H'}, simple=False)
# show_reaction('squaramide', {'R1':'H', 'R2':'Ph', 'Rc1':'tBu', 'Rc2':'Ph', 'Rch':'O', 'Rch2':'O'})
# show_reaction('urea_tBu_Ph', {'R1':'H', 'R2':'tBu', 'Rch':'O'})
res = get_reaction_results('urea_tBu_Ph', {'R1':'H', 'R2':'tBu', 'Rch':'O'})
get_reaction_profile(res)

# reaction = Reaction('urea_tBu_Ph', {'R1':'H', 'R2':'tBu', 'Rch':'O'})
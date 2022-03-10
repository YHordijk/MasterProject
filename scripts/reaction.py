import paths, os
import results_database3 as results_database
from utility import hartree2kcalmol as h2k
import struct_generator2, job_results3
import matplotlib.pyplot as plt
import itertools

try: import mol_viewer2
except: raise


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



def get_reaction_results(template, substituents, functional='BLYP-D3(BJ)', basis='TZ2P', numerical_quality='Good', phase='vacuum', silent=False):
	if not silent: print(f'Searching for reactions matching: reaction={template}, ' + ', '.join(f'{R}={sub}' for R, sub in substituents.items()))
	results = job_results3.all_results
	expected_sp = struct_generator2.get_all_stationary_points(template)
	#filter the results
	results = [r for r in results if r.reaction == template] #first filter by template
	results = [r for r in results if all(sub == substituents[R] for R, sub in r.substituents.items())]
	results = [r for r in results if all((r.functional == functional, r.basis == basis, r.numerical_quality == numerical_quality))]
	results = [r for r in results if r.phase == phase]
	results = [r for r in results if r.stationary_point in expected_sp]
	
	if not silent:
		print(f'Found {len(results)} results')
		for r in results:
			print(f'\t{r.stationary_point}, {r.status}')
	return results


def get_reaction_profile(res):
	def plot_path(energies, label=None):
		plt.plot(range(len(energies)), energies, label=label)


	reaction = res[0].reaction
	res = {r.stationary_point: r for r in res}

	fig, ax = plt.subplots(1,1) 

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
		sps = ['R', 'RC', 'TS', 'P1C', 'P2C', 'P']

		trends = {'(R)':Renergies, '(S)':Senergies}

	if reaction == 'no_catalyst':
		energies = [res['sub'].gibbs + res['H'].gibbs + res['rad'].gibbs,
					res['TS'].gibbs + res['H'].gibbs,
					res['P1'].gibbs + res['H'].gibbs,
					res['P2'].gibbs]

		trends = {'':energies}

	[plot(t, l) for t, l in trends.items()]
	plt.gca.set_xticks(range(len(sps)))
	plt.gca.set_xticklabels(sps)
	plt.legend()
	plt.show()


def get_all_reactions(res_dir=paths.results):
	res = job_results3.get_all_results(res_dir=res_dir)
	reaction_names = set(r.reaction for r in res)
	reactions = []
	for rn in reaction_names:
		all_subs = struct_generator2.get_all_substituents(rn)
		bases = set(r.basis for r in res if not r.basis is None)
		functionals = set(r.functional for r in res if not r.functional is None)
		qualities = set(r.numerical_quality for r in res if not r.numerical_quality is None)
		phases = set(r.phase for r in res if not r.phase is None)
		used_subs = {s: set(r.substituents[s] for r in res if s in r.substituents) for s in all_subs}

		#filter based on settings:
		settings_filtered = []
		for basis in bases:
			for functional in functionals:
				for quality in qualities:
					for phase in phases:
						subs_filtered = []
						sub_products = list(itertools.product(*used_subs.values()))

						for product in sub_products:
							substituents = {sn:p for sn, p in zip(used_subs.keys(), product)}
							R = Reaction(rn, substituents, functional, basis, quality)
							if R.complete:
								reactions.append(R)	

	return reactions


class Reaction:
	def __init__(self, reaction, substituents, functional='BLYP-D3(BJ)', basis='TZ2P', numerical_quality='Good', phase='vacuum'):
		self.reaction = reaction 
		self.substituents = substituents
		self.functional = functional
		self.basis = basis 
		self.numerical_quality = numerical_quality
		self.phase = phase
		self.results = get_reaction_results(self.reaction, self.substituents, functional, basis, numerical_quality, silent=True)

	def __repr__(self):
		return f'{self.reaction} ({", ".join(f"{R}={s}" for R, s in self.substituents.items())}) @ ({self.basis}/{self.functional}, numerical_quality={self.numerical_quality})'

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
	def all_stationary_points_present(self):
		return len(self.missing_stationary_points) == 0

	@property
	def all_energies_present(self):
		profile = self.get_reaction_pathways()
		return all(all(e is not None for e in energy.values()) for energy in profile.values())

	@property
	def expected_stationary_points(self):
		return struct_generator2.get_all_stationary_points(self.reaction)
	
	@property
	def missing_stationary_points(self):
		expected_sp = self.expected_stationary_points
		return expected_sp - set(r.stationary_point for r in self.results)

	@property
	def complete(self):
		finished = all(r.status in ['Success', 'Warning', 'Failed', 'Error'] for r in self.results)
		if not finished: return False
		all_present = self.all_stationary_points_present
		if not all_present: return False
		try:
			all_energies = self.all_energies_present
			return all_energies
		except:
			return False

	@property
	def incomplete_reason(self):
		reasons = []
		for sp in self.missing_stationary_points:
			reasons.append(f'Missing stationary point: {sp}')
		for r in self.results:
			if not r.status in ['Success', 'Warning', 'Failed', 'Error']:
				reasons.append(f'Run not finished: {r.stationary_point}')
				continue
			if r.gibbs is None:
				reasons.append(f'Run does not have Gibbs energy: {r.stationary_point}')
		return reasons

	def show(self, simple=False):
		mol_viewer2.show_results(self.results, simple=simple)

	def get_reaction_pathways(self, type='gibbs'):
		res = self.results
		res = {r.stationary_point: r for r in res}
		if type == 'gibbs':
			if self.reaction == 'urea_tBu_Ph':
				Renergies = {'R':  res['sub_cat_complex'].gibbs 	+ res['rad'].gibbs,
							 'TS':  res['TSR'].gibbs,
							 'P': res['P1R_cat_complex'].gibbs}

				Senergies = {'R':  res['sub_cat_complex'].gibbs 	+ res['rad'].gibbs,
							 'TS':  res['TSS'].gibbs,
							 'P': res['P1S_cat_complex'].gibbs}

				trends = {'(R)':Renergies, '(S)':Senergies}


			elif self.reaction == 'no_catalyst':
				energies = {'R': 	res['sub'].gibbs + res['rad'].gibbs,
							'TS': 	res['TS'].gibbs,
							'P': 	res['P1'].gibbs}

				trends = {'':energies}

			elif self.reaction == 'achiral_catalyst':
				energies = {'R':	res['sub_cat_complex'].gibbs + res['rad'].gibbs,
							'TS':	res['TS'].gibbs,
							'P':	res['P1_cat_complex'].gibbs}

				trends = {'':energies}

		elif type == 'bond':
			if self.reaction == 'urea_tBu_Ph':
				Renergies = {'R':  res['sub_cat_complex'].energy 	+ res['rad'].energy,
							 'TS':  res['TSR'].energy,
							 'P': res['P1R_cat_complex'].energy}

				Senergies = {'R':  res['sub_cat_complex'].energy 	+ res['rad'].energy,
							 'TS':  res['TSS'].energy,
							 'P': res['P1S_cat_complex'].energy}

				trends = {'(R)':Renergies, '(S)':Senergies}


			elif self.reaction == 'no_catalyst':
				energies = {'R': 	res['sub'].energy + res['rad'].energy,
							'TS': 	res['TS'].energy,
							'P': 	res['P1'].energy}

				trends = {'':energies}

			elif self.reaction == 'achiral_catalyst':
				energies = {'R':	res['sub_cat_complex'].energy + res['rad'].energy,
							'TS':	res['TS'].energy,
							'P':	res['P1_cat_complex'].energy}

				trends = {'':energies}

		else: NotImplemented
		return trends

	def get_activation_energy(self, type='gibbs'):
		profile = self.get_reaction_pathways(type=type)
		trend_names = profile.keys()

		if self.reaction in ['urea_tBu_Ph', 'achiral_catalyst', 'no_catalyst']:
			act = {n:profile[n]['TS'] - profile[n]['R'] for n in trend_names}
		else: NotImplemented
		return act

	def plot_reaction_profile(self, name=None, format=''):
		def plot_path(energies, label=None):
			plt.plot(range(len(energies)), energies, format, label=label)

		profile = self.get_reaction_pathways()
		trend_names = [f'{name} {k}' for k in profile.keys()]

		if self.reaction in ['urea_tBu_Ph', 'achiral_catalyst', 'no_catalyst']:
			order = ['R', 'TS', 'P']
		else: NotImplemented

		energies = [[t[o] for o in order] for t in profile.values()]
		energies = [[h2k(e-energy[0]) for e in energy] for energy in energies]

		[plot_path(energy, label=l) for l, energy in zip(trend_names, energies)]
		plt.gca().set_xticks(range(len(order)))
		plt.gca().set_xticklabels(order)
		plt.legend()

	def print_energies(self, relative=True, tabs=0):
		profile = self.get_reaction_pathways()
		trend_names = profile.keys()

		if self.reaction in ['urea_tBu_Ph', 'achiral_catalyst', 'no_catalyst']:
			order = ['R', 'TS', 'P']
		else: NotImplemented

		energies = [[t[o] for o in order] for t in profile.values()]
		energies = [[str(round(h2k(e-energy[0]),2)) for e in energy] for energy in energies]
		energy_len = max(max(len(e) for e in energy) for energy in energies)
		trend_len = max(len(n) for n in trend_names)
		header = '\t'*tabs + ' '*(trend_len+1) + ' '.join([o.center(energy_len) for o in order])
		print(header)
		for name, trend in zip(trend_names, energies):
			print(f'{name} {" ".join(e.rjust(energy_len) for e in trend)}')


def get_all_reactions(silent=True):
	reactions = []
	for R1 in ['H', 'F', 'Cl', 'Br', 'I', 'OMe', 'Me', 'NH2']:
		# for R2 in ['H', 'tBu', 'Ph']:
		for R2 in ['H', 'tBu', 'Ph', 'o-FPh', 'm-FPh', 'p-FPh']:
			for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']:
				r = Reaction('achiral_catalyst', {'R1':R1, 'R2':R2, 'Rcat':Rcat})
				if r.complete:
					reactions.append(r)
				else:
					if not silent: print('Reaction not complete', r, r.incomplete_reason)
	return reactions


all_reactions = get_all_reactions()



if __name__ == '__main__':
	plt.figure()
	colors = ['b', 'r', 'g']
	rxns = [Reaction('no_catalyst', {'R1':'H', 'R2':R2}, functional='OLYP', numerical_quality='VeryGood') for R2 in ['H', 'Ph', 'tBu']]
	for i, rxn in enumerate(rxns):
		if not rxn.complete: continue
		rxn.plot_reaction_profile(name=f'R2={rxn.substituents["R2"]} OLYP', format='-'+colors[i])
		# print(rxn, h2k(rxn.get_activation_energy()['']))
	rxns = [Reaction('no_catalyst', {'R1':'H', 'R2':R2}) for R2 in ['H', 'Ph', 'tBu']]
	for i, rxn in enumerate(rxns):
		if not rxn.complete: continue
		rxn.plot_reaction_profile(name=f'R2={rxn.substituents["R2"]} BLYP-D3(BJ)', format='--'+colors[i])
		# print(rxn, h2k(rxn.get_activation_energy()['']))
	plt.title('No catalyst (R1=H)')
	plt.xlabel(r'$\xi$')
	plt.ylabel(r'$\Delta G \quad (kcal/mol)$')
	# plt.show()

	plt.figure()
	rxn_no_cat = Reaction('no_catalyst', {'R1':'H', 'R2':'Ph'})
	rxns = [Reaction('achiral_catalyst', {'R1':'H', 'R2':'Ph', 'Rcat':Rcat}) for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']]
	for rxn in rxns:
		if not rxn.complete: continue
		rxn.plot_reaction_profile(name=f'Cat={rxn.substituents["Rcat"]}')
		# print(rxn, h2k(rxn.get_activation_energy()['']))
	rxn_no_cat.plot_reaction_profile(name=f'No cat', format='--k')
	plt.title('Achiral catalyst (R1=H, R2=Ph) @BLYP-D3(BJ)/TZ2P (Good)')
	plt.xlabel(r'$\xi$')
	plt.ylabel(r'$\Delta G \quad (kcal/mol)$')
	# plt.show()

	plt.figure()
	rxn_no_cat = Reaction('no_catalyst', {'R1':'H', 'R2':'tBu'})
	rxns = [Reaction('achiral_catalyst', {'R1':'H', 'R2':'tBu', 'Rcat':Rcat}) for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']]
	for rxn in rxns:
		if not rxn.complete: continue
		rxn.plot_reaction_profile(name=f'Cat={rxn.substituents["Rcat"]}')
		# print(rxn, h2k(rxn.get_activation_energy()['']))
	rxn_no_cat.plot_reaction_profile(name=f'No cat', format='--k')
	rxn_no_cat.print_energies()
	plt.title('Achiral catalyst (R1=H, R2=tBu) @BLYP-D3(BJ)/TZ2P (Good)')
	plt.xlabel(r'$\xi$')
	plt.ylabel(r'$\Delta G \quad (kcal/mol)$')
	plt.show()


	rxn = Reaction('achiral_catalyst', {'Rcat':'I2', 'R1':'H', 'R2':'Ph'})
	rxn.show()
	# rxn = Reaction('achiral_catalyst', {'Rcat':'TiCl4', 'R1':'H', 'R2':'Ph'})
	# rxn.show()

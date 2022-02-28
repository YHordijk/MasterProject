import reaction, utility
import matplotlib.pyplot as plt


def achiral_catalyst(R1='H', R2='Ph'):
	rxn_no_cat = reaction.Reaction('no_catalyst', {'R1':R1, 'R2':R2})
	rxns = [reaction.Reaction('achiral_catalyst', {'R1':R1, 'R2':R2, 'Rcat':Rcat}) for Rcat in ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']]
	for rxn in rxns:
		if not rxn.complete: continue
		rxn.plot_reaction_profile(name=f'Cat={rxn.substituents["Rcat"]}')
		# print(rxn, h2k(rxn.get_activation_energy()['']))
	rxn_no_cat.plot_reaction_profile(name=f'No cat', format='--k')
	plt.title(f'Achiral catalyst (R1={R1}, R2={R2}) @BLYP-D3(BJ)/TZ2P (Good)')
	plt.xlabel(r'$\xi$')
	plt.ylabel(r'$\Delta G \quad (kcal/mol)$')


def achiral_catalyst_table(R2='Ph'):
	R1s = ['H', 'F', 'Cl', 'Br', 'I']
	Rcats = ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']

	rxns = []
	for R1 in R1s:
		for Rcat in Rcats:
			rxns.append(reaction.Reaction('achiral_catalyst', {'R1':R1, 'R2':R2, 'Rcat':Rcat}))

	for r in rxns:
		if not r.complete:
			print(r, r.incomplete_reason)
	Eacts = [r.get_activation_energy()[''] for r in rxns]
	print(Eacts)

	print(f'E_act table: R2={R2}')
	print('  '.join(R1s))


def achiral_catalyst_Eact(R2='Ph'):
	R1s = ['H', 'F', 'Cl', 'Br', 'I']
	Rcats = ['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']


	rxns = {}
	for R1 in R1s:
		for Rcat in Rcats:
			rxns[(R1, Rcat)] = reaction.Reaction('achiral_catalyst', {'R1':R1, 'R2':R2, 'Rcat':Rcat})

	for s, r in rxns.items():
		if not r.complete:
			print(r, r.incomplete_reason)
	Eacts = {s:r.get_activation_energy()[''] for s, r in rxns.items() if r.complete}

	for R1 in R1s:
		trend = {s[1]: E for s, E in Eacts.items() if s[0] == R1}
		trend = list(sorted(trend.items(), key=lambda x: Rcats.index(x[0])))
		trend = [utility.hartree2kcalmol(t[1]) for t in trend]
		plt.plot(trend, label=R1)

	plt.gca().set_xticks(range(len(Rcats)))
	plt.gca().set_xticklabels(Rcats)
	plt.ylabel(r'$\Delta_{act} \quad (kcal/mol)$')
	plt.xlabel(r'Catalyst')
	plt.title(f'Activation energy trends R2={R2}')
	plt.legend()



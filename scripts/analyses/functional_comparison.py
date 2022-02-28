import reaction
import matplotlib.pyplot as plt



def main():
	colors = ['b', 'r', 'g']
	rxns = [reaction.Reaction('no_catalyst', {'R1':'H', 'R2':R2}, functional='OLYP', numerical_quality='VeryGood') for R2 in ['H', 'Ph', 'tBu']]
	for i, rxn in enumerate(rxns):
		if not rxn.complete: continue
		rxn.plot_reaction_profile(name=f'R2={rxn.substituents["R2"]} OLYP', format='-'+colors[i])
		# print(rxn, h2k(rxn.get_activation_energy()['']))
	rxns = [reaction.Reaction('no_catalyst', {'R1':'H', 'R2':R2}) for R2 in ['H', 'Ph', 'tBu']]
	for i, rxn in enumerate(rxns):
		if not rxn.complete: continue
		rxn.plot_reaction_profile(name=f'R2={rxn.substituents["R2"]} BLYP-D3(BJ)', format='--'+colors[i])
		# print(rxn, h2k(rxn.get_activation_energy()['']))
	plt.title('No catalyst (R1=H)')
	plt.xlabel(r'$\xi$')
	plt.ylabel(r'$\Delta G \quad (kcal/mol)$')



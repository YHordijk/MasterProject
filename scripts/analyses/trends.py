import reaction, utility
import matplotlib.pyplot as plt




def Eact_EDA():
	rxns = reaction.all_reactions
	#we now look at only achiral catalysts
	rxns = [r for r in rxns if r.reaction == 'achiral_catalyst']
	Eact = [r.get_activation_energy()[''] for r in rxns]
	sub_res = []
	for rxn in rxns:
		res = [res for res in rxn.results if res.stationary_point == 'sub_cat_complex'][0]
		if res.data['EDA'] is None:
			print(res.substituents)
		print(res.data['EDA'])
		sub_res.append(res)

	pauli = [r.data['EDA']['pauli'] for r in sub_res]
	bonding = [r.data['EDA']['bonding'] for r in sub_res]
	elstat = [r.data['EDA']['elstat'] for r in sub_res]
	oi = [r.data['EDA']['orbital interaction'] for r in sub_res]

	plt.suptitle('EDA trends')

	plt.subplot(2,2,1)
	plt.scatter(pauli, Eact)
	plt.xlabel('Pauli repulsion (ha)')
	plt.ylabel('Eact (ha)')

	plt.subplot(2,2,2)
	plt.scatter(bonding, Eact)
	plt.xlabel('Interaction energy (ha)')
	plt.ylabel('Eact (ha)')

	plt.subplot(2,2,3)
	plt.scatter(elstat, Eact)
	plt.xlabel('Elstat interaction (ha)')
	plt.ylabel('Eact (ha)')

	plt.subplot(2,2,4)
	plt.scatter(oi, Eact)
	plt.xlabel('Orbital interaction (ha)')
	plt.ylabel('Eact (ha)')

	plt.show()


	plt.suptitle('Atom properties trends')

	mulliken = [r.data['GO']['Mulliken charge'][int(r.data['info']['active atom index'])] for r in sub_res]
	elec_dens = [r.data['GO']['electron dens at nuclei'][int(r.data['info']['active atom index'])] for r in sub_res]
	hirshfeld = [r.data['GO']['Hirshfeld charge'][int(r.data['info']['active atom index'])] for r in sub_res]
	voronoi = [r.data['GO']['Voronoi charge'][int(r.data['info']['active atom index'])] for r in sub_res]
	plt.subplot(2,2,1)
	plt.scatter(mulliken, Eact)
	plt.xlabel('Mulliken charge')
	plt.ylabel('Eact (ha)')

	plt.subplot(2,2,2)
	plt.scatter(elec_dens, Eact)
	plt.xlabel('Elec. density at nucleus')
	plt.ylabel('Eact (ha)')

	plt.subplot(2,2,3)
	plt.scatter(hirshfeld, Eact)
	plt.xlabel('Hirshfeld charge')
	plt.ylabel('Eact (ha)')

	plt.subplot(2,2,4)
	plt.scatter(voronoi, Eact)
	plt.xlabel('Voronoi charge')
	plt.ylabel('Eact (ha)')

	plt.show()
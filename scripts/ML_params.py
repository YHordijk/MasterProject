import paths, os, reaction, csv, job_results3
import matplotlib.pyplot as plt
import matplotlib.cm as cm

join = os.path.join

def help(parameter):
	explanation = {
	'EDA_pauli': 	'Total Pauli repulsion term between substrate and catalyst fragments, unit=[ha]',
	'EDA_bonding': 	'Bonding energy difference between substrate and catalyst fragments, unit=[ha]',
	'EDA_elstat': 	'Electrostatic repulsion term between substrate and catalyst fragments, unit=[ha]',
	'EDA_oi': 		'Orbital interaction term between substrate and catalyst fragments, unit=[ha]',
	'AA_mulliken': 	'Active atom Mulliken charge on substrate/catalyst complex',
	'AA_eldens': 	'Active atom electron density catalyst complex',
	'AA_hirshfeld': 'Active atom Hirshfeld charge on substrate/catalyst complex',
	'AA_voronoi': 	'Active atom Voronoi charge on substrate/catalyst complex',
	'HOMO_contr': 	'Highest contribution of active atom C[2pz] AO on occupied MOs',
	'HOMO_energy': 	'Energy of occupied MO with highest active atom C[2pz] contribution, unit=[ha]',
	'LUMO_contr': 	'Highest contribution of active atom C[2pz] AO on virtual MOs',
	'LUMO_energy': 	'Energy of virtual MO with highest active atom C[2pz] contribution, unit=[ha]',
	}[parameter]

	return explanation


def define_parameters():
	parameters = [
		'EDA_pauli',
		'EDA_bonding',
		'EDA_elstat',
		'EDA_oi',
		'AA_mulliken',
		'AA_eldens',
		'AA_hirshfeld',
		'AA_voronoi',
		'HOMO_contr',
		'HOMO_energy',
		'LUMO_contr',
		'LUMO_energy']

	return parameters


def get_data(results, parameter):
	parameters = define_parameters()
	assert parameter in parameters

	if parameter == 'EDA_pauli':
		return [r.data['EDA']['pauli'] for r in results]
	elif parameter == 'EDA_bonding':
		return [r.data['EDA']['bonding'] for r in results]
	elif parameter == 'EDA_elstat':
		return [r.data['EDA']['elstat'] for r in results]
	elif parameter == 'EDA_oi':
		return [r.data['EDA']['orbital interaction'] for r in results]

	elif parameter == 'AA_mulliken':
		return [r.data['GO']['Mulliken charge'][r.active_atom] for r in results]
	elif parameter == 'AA_eldens':
		return [r.data['GO']['electron dens at nuclei'][r.active_atom] for r in results]
	elif parameter == 'AA_hirshfeld':
		return [r.data['GO']['Hirshfeld charge'][r.active_atom] for r in results]
	elif parameter == 'AA_voronoi':
		return [r.data['GO']['Voronoi charge'][r.active_atom] for r in results]

	elif parameter == 'HOMO_contr':
		return [r.data['SP']['occupied cont'] for r in results]
	elif parameter == 'HOMO_energy':
		return [r.data['SP']['occupied energy'] for r in results]
	elif parameter == 'LUMO_contr':
		return [r.data['SP']['virtual cont'] for r in results]
	elif parameter == 'LUMO_energy':
		return [r.data['SP']['virtual energy'] for r in results]


def generate_params(outfile=paths.training_data):
	rxns = reaction.all_reactions
	#we now look at only achiral catalysts
	rxns = [r for r in rxns if r.reaction == 'achiral_catalyst']
	Eact = [r.get_activation_energy()[''] for r in rxns]
	sub_res = []
	for rxn in rxns:
		res = [res for res in rxn.results if res.stationary_point == 'sub_cat_complex'][0]
		sub_res.append(res)
	respaths = [os.path.relpath(r.path, paths.results) for r in sub_res]

	parameters = define_parameters()

	print(f'Generating {len(parameters)} parameters for {len(sub_res)} reactions')

	data = {}
	for p in parameters:
		print('\t' + p, end='')
		try:
			data[p] = get_data(sub_res, p)
			print('✔️')
		except Exception as e:
			print('❌')
			print('\t\t', type(e), e)

	meta = ['path', 'Eact']
	meta_data = {
		'path': respaths,
		'Eact': Eact
	}

	with open(outfile, 'w+', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(['sep=,'])
		writer.writerow(meta + parameters)

		for i in range(len(sub_res)):
			row = []
			for p in meta:
				row.append(meta_data[p][i])
			for p in parameters:
				row.append(data[p][i])
			writer.writerow(row)


def read_data(infile=paths.training_data):
	data = []
	with open(infile, newline='') as file:
		next(file)
		reader = csv.DictReader(file)
		for row in reader:
			data.append(row)

	return data


def plot_trend(data, param, colors=None, labels=None, unit=''):
	used_labels = set()
	for i, d in enumerate(data):
		x = float(d[param])
		y = float(d['Eact'])
		if colors is not None: c = colors[i]
		else: c = None
		if labels is not None: 
			l = labels[i]
			if l in used_labels: 
				l = None
			else:
				used_labels.add(l)
		else: 
			l = None
		
		plt.scatter(x, y, color=c, label=l)

	plt.xlabel(param)
	plt.ylabel('Eact')
	

def get_results(data):
	respaths = [d['path'] for d in data]
	res = [r for r in job_results3.all_results if os.path.relpath(r.path, paths.results) in respaths]
	res = sorted(res, key=lambda r: respaths.index(os.path.relpath(r.path, paths.results)))
	return res



if __name__ == '__main__':
	generate_params()
	print('Explanation of parameters:')
	pl = max(len(p) for p in define_parameters())
	for p in define_parameters():
		print(f'\t{p.ljust(pl+4)} {help(p)}')

	cmap = cm.get_cmap('Set1')
	data = read_data()
	res = get_results(data)
	all_R1s = list(set(r.substituents['R1'] for r in res))
	colors = cmap([all_R1s.index(r.substituents['R1'])/len(all_R1s) for r in res])
	labels = [r.substituents['R1'] for r in res]

	plt.subplot(2,2,1)
	plot_trend(data, 'EDA_pauli', colors=colors, labels=labels)
	plt.subplot(2,2,2)
	plot_trend(data, 'EDA_bonding', colors=colors, labels=labels)
	plt.subplot(2,2,3)
	plot_trend(data, 'EDA_elstat', colors=colors, labels=labels)
	plt.subplot(2,2,4)
	plot_trend(data, 'EDA_oi', colors=colors, labels=labels)

	plt.legend()
	plt.show()
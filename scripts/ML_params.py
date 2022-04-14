import paths, os, reaction, csv, job_results3
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from utility import hartree2eV as h2e
from utility import hartree2kcalmol as h2k
from scipy import stats
import numpy as np
from regression import MLutils

join = os.path.join

def help(parameter):
	explanation = {
	'EDA_pauli': 	'Total Pauli repulsion term between substrate and catalyst fragments, unit=[ha]',
	'EDA_bonding': 	'Bonding energy difference between substrate and catalyst fragments, unit=[ha]',
	'EDA_elstat': 	'Electrostatic repulsion term between substrate and catalyst fragments, unit=[ha]',
	'EDA_oi': 		'Orbital interaction term between substrate and catalyst fragments, unit=[ha]',
	'AA_mulliken': 	'Active atom Mulliken charge on substrate/catalyst complex, unit=[a.u.]',
	'AA_eldens': 	'Active atom electron density catalyst complex, unit=[a.u.]',
	'AA_hirshfeld': 'Active atom Hirshfeld charge on substrate/catalyst complex, unit=[a.u.]',
	'AA_voronoi': 	'Active atom Voronoi charge on substrate/catalyst complex, unit=[a.u.]',
	'HOMO_coeff': 	'Highest contribution of active atom C[2pz] AO on occupied MOs',
	'HOMO_energy': 	'Energy of occupied MO with highest active atom C[2pz] contribution, unit=[ha]',
	'LUMO_coeff': 	'Highest contribution of active atom C[2pz] AO on virtual MOs',
	'LUMO_energy': 	'Energy of virtual MO with highest active atom C[2pz] contribution, unit=[ha]',
	}[parameter]

	return explanation


def define_parameters():
	#comment out the param you dont want to use
	parameters = [
		'EDA_pauli',
		'EDA_bonding',
		'EDA_elstat',
		'EDA_oi',
		'AA_mulliken',
		'AA_eldens',
		'AA_hirshfeld',
		'AA_voronoi',
		'HOMO_coeff',
		'HOMO_energy',
		'LUMO_coeff',
		'LUMO_energy']

	return parameters


def outlier_idx(columns, z_thresh=2.5):
	oi = []
	n = columns.shape[1] #number of variables
	for i in range(n):
		zs = np.abs(stats.zscore(columns[:,i]))
		oi = oi + list(np.where(zs >= z_thresh)[0])
	return list(set(oi))


def generate_data(outfile=paths.training_data, Eact_type='gibbs'):
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

		elif parameter == 'HOMO_coeff':
			return [r.data['SP']['occupied cont'] for r in results]
		elif parameter == 'HOMO_energy':
			return [r.data['SP']['occupied energy'] for r in results]
		elif parameter == 'LUMO_coeff':
			return [r.data['SP']['virtual cont'] for r in results]
		elif parameter == 'LUMO_energy':
			return [r.data['SP']['virtual energy'] for r in results]

	rxns = reaction.all_reactions
	#we now look at only achiral catalysts
	rxns = [r for r in rxns if r.reaction == 'achiral_catalyst']
	# for r in rxns:
	# 	print(r)
	Eact = [r.get_activation_energy(type='bond')[''] for r in rxns]
	Gact = [r.get_activation_energy(type='gibbs')[''] for r in rxns]
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
			# print('\t\t', type(e), e)

	meta = ['path', 'Eact', 'Gact']
	meta_data = {
		'path': respaths,
		'Eact': Eact,
		'Gact': Gact
	}

	with open(outfile, 'w+', newline='') as file:
		writer = csv.writer(file)
		writer.writerow(['sep=,'])
		writer.writerow(meta + parameters)

		for i in range(len(sub_res)):
			try:
				row = []
				for p in meta:
					row.append(meta_data[p][i])
				for p in parameters:
					row.append(data[p][i])
				writer.writerow(row)
			except Exception as e:
				print(e)
				continue


def read_data(infile=paths.training_data):
	data = []
	with open(infile, newline='') as file:
		next(file)
		reader = csv.DictReader(file)
		for row in reader:
			data.append(row)
	return data


def plot_trend(param, colors=None, labels=None, 
		xunit='', yunit='kcal/mol', 
		xscale=1, yscale=h2k(1),
		fit=True):

	used_labels = set()
	X = get_column(param)
	Y = get_column('Eact')
	
	# oi = outlier_idx(Y) + outlier_idx(X)
	# X = np.delete(X, oi)
	# Y = np.delete(Y, oi)
	# colors = np.delete(colors, oi, axis=0)
	# labels = np.delete(labels, oi)

	# X, X_norm_data = MLutils.prepare_data(X, False, True)
	# Y, Y_norm_data = MLutils.prepare_data(Y)

	# X = MLutils.unprepare_data(X, X_norm_data)

	for i in range(X.size):
		x = X[i] * xscale
		y = Y[i] * yscale
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

	if xunit != '':
		plt.xlabel(f'{param} ({xunit})')
	else:
		plt.xlabel(f'{param}')
	plt.ylabel(f'Eact ({yunit})')
	

def get_results():
	respaths = [d['path'] for d in data]
	res = [r for r in job_results3.all_results if os.path.relpath(r.path, paths.results) in respaths]
	res = sorted(res, key=lambda r: respaths.index(os.path.relpath(r.path, paths.results)))
	return res


def get_column(param):
	c = np.array([float(d[param]) for d in data])
	return c.reshape(-1,1)

def get_all_columns():
	'''
	Return a matrix with m rows and n columns 
	for m datapoints and n dependent variables
	'''
	ps = define_parameters()
	columns = []
	for p in ps:
		column = get_column(p)
		columns.append(column)
	
	c = np.array([[float(d[p]) for d in data] for p in ps]).T
	return c


def tf_dataset():
	import tensorflow.data as tfdata
	print(data)


generate_data()
data = read_data()

if __name__ == '__main__':
	print('Explanation of parameters:')
	pl = max(len(p) for p in define_parameters())
	for p in define_parameters():
		print(f'\t{p.ljust(pl+4)} {help(p)}')

	cmap = cm.get_cmap('Set1')
	
	res = get_results()
	# all_R1s = list(set(r.substituents['R1'] for r in res))
	all_R2s = list(set(r.substituents['R2'] for r in res))
	# colors = cmap([all_R1s.index(r.substituents['R1'])/len(all_R1s) for r in res])
	colors = cmap([all_R2s.index(r.substituents['R2'])/len(all_R2s) for r in res])
	# labels = [r.substituents['R1'] for r in res]
	labels = [r.substituents['R2'] for r in res]

	plt.subplot(2,2,1)
	plot_trend('EDA_oi',  colors=colors, labels=labels, xunit='Ha')
	plt.subplot(2,2,2)
	plot_trend('EDA_pauli', colors=colors, labels=labels, xunit='Ha', xscale=1)
	plt.subplot(2,2,3)
	plot_trend('EDA_elstat',  colors=colors, labels=labels, xunit='Ha')
	plt.subplot(2,2,4)
	plot_trend('EDA_bonding', colors=colors, labels=labels, xunit='Ha', xscale=1)
	plt.legend()
	plt.show()

	plt.subplot(2,2,1)
	plot_trend('HOMO_coeff',  colors=colors, labels=labels)
	plt.subplot(2,2,2)
	plot_trend('HOMO_energy', colors=colors, labels=labels, xunit='eV', xscale=1)
	plt.subplot(2,2,3)
	plot_trend('LUMO_coeff',  colors=colors, labels=labels)
	plt.subplot(2,2,4)
	plot_trend('LUMO_energy', colors=colors, labels=labels, xunit='eV', xscale=1)
	plt.legend()
	plt.show()

	plt.subplot(2,2,1)
	plot_trend('AA_mulliken', colors=colors, labels=labels, xunit='a.u.')
	plt.subplot(2,2,2)
	plot_trend('AA_eldens',   colors=colors, labels=labels, xunit='a.u.')
	plt.subplot(2,2,3)
	plot_trend('AA_hirshfeld',colors=colors, labels=labels, xunit='a.u.')
	plt.subplot(2,2,4)
	plot_trend('AA_voronoi',  colors=colors, labels=labels, xunit='a.u.')
	plt.legend()
	plt.show()
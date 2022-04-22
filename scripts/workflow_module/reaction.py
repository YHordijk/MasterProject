import os, paths
import workflow_module.results as wfr
import workflow_module.stationary_point as wfsp
import struct_generator2
import scm.plams as plams


join = os.path.join
ls = os.listdir
isdir = os.path.isdir
exists = os.path.exists



def get_reaction(calc_dir=paths.calculations2, **kwargs):
	settings = {
	'reaction' 			: kwargs.get('reaction'),
	'substituents' 		: kwargs.get('substituents'),
	'functional' 		: kwargs.get('functional', 'BLYP-D3(BJ)'),
	'basis' 			: kwargs.get('basis', 'TZ2P'),
	'numerical_quality' : kwargs.get('quality', 'Good'),
	'frozen_core' 		: kwargs.get('frozen_core', 'None'),
	'phase' 			: kwargs.get('phase', 'Dichloromethane'),
	}

	matched = [r for r in wfr.results if r.match(**settings)]
	stationary_points = {sp: [] for sp in set(match['stationary_point'] for match in matched)}
	[stationary_points[match['stationary_point']].append(match) for match in matched]
	stationary_points = [wfsp.StationaryPoint(rs) for rs in stationary_points.values()]
	return Reaction(stationary_points)
	# #sort into stationary points
	# res_for_sp = {sp}
	# for match in matched:
	
	# for match in matched:
	# 	print(match['task'], match['stationary_point'])


class Reaction:
	'''
	Container class for wfsp.StationaryPoint objects
	each reaction holds a number of these objects and interfaces them
	'''
	def __init__(self, stationary_points=[]):
		self.stationary_points = stationary_points
		self.stationary_points_dict = {sp['stationary_point']: sp for sp in stationary_points}

	def __repr__(self):
		s = f'{self["name"]} [{", ".join(sp["stationary_point"] for sp in self.stationary_points)}]'
		return s

	def __getitem__(self, key):
		if key in ['name', 'phase']:
			return self.stationary_points[0][key]
		if key in ['substituents']:
			subs = {}
			for sp in self.stationary_points:
				for R, sub in sp['substituents'].items():
					if not R in subs:
						subs[R] = sub
			return subs

	def activation_energy(self, energy_type='gibbs_COSMORS', unit='kcal/mol'):
		'''
		gets the reaction energy as defined in the meta.info of the reaction template
		'''

		meta_info = struct_generator2.get_meta_info(self['name'])
		order = [m  for m in meta_info if m[0] == 'activation_energy'][0][1:]
		sign = [{'+':1, '-':-1}[o[0]] for o in order]
		sps = [o[1:] for o in order]
		energies = [self.stationary_points_dict[sp][energy_type] for sp in sps]

		return plams.Units.convert(sum(si*en for si, en in zip(sign, energies)), 'ha', unit)

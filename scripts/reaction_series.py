import reaction 
import itertools
import matplotlib.pyplot as plt

class Series:
	def __init__(self, template, substituents, functional='BLYP-D3(BJ)', basis='TZ2P', quality='Good', phase='vacuum'):
		self.template = template
		self.substituents = substituents
		self.functional = functional
		self.basis = basis
		self.quality = quality
		self.phase = phase

		self.identifiers = self.get_identifiers()
		self.variable = self.get_variable()
		assert len(self.variable) > 1, 'Too many variables!'
		self.reactions = self.get_reactions()

	def get_identifiers(self):
		ide = {}
		ide['template'] = self.template 
		for R, s in self.substituents.items():
			ide[R] = s
		ide['functional'] = self.functional
		ide['basis'] = self.basis
		ide['quality'] = self.quality
		ide['phase'] = self.phase
		return ide

	def get_variable(self):
		for i, v in self.identifiers.items():
			if type(v) is list:
				return i

	def get_reactions(self):
		#get product first
		reactions = {}
		prod = itertools.product(*[[[v],v][i == self.variable] for i, v in self.identifiers.items()])
		for p in prod:
			p_dict = {i:x for i, x in zip(self.identifiers.keys(), p)}
			reactions[p] = reaction.Reaction(p_dict['template'], {R:s for R, s in p_dict.items() if R.startswith('R')}, functional=p_dict['functional'], basis=p_dict['basis'], numerical_quality=p_dict['quality'])
		return reactions

	def plot_activation_energy(self, enantiomer=''):
		Eact = {p: r.get_activation_energy()[enantiomer] for p, r in self.reactions.items()}
		main_keys = self.identifiers[self.variable]
		main_Eact = {v:[] for v in main_keys}
		for p, E in Eact.items():
			p_dict = {i:x for i, x in zip(self.identifiers.keys(), p)}
			main_Eact[p_dict[self.variable]] = E

		plt.plot([])
		plt.gca().set_xticks(range(len(self.identifiers[self.variable])))
		plt.gca().set_xticklabels(self.identifiers[self.variable])
		plt.legend()
		




if __name__ == '__main__':
	s = Series('achiral_catalyst', {'R1':'H', 'R2':'Ph', 'Rcat':['AlF3', 'BF3', 'I2', 'SnCl4', 'TiCl4', 'ZnCl2']})
	s.plot_activation_energy()
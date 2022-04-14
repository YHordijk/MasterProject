import scm.plams as plams



class Workflow:
	def __init__(self):
		...


	def __call__(self, mol, rundir, **kwargs):
		workflow_settings = self.parse_settings(**kwargs)

		res1 = self.optimize_mol(mol, workflow_settings)
		res2 = self.


	def parse_settings(**kwargs):
		settings = {}
		settings['solvent']     = kwargs.get('solvent', 'Dichloromethane')
		settings['functional']  = kwargs.get('functional', 'BLYP-D3(BJ)')
		settings['basis']       = kwargs.get('basis', 'TZ2P')
		settings['quality']     = kwargs.get('quality', 'Good')
		settings['relativity']  = kwargs.get('relativity', 'Scalar')
		settings['frozen core'] = kwargs.get('frozen core', 'None')

		return settings



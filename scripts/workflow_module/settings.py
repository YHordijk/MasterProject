import scm.plams as plams
import numpy as np
import paths, os



def _get_COSMO_radii(preset='Klamt'):
	radii = {}
	with open(os.path.join(paths.COSMO_radii, preset)) as ref:
		for line in ref.readlines():
			a, r = line.strip().split(',')
			radii[a] = float(r)

	return radii


def get_TSRC_preopt_settings(**kwargs):
	s = get_opt_settings(**kwargs)
	s.input.ams.task = 'GeometryOptimization'
	s.input.ams.TransitionStateSearch = None
	a, b = kwargs.get('TSRC idx')
	mol = kwargs.get('mol')
	aa, ab = mol.atoms[a-1], mol.atoms[b-1]
	d = np.linalg.norm(np.array(aa.coords) - np.array(ab.coords))
	
	s.input.ams.Constraints.Distance = f'{a} {b} {d}'

	return s


def get_opt_settings(**kwargs):
	s = plams.Settings()

	if kwargs.get('opt task') == 'TSRC':
		s.input.ams.task = 'TransitionStateSearch'
		TSRC_idx = kwargs.get('TSRC idx')
		s.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = f'{TSRC_idx[0]} {TSRC_idx[1]} 1.0'
	else:
		s.input.ams.task = 'GeometryOptimization'

	if kwargs.get('functional') == 'BLYP-D3(BJ)':
		s.input.adf.XC.GGA = 'BLYP'
		s.input.adf.XC.Dispersion = 'GRIMME3 BJDAMP'

	# s.input.ams.Properties.PESPointCharacter = 'Yes'
	s.input.adf.basis.Type = kwargs.get('basis')
	s.input.adf.basis.Core = kwargs.get('frozen core')
	s.input.adf.NumericalQuality = kwargs.get('quality')
	s.input.adf.SYMMETRY = 'NOSYM'

	if kwargs.get('radical'):
		s.input.adf.SpinPolarization = kwargs.get('spin polarization', 1.0)
		s.input.adf.Unrestricted = kwargs.get('unrestricted', 'Yes')

	if kwargs.get('solvent') is not None:
		atoms = kwargs.get('mol').atoms
		radii = _get_COSMO_radii('Klamt')

		solv = kwargs.get('solvent')
		solvation_block = {
		    'surf': 'Delley',
		    'solv': f'name={solv} cav0=0.0 cav1=0.0067639',
		    'charged': 'method=Conj',
		    'c-mat': 'POT',
		    'scf': 'Var All',
		    # 'CSMRSP': '',
		    'radii': {a.symbol: radii[a.symbol] for a in atoms}}

		s.input.adf.solvation = solvation_block

	return s


def get_COSMORS_settings(**kwargs):
	s = plams.Settings()
	s.input.property._h = 'ACTIVITYCOEF'
	compounds = [plams.Settings(), plams.Settings()]
	solvent_path = os.path.join(paths.COSMO_rkf, kwargs['solvent'] + '.coskf')
	compounds[0]._h = solvent_path
	compounds[1]._h = kwargs['solute_path']
	compounds[0].frac1 = 1
	compounds[1].frac1 = 0
	s.input.temperature = kwargs.get('temperature', '298.15')
	s.input.compound = compounds
	return s


def get_FREQ_settings(**kwargs):
	s = get_opt_settings(**kwargs)
	s.input.ams.task = 'SinglePoint'
	s.input.ams.TransitionStateSearch = None
	s.input.ams.Properties.NormalModes = 'Yes'
	s.input.ams.NormalModes.ReScanFreqRange = '-1000000000.0 20.0'
	s.input.ams.PESPointCharacter.NegativeFrequenciesTolerance = '-30'
	return s
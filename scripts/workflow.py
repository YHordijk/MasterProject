import scm.plams as plams
import workflow_module.settings as wfs
import workflow_module.results as wfr
import workflow_module.reaction as wfreact
import workflow_module.stationary_point as wfsp
import paths, struct_generator2, os, utility, sys
import workflow_module.task as wft
import matplotlib.pyplot as plt
import utility
import numpy as np



def sort_substituents(subs):
	sorted_Rnames = list(sorted(subs))
	sorted_R = [subs[R] for R in sorted_Rnames]
	return sorted_R



class Workflow:
	def __init__(self, 
				rerun_calculations=False, 
				calc_dir=paths.calculations2,
				skip_TS_preopt=False,
				frequencies=False,
				**kwargs):
		'''
		This workflow performs a few steps
		1. GO or TSRC for all stationary points
		1a. For TSRC perform a contrained geometry optimization beforehand to improve stability
		2. COSMORS calculation for all stationary points
		3. EDA for substrate complex
		4. aligned SP for substrate complex
		5. optionally frequency calculations for all mols
		
		The filetree will look like this
		calc_dir
		|-- {reaction}.{sorted_substituents}
		|-- -- {stationary_point}
		|-- -- -- {task_name}.{index}
		|-- -- -- -- job.hash
		|-- -- -- -- job.info
		|-- -- -- -- job.input
		|-- -- -- -- SUCCESS #if job was successful
		|-- -- -- -- input.xyz
		|-- -- -- -- output.xyz
		'''
		self.workflow_settings = self.parse_settings(**kwargs)
		self.workflow_settings['rerun calculations'] = rerun_calculations
		self.workflow_settings['calc dir'] = calc_dir
		self.workflow_settings['skip_TS_preopt'] = skip_TS_preopt
		self.workflow_settings['frequencies'] = frequencies


	def run(self, **kwargs):
		name = f'{self.workflow_settings["reaction"]}.{"_".join(sort_substituents(self.workflow_settings["substituents"]))}'
		rundir = os.path.join(self.workflow_settings['calc dir'], name)
		self.workflow_settings['run dir'] = rundir

		for mol, mol_settings in self.setup(rundir):
			#we will preoptimize the TS using a constrained GO calculation to let the system relax a little
			res_GO 		= self.optimize_mol(mol, **mol_settings)
			### CALCULATE THE ACTIVATION COEFFICIENT WHICH ALSO CALCULATES THE GIBBS FREE ENERGY
			res_COSMORS = self.COSMORS(res_GO.get_solute_path(), **mol_settings)
			### OPTIONALLY CALCULATE FREQUENCIES
			if self.workflow_settings['frequencies']:
				res_freq    = self.frequencies(res_GO.get_main_molecule(), **mol_settings)



	def setup(self, rundir):
		os.makedirs(rundir, exist_ok=True)

		mols = struct_generator2.generate_stationary_points(self.workflow_settings['reaction'], self.workflow_settings['substituents'])
		mol_settings = []
		for name, mol in mols.items():
			s = {}
			s['stationary point'] = name
			s['radical'] = mol.radical
			s['enantiomer'] = mol.enantiomer
			s['TSRC idx'] = mol.TSRC_idx
			s['opt task'] = mol.task
			if mol.get_frags:
				s['frag1 idx'] = mol.frag1idx
				s['frag2 idx'] = mol.frag2idx
			if hasattr(mol, 'active_atom_idx'):
				s['active atom idx'] = mol.active_atom_idx
			else:
				s['active atom idx'] = None
			s['plane idx'] = mol.plane_idx
			s['align idx'] = mol.align_idx
			s['center idx'] = mol.center_idx
			mol_settings.append(s)

		xyz_already_exists = {name: os.path.exists(os.path.join(rundir, name + '.xyz')) for name in struct_generator2.get_all_stationary_points(self.workflow_settings['reaction'])}
		[mol.write(os.path.join(rundir, name + '.xyz')) for name, mol in mols.items() if not xyz_already_exists[name]]
		mols = [plams.Molecule(os.path.join(rundir, name + '.xyz')) for name in mols.keys()]
		return zip(mols, mol_settings)


	def parse_settings(self, **kwargs):
		settings = {}
		settings['reaction']     = kwargs.get('reaction', None) #string
		settings['substituents'] = kwargs.get('substituents', None) #dictionary string[string]
		settings['solvent']      = kwargs.get('solvent', 'Dichloromethane') #string
		settings['functional']   = kwargs.get('functional', 'BLYP-D3(BJ)') #string
		settings['basis']        = kwargs.get('basis', 'TZ2P') #string
		settings['quality']      = kwargs.get('quality', 'Good') #string
		settings['relativity']   = kwargs.get('relativity', 'Scalar') #string
		settings['frozen core']  = kwargs.get('frozen_core', 'None') #string
		return settings


	def optimize_mol(self, mol, **mol_settings):
		### FIRST OPTIMIZE THE MOLECULE
		### IF WE ARE DOING TSRC THEN WE WANT TO PRE-OPTIMIZE USING CONSTRAINED GEO-OPT
		if mol_settings['opt task'] == 'TSRC' and not self.workflow_settings['skip_TS_preopt']:
			s = wfs.get_TSRC_preopt_settings(mol=mol, **self.workflow_settings, **mol_settings)
			res = wft.Optimization(mol, s, 'constrained_GO', **self.workflow_settings, **mol_settings).run()
			mol = res.get_main_molecule()

		### DO THE HIGH-LEVEL GEO-OPT/TSRC
		s = wfs.get_opt_settings(mol=mol, **self.workflow_settings, **mol_settings)
		res = wft.Optimization(mol, s, 'GO', **self.workflow_settings, **mol_settings).run()
		return res


	def frequencies(self, mol, **mol_settings):
		s = wfs.get_FREQ_settings(mol=mol, **self.workflow_settings, **mol_settings)
		res = wft.Frequencies(mol, s, 'FREQ', **self.workflow_settings, **mol_settings).run()
		return res


	def COSMORS(self, solute_path, **mol_settings):
		s = wfs.get_COSMORS_settings(solute_path=solute_path, **self.workflow_settings, **mol_settings)
		res = wft.COSMORS(s, 'COSMORS', **self.workflow_settings, **mol_settings).run()
		return res


if __name__ == '__main__':
	if len(sys.argv) > 1:
		arg_dict 			= {a.split('=')[0][1:].strip(): a.split('=')[1].strip() for a in sys.argv[1:]}

		reaction 			= arg_dict.get('reaction', 'achiral_radical_addition')
		substituents 		= {R:sub for R, sub in arg_dict.items() if R.startswith('R')}
		basis 				= arg_dict.get('basis', 'TZ2P')
		functional 			= arg_dict.get('functional', 'BLYP-D3(BJ)')
		frozen_core 		= arg_dict.get('frozen_core', 'None')
		quality 			= arg_dict.get('quality', 'Good')
		frequencies 		= bool(arg_dict.get('frequencies', False))
		skip_TS_preopt 		= bool(arg_dict.get('skip_TS_preopt', False))
		rerun_calculations 	= bool(arg_dict.get('rerun_calculations', False))
		relativity 			= arg_dict.get('relativity', 'Scalar')
		solvent 			= arg_dict.get('solvent', 'Dichloromethane')

		workflow = Workflow(reaction=reaction, solvent=solvent,
						    substituents=substituents, relativity=relativity,
						    basis=basis, frozen_core=frozen_core, quality=quality,
							frequencies=frequencies, skip_TS_preopt=skip_TS_preopt, rerun_calculations=rerun_calculations)
		workflow.run()


	else:
		# workflow = Workflow(reaction='achiral_radical_addition', 
		# 				    substituents={'R1':'p-HOPh', 'R2':'tBu', 'Rcat':'ZnCl2'}, 
		# 				    basis='TZ2P', frozen_core='None', quality='Good',
		# 					frequencies=True, skip_TS_preopt=False, rerun_calculations=False)
		# workflow.run()

		
		# r = wfreact.get_reaction(reaction='achiral_radical_addition', 
		# 				 substituents={'R1':'Me', 'R2':'Et', 'Rcat':'AlF3'}, 
		# 				 basis='TZ2P', frozen_core='None', quality='Good')

		# print(r.activation_energy())
		# print(r.activation_energy('gibbs'))


		# join = os.path.join


		# c = paths.calculations2
		# jobs = []
		# for rundir in os.listdir(c):
		# 	# print(rundir)
		# 	job = []
		# 	subs = rundir.split('.')[1]
		# 	subs = {'R1':subs.split('_')[0], 'R2':subs.split('_')[1], 'Rcat':subs.split('_')[2]}
		# 	job.append(subs)
		# 	for spdir in os.listdir(join(c, rundir)):
		# 		if not os.path.isdir(join(c, rundir, spdir)): continue
		# 		for jobdir in os.listdir(join(c, rundir, spdir)):

		# 			with open(join(c, rundir, spdir, jobdir, 'job.info')) as info:
		# 				lines = info.readlines()
		# 				for line in lines:
		# 					if line.startswith('phase'):
		# 						phase = line.split('=')[1].strip()
		# 						if not phase in job:
		# 							job.append(phase)
		# 						break
		# 	jobs.append(job)

		# for job in jobs:
		# 	for solvent in job[1:]:
		# 		print(job[0], solvent)
		# 		workflow = Workflow(reaction='achiral_radical_addition', solvent=solvent,
		# 						    substituents={'R1':job[0]['R1'], 'R2':job[0]['R2'], 'Rcat':job[0]['Rcat']}, 
		# 						    basis='TZ2P', frozen_core='None', quality='Good',
		# 							frequencies=False, skip_TS_preopt=False, rerun_calculations=False)
		# 		workflow.run()




		match_dict = lambda a, b: all(a[k] == b.get(k, None) for k in a.keys())

		Gwater = np.array([sp['gibbs'] for sp in wfsp.stationary_points if sp['phase'] == 'Water'])
		Gcosmowater = np.array([sp['gibbs_COSMORS'] for sp in wfsp.stationary_points if sp['phase'] == 'Water'])
		Gdcm = np.array([sp['gibbs'] for sp in wfsp.stationary_points if sp['phase'] == 'Dichloromethane'])
		Gcosmodcm = np.array([sp['gibbs_COSMORS'] for sp in wfsp.stationary_points if sp['phase'] == 'Dichloromethane'])


		
		validwater = []
		for c, s in zip(Gcosmowater, Gwater):
			if not c is None and not s is None:
				validwater.append((c, s))
		validwater = np.array(validwater) * utility.hartree2kcalmol(1)


		validdcm = []
		for c, s in zip(Gcosmodcm, Gdcm):
			if not c is None and not s is None:
				validdcm.append((c, s))
		validdcm = np.array(validdcm) * utility.hartree2kcalmol(1)

		ma = max(validwater.max(), validdcm.max())*0.9
		mi = min(validwater.min(), validdcm.min())*1.1

		plt.subplot(1,2,1)
		plt.suptitle('Statistical analysis vs COSMO-RS for single molecules')
		plt.gca().set_title('Correlation')
		plt.scatter(validwater[:,0], validwater[:,1], label='Water')
		plt.scatter(validdcm[:,0], validdcm[:,1],   label='Dichloromethane')
		plt.plot((mi,ma), (mi,ma), 'k--')
		plt.gca().set_xlim(mi, ma)
		plt.gca().set_ylim(mi, ma)
		plt.legend()
		plt.xlabel(r'$\Delta G_{ST}$ (kcal/mol)')
		plt.ylabel(r'$\Delta G_{CRS}$ (kcal/mol)')
	
		plt.subplot(1,2,2)
		plt.gca().set_title('Residuals')
		plt.hist(validwater[:,0]-validwater[:,1], label='Water', fc='blue', alpha=0.5, density=False)
		plt.hist(validdcm[:,0]-validdcm[:,1], label='Dichloromethane', fc='orange', alpha=0.5, density=False)
		plt.xlabel(r'$\Delta\Delta G$ (kcal/mol)')
		plt.ylabel('Count')
		plt.legend()
		plt.show()


		reactions = []
		for sp in wfsp.stationary_points:
			reactions.append(wfreact.get_reaction(reaction=sp['reaction'], substituents=sp['substituents'], phase=sp['phase']))

		reactions = set(reactions)

		datawater = []
		datadcm = []
		for r in reactions:
			if r['phase'] == 'Water':
				try:
					datawater.append((r.activation_energy('gibbs_COSMORS'), r.activation_energy('gibbs')))
				except:
					pass
			elif r['phase'] == 'Dichloromethane':
				try:
					datadcm.append((r.activation_energy('gibbs_COSMORS'), r.activation_energy('gibbs')))
				except:
					pass

		datawater = np.array(datawater)
		datadcm = np.array(datadcm)

		ma = max(datawater.max(), datadcm.max())*1.1
		mi = min(datawater.min(), datadcm.min())*1.1

		plt.subplot(1,2,1)
		plt.gca().set_title('Correlation')
		plt.suptitle('Statistical analysis vs COSMO-RS for activation energies')
		plt.scatter(datawater[:,1], datawater[:,0], label='Water')
		plt.scatter(datadcm[:,1], datadcm[:,0], label='Dichloromethane')
		plt.plot((mi,ma), (mi,ma), 'k--')
		plt.legend()
		plt.xlabel(r'$\Delta G^‡_{ST}$ (kcal/mol)')
		plt.ylabel(r'$\Delta G^‡_{CRS}$ (kcal/mol)')
		plt.gca().set_xlim(mi, ma)
		plt.gca().set_ylim(mi, ma)


		plt.subplot(1,2,2)
		plt.gca().set_title('Residuals')
		plt.hist(datawater[:,1]-datawater[:,0], label='Water', fc='blue', alpha=0.5, density=False)
		plt.hist(datadcm[:,1]-datadcm[:,0], label='Dichloromethane', fc='orange', alpha=0.5, density=False)
		plt.xlabel(r'$\Delta\Delta G^‡$ (kcal/mol)')
		plt.ylabel('Count')
		plt.legend()
		plt.show()
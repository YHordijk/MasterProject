import struct_generator2, paths, os, utility
import job_runner.JRutils as JRutils
import scm.plams as plams

join = os.path.join



def GO(calc_path, phase='DCM', functional='BLYP-D3(BJ)', basis='TZ2P', numerical_quality='Good'):

	assert basis in ['SZ', 'DZ', 'DZP', 'TZP', 'TZ2P', 'QZ4P']
	assert functional in ['OLYP', 'BLYP-D3(BJ)', 'M06L']
	assert numerical_quality in ['Basic', 'Normal', 'Good', 'VeryGood', 'Excellent']


	def get_task():
		if mol_info['TSRC_idx'] == 'N/A':
			return 'GO'
		else:
			return 'TSRC'

	def write_run_info(run_dir):
		with open(join(run_dir, 'run.info'), 'w+') as info:
			info.write(f'task={get_task()}\n')
			info.write(f'reaction={mol_info["reaction"]}\n')
			info.write(f'sorted_substituents={mol_info["sorted_substituents"]}\n')
			info.write(f'stationary_point={mol_info["stationary_point"]}\n')
			info.write(f'TSRC_idx={mol_info["TSRC_idx"]}\n')
			info.write(f'radical={mol_info["radical"]}\n')
			info.write(f'functional={functional}\n')
			info.write(f'basis={basis}\n')
			info.write(f'numerical_quality={numerical_quality}\n')
			info.write(f'phase={phase}\n')
			info.write(f'unique={JRutils.get_unique(mol_info)}\n')

	def runscript_path(run_dir):
		return join(run_dir, f'{get_task()}.{mol_info["reaction"]}.{mol_info["sorted_substituents"]}.{mol_info["stationary_point"]}')

	def write_runscript(run_dir):
		RADICAL = ''
		if mol_info['radical']:
			RADICAL = '\n  SpinPolarization 1\n  Unrestricted Yes\n'

		mol = plams.Molecule(join(calc_path, 'init.xyz'))
		ATOMS = ''
		for atom in mol.atoms:
			ATOMS += f'    {atom.symbol:<2}\t{atom.coords[0]:=11.8f}\t{atom.coords[1]:=11.8f}\t{atom.coords[2]:=11.8f}\n'

		TSRC = ''
		if get_task() == 'TSRC':
			TSRC_idx = [int(i) for i in mol_info['TSRC_idx'].split('_')]
			TSRC = f'    Distance {TSRC_idx[0]} {TSRC_idx[1]} -1\n'
		SLURM_SUBMIT_DIR = runscript_path(run_dir)

		CORES = '64'
		NODES = '1'
		if mol_info['radical'] and mol_info['reaction'] in ['urea_tBu_Ph', 'squaramide']:
			CORES = '64'
			NODES = '2'

		FUNCTIONAL = {
			'OLYP': 'GGA OLYP',
			'BLYP-D3(BJ)': 'GGA BLYP\n    DISPERSION GRIMME3 BJDAMP',
			'M06L': 'MetaGGA M06L'
		}[functional]


		blocks = {
				'[RADICAL]':          RADICAL, 
				'[ATOMS]':            ATOMS, 
				'[TSRC]':             TSRC, 
				'[SLURM_SUBMIT_DIR]': SLURM_SUBMIT_DIR,
				'[CORES]':            CORES,
				'[NODES]':            NODES,
				'[BASIS]':            basis,
				'[FUNCTIONAL]':       FUNCTIONAL,
				'[NUMERICALQUALITY]': numerical_quality,
				}

		with open(JR_template) as temp:
			temp_lines = temp.readlines()
			with open(SLURM_SUBMIT_DIR, 'w+') as run:
				for line in temp_lines:
					for block in blocks:
						if block in line:
							line = line.replace(block, blocks[block])

					run.write(line)


	mol_info = JRutils.read_mol_info(calc_path)
	mol_info['task'] = get_task()
	mol_info['phase'] = phase
	mol_info['functional'] = functional
	mol_info['basis'] = basis
	mol_info['numerical_quality'] = numerical_quality

	GO_run_dir = JRutils.get_run_dir(calc_path, 'GO', mol_info)
	if GO_run_dir is None:
		print('collision')
		return

	os.makedirs(GO_run_dir, exist_ok=True)
	JR_template = JRutils.get_job_runner_template(get_task())
	write_run_info(GO_run_dir)
	write_runscript(GO_run_dir)

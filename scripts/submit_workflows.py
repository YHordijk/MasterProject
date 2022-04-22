
import os, sys, paths
import numpy as np


job_run_location = os.path.join(paths.scripts, 'workflow_runscripts')

Rs = ['Bz', 'Et', 'Me', 'm-FPh', 'o-FPh', 'p-FPh', 'm-HOPh', 'o-HOPh', 'm-HOPh', 'NH2', 'NHMe', 'NMe2', 'OEt', 'OMe', 'Ph', 'tBu']
cats = ['AlF3', 'I2', 'BF3', 'SnCl4', 'TiCl4', 'ZnCl2']

for i in range(15):
	settings = {
	'reaction': 'achiral_radical_addition',
	'R1': np.random.choice(Rs),
	'R2': np.random.choice(Rs),
	'Rcat': np.random.choice(cats),
	'frequencies': True,
	'solvent': 'Water',
	}

	args = ' '.join([f'-{key}={str(value)}' for key, value in settings.items()])
	print(args)

	runscript = rf'''#!/bin/bash
#SBATCH -N 2
#SBATCH -t 240:00:00
#SBATCH --ntasks-per-node=64
#SBATCH --partition=tc

echo "== Starting workflow at $(date)"
echo "== Job ID: ${{SLURM_JOBID}}"
echo "== Node list: ${{SLURM_NODELIST}}"
echo "== Submit dir. : ${{SLURM_SUBMIT_DIR}}"
echo "== Scratch dir. : ${{TMPDIR}}"

module load ams
module list

$AMSBIN/amspython {paths.scripts}/workflow.py {args}
	'''

	file = os.path.join(job_run_location, f'{settings["reaction"]}.{".".join([settings[R] for R in settings.keys() if R.startswith("R")])}.{settings["solvent"]}')
	with open(file, 'w+') as r:
		r.write(runscript)

	start_script = f"""#!/bin/bash
cd {job_run_location}
sbatch {file}
"""     
	# print(start_script)
	os.system(start_script)
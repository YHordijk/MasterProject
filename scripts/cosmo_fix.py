
import os, sys, paths
join = os.path.join


c = paths.calculations2


for rundir in os.listdir(c):
	for spdir in os.listdir(join(c, rundir)):
		if not os.path.isdir(join(c, rundir, spdir)): continue
		for jobdir in os.listdir(join(c, rundir, spdir)):
			if not jobdir.startswith('COSMORS'): continue
			with open(join(c, rundir, spdir, jobdir, 'COSMORS', 'COSMORS.run'), 'w+') as run:
				run.write('''#!/bin/sh
#SBATCH --partition=tc

$AMSBIN/crs <"COSMORS.in" > COSMORS.out
''')
			with open(join(c, rundir, spdir, jobdir, 'COSMORS', 'COSMORS.in')) as inf:
				lines = inf.readlines()
				lines.pop()
				lines.pop()
				lines.append('temperature 298.15\n')
				print(lines)

			with open(join(c, rundir, spdir, jobdir, 'COSMORS', 'COSMORS.in'), 'w+') as inf:
				inf.writelines(lines)

			file = join(c, rundir, spdir, jobdir, 'COSMORS', 'COSMORS.run')
			start_script = f"""#!/bin/bash
cd {os.path.dirname(file)}
sbatch {os.path.basename(file)}
"""     
			print(start_script)
			os.system(start_script)



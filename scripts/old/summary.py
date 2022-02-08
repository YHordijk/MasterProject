'''
This scripts provide API for summarizing calculations and outputs summary files
'''

import os
import scm.plams as plams
import paths
import results_database as rdb

def summarize_calculation(calc):
	'''
	Method that extracts usefull data from a calculation.
	To be used in DatabaseManager class
	'''

	res_dir = os.path.join(paths.results, f'{calc["id"]}.{calc["name"]}')
	try: os.mkdir(res_dir)
	except: pass
	summ = open(os.path.join(res_dir, 'summary.txt'), 'w+')

	summ.write(f'{"="*8} SUMMARY OF CALCULATION {calc["name"]} (ID = {calc["id"]}) {"="*8}\n')
	results = {}

	kf = plams.KFFile(os.path.join(paths.master, calc['directory'], calc['adfrkf']))
	results['termination'] = kf.read('General', 'termination status')
	results['runtype'] = kf.read('General', 'runtype')

	summ.write(f'==== GENERAL\n')
	summ.write(f'Termination status : {results["termination"]}\n')
	summ.write(f'Task:                {results["runtype"]}\n')


	summ.close()

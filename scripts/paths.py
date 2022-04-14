import os, sys


master = 				os.path.dirname(os.path.dirname(__file__))
driveletter =			master.split(':')[0] + ':'

scripts = 				os.path.join(master, 'scripts')
results = 				os.path.join(master, 'results')
training_data =			os.path.join(results, 'params.csv')
calculations = 			os.path.join(master, 'calculations')
calculations2 = 		os.path.join(master, 'calculations2')
results_table = 		os.path.join(results, 'master_results_table.csv')
results_table_pretty = 	os.path.join(results, 'results_formatted.xlsx')

resources = 			os.path.join(master, 'resources')
JR_templates =			os.path.join(resources, 'job_runner_templates')
SGT = 					os.path.join(resources, 'struct_generator_templates')
SGT_substituents = 		os.path.join(SGT, 'substituents')
input_xyz = 			os.path.join(resources, 'input_xyz')
opt_xyz = 				os.path.join(resources, 'opt_xyz')

__all__ = [master, scripts, results, resources, input_xyz, results_table, calculations, calculations2, SGT, SGT_substituents, results_table_pretty, training_data]
def check_paths():
	if not all(map(os.path.exists, __all__)):
		print('[paths.py]: Warning, not all paths exist.')
		for p in __all__:
			if not os.path.exists(p):
				print(p)
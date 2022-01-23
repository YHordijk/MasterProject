import os, sys


master = 				os.path.dirname(os.path.dirname(__file__))
driveletter =			master.split(':')[0] + ':'

scripts = 				os.path.join(master, 'scripts')
results = 				os.path.join(master, 'results')
calculations = 			os.path.join(master, 'calculations')
results_table = 		os.path.join(results, 'master_results_table.csv')
results_table_pretty = 	os.path.join(results, 'results_formatted.xlsx')

id_list = 				os.path.join(scripts, 'id_list')
resources = 			os.path.join(master, 'resources')
SGT = 					os.path.join(resources, 'struct_generator_templates')
SGT_substituents = 		os.path.join(SGT, 'substituents')
input_xyz = 			os.path.join(resources, 'input_xyz')
opt_xyz = 				os.path.join(resources, 'opt_xyz')

__all__ = [master, scripts, results, resources, input_xyz, results_table, calculations, SGT, SGT_substituents, id_list, results_table_pretty]
if not all(map(os.path.exists, __all__)):
	print('[paths.py]: Warning, not all paths exist.')
	for p in __all__:
		if not os.path.exists(p):
			print(p)
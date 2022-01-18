import os, sys
from types import SimpleNamespace





master = os.path.dirname(os.path.dirname(__file__))
scripts = os.path.join(master, 'scripts')
results = os.path.join(master, 'results')
calculations = os.path.join(master, 'calculations')
results_table = os.path.join(results, 'master_results_table.csv')
resources = os.path.join(master, 'resources')
xyz = os.path.join(resources, 'xyz')

__all__ = [master, scripts, results, resources, xyz, results_table, calculations]
if not all(map(os.path.exists, __all__)):
	print('[paths.py]: Warning, not all paths exist.')
import job_results3
from regression.kernels import distance2
import numpy as np
import matplotlib.pyplot as plt
# import
'''This module checks if vibrations are correct; 
i.e. do they have largest contibution along correct axis
Do they have largest one imaginairy mode?'''


def check(result=None, path=None, draw_reaction=True):
	if result is None:
		result = job_results3.Result(path)
	mol = result.get_mol()
	nm = mol.normalmode
	if nm is None:
		reason = 'FAILED: No imaginary modes'
		return reason
	nm = np.array(nm)
	nm = nm.reshape(nm.size//3,3)
	#Distance matrix will show largest element 
	#for largest bond distance change
	D = distance2(nm, nm)
	#argmax is the scalar a1 + a2*N where a1, a2 are the atom indices of interest and N is the number of atoms
	amax = np.argmax(D)
	#decode argmax
	a1, a2 = amax%result.natoms, amax//result.natoms
	TSRC_idx = result.data['info'].get('TSRC_idx')
	if TSRC_idx is None:
		TSRC_idx = [int(i) for i in result.data['info']['TSRC'].split('_')]
	if a1+1 in TSRC_idx and a2+1 in TSRC_idx:
		return True
	else:
		reason = f'FAILED: Largest imaginary mode component on atoms {a1+1} and {a2+1} instead of expected {TSRC_idx[0]} and {TSRC_idx[1]}'
		return reason

	# print(D)
	# plt.imshow(D)
	# plt.show()



if __name__ == '__main__':
	results = job_results3.all_results
	results = job_results3.filter(results, stationary_point='TS')
	for res in results:
		try:
			if not check(res) == True:
				print('Failure in', res.path, ':')
				print(check(res))
		except:
			print(res.path)
			raise
		# print(check(res), res)
	# check(path=r"D:\Users\Yuman\Desktop\MasterProject\results\achiral_catalyst.Br_Ph_BF3.vacuum\TS")
import job_results3, os, utility
join = os.path.join


def check_EDA(result=None, path=None, draw_reaction=False, print_reason=False, suppress=[]):
	if result is None:
		result = job_results3.Result(path)
	reasons = []
	eda_out = result.data['files'].get('EDA log', None)
	if eda_out is None:
		if not 0 in suppress: reasons.append((0, 'EDA job not finished'))
	else:
		with open(join(result.calc_path, eda_out)) as log:
			for line in log.readlines():
				if 'ERROR: Engine ADF failed to solve for single point calculation.' in line:
					if not 1 in suppress: reasons.append((1, 'SCF failed to converge'))
				elif 'ERROR' in line:
					if not -1 in suppress: reasons.append((-1, f'Unknown error: {" ".join(line.split()[3:])}'))
					
	if print_reason:
		for r in reasons:
			print('\t',f'code {r[0]}: {r[1]}')

	return len(reasons) == 0


def check_SP(result=None, path=None, draw_reaction=False, print_reason=False, suppress=[]):
	if result is None:
		result = job_results3.Result(path)
	reasons = []
	SP_log = result.data['files'].get('SP log', None)
	if SP_log is None:
		if not 0 in suppress: reasons.append((0, 'SP job not finished'))
	else:
		with open(join(result.calc_path, SP_log)) as log:
			for line in log.readlines():
				if 'ERROR: Engine ADF failed to solve for single point calculation.' in line:
					if not 1 in suppress: reasons.append((1, 'SCF failed to converge'))
				elif 'ERROR' in line:
					if not -1 in suppress: reasons.append((-1, f'Unknown error: {" ".join(line.split()[3:])}'))
					
	if print_reason:
		for r in reasons:
			print('\t',f'code {r[0]}: {r[1]}')

	return len(reasons) == 0


if __name__ == '__main__':
	results = job_results3.all_results
	results = job_results3.filter(results, stationary_point='sub_cat_complex', template='achiral_catalyst')
	print('Checking EDA calculations:')
	for i, res in enumerate(results):
		utility.loading_bar(i, len(results)-1)
		if check_EDA(res) != True:
			print('Failure in', res.path, ':')
			check_EDA(res, draw_reaction=False, print_reason=True)

	print('Checking SP calculations:')
	for i, res in enumerate(results):
		utility.loading_bar(i, len(results)-1)
		if check_SP(res) != True:
			print('Failure in', res.path, ':')
			check_SP(res, draw_reaction=False, print_reason=True)

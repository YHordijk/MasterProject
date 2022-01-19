import paths

def get_free_ID():
	with open(paths.results_table, 'r') as table:
		IDs = [int(line.split(',')[0]) for line in table.readlines()]
		if len(IDs) > 0:
			for i in range(max(IDs)):
				if not i+1 in IDs:
					return i+1
		else:
			return 1
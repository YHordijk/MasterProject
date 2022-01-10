import paths


## this script implements a way to retrieve results from calculations done using AMS
## Each calculation result must have a unique ID used to retrieve it from the master results table
class AMS_results:
	def __init__(self, ID):
		self.ID = ID
		self.load_result(ID)

	def load_result(self, ID):
		with open(paths.results_table, 'r') as table:
			print(table.readlines())


if __name__ == '__main__':
	r = AMS_results(0)
	
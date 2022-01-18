import paths
import csv
import os
import summary


def create_database_from_calculations(calc_path=paths.calculations, out_path=paths.results_table):
	fields = ['id', 
	          'name', 
	          'directory', 
	          'logfile', 
	          'outfile', 
	          'adfrkf', 
	          'amsrkf']

	data = []
	dirs = []
	for file in os.listdir(calc_path):
		p = os.path.join(calc_path, file)
		if os.path.isdir(p):
			dirs.append(os.path.relpath(p, paths.master))

	for d in dirs:
		data_dict = {}

		data_dict['id'] = os.path.basename(d).split('.')[0]
		data_dict['name'] = '.'.join(os.path.basename(d).split('.')[1:])
		data_dict['directory'] = d



		for file in os.listdir(os.path.join(paths.master, d)):
			if file.endswith('.log'): data_dict['logfile'] = file
			if file.endswith('.out') and not file.startswith('slurm'): data_dict['outfile'] = file
			if file.endswith('.adf.rkf'): data_dict['adfrkf'] = file
			if file.endswith('.ams.rkf'): data_dict['amsrkf'] = file

		summary.summarize_calculation(data_dict)

		data.append(data_dict)

	writer = csv.writer(open(out_path, 'w+', newline=''))
	writer.writerow(['sep=,'])
	writer.writerow(fields)
	for d in data:
		l = [i[1] for i in sorted(d.items(), key=lambda p: fields.index(p[0]))]
		writer.writerow(l)




class DatabaseManager:
	def __init__(self, master_table_path=paths.results_table):
		self.master_table_path = master_table_path
		self.database_reader = csv.reader(open(self.master_table_path, newline=''))
		self.read_data()


	def __enter__(self):
		return self

	def __exit__(self, *args):
		return None

	#database management code
	def read_data(self):
		self.database_reader.__next__()
		self.fieldnames = self.database_reader.__next__()
		self.data = []	
		for d in self.database_reader:
			self.data.append({f:d[i] for i, f in enumerate(self.fieldnames)})

		self.data_dicts = {}
		for i, f in enumerate(self.fieldnames):
			self.data_dicts[f] = [d[f] for d in self.data]


	def get_by_id(self, id):
		id = str(id)
		if self.database_reader is None:
			raise Exception('Database not loaded, use "with"-statement syntax or call __enter__ and __exit__ manually')
		if not id in self.data_dicts['id']:
			print(Exception(f'Calculation with id {id} not found in table {self.master_table_path}'))
			return None

		for d in self.data:
			if d['id'] == id:
				return d 





if __name__ == '__main__':
	create_database_from_calculations()
	with DatabaseManager() as dbm:
		print(dbm.get_by_id(3))

import paths
import csv
import os
import summary, job_results


# def create_database_from_calculations(calc_path=paths.calculations, out_path=paths.results_table):
# 	fields = ['id', 
# 	          'name', 
# 	          'directory', 
# 	          'logfile', 
# 	          'outfile', 
# 	          'adfrkf', 
# 	          'amsrkf']

# 	data = []
# 	dirs = []
# 	for file in os.listdir(calc_path):
# 		p = os.path.join(calc_path, file)
# 		if os.path.isdir(p):
# 			dirs.append(os.path.relpath(p, paths.master))

# 	for d in dirs:
# 		data_dict = {}

# 		data_dict['id'] = os.path.basename(d).split('.')[0]
# 		data_dict['name'] = '.'.join(os.path.basename(d).split('.')[1:])
# 		data_dict['directory'] = d



# 		for file in os.listdir(os.path.join(paths.master, d)):
# 			if file.endswith('.log'): data_dict['logfile'] = file
# 			if file.endswith('.out') and not file.startswith('slurm'): data_dict['outfile'] = file
# 			if file.endswith('.adf.rkf'): data_dict['adfrkf'] = file
# 			if file.endswith('.ams.rkf'): data_dict['amsrkf'] = file

# 		# summary.summarize_calculation(data_dict)

# 		data.append(data_dict)

# 	writer = csv.writer(open(out_path, 'w+', newline=''))
# 	writer.writerow(['sep=,'])
# 	writer.writerow(fields)
# 	for d in data:
# 		l = [i[1] for i in sorted(d.items(), key=lambda p: fields.index(p[0]))]
# 		writer.writerow(l)

def create_database_from_calculations(calc_path=paths.calculations, out_path=paths.results_table):
	fields = ['id',
			  'status',
			  'task',
	          'reaction',
	          'stationary_point', 
	          'R1', 'R2',
	          'radical',
	          'directory',
	          'flags',
	          ]

	
	dirs = []
	for reaction_dir in os.listdir(calc_path):
		p = os.path.join(calc_path, reaction_dir)
		if not os.path.isdir(p): continue
		for calc_dir in os.listdir(p):
			calc_dir = os.path.join(p, calc_dir)
			if os.path.isdir(calc_dir):
				dirs.append(os.path.relpath(calc_dir, paths.master))

	data = []
	for d in dirs:
		data_dict = job_results.get_results(os.path.join(paths.master, d))
		data.append(data_dict)

	#open csv in normal writing mode for first two lines
	writer = csv.writer(open(out_path, 'w+', newline=''))
	writer.writerow(['sep=,'])
	writer.writerow([f.upper() for f in fields])
	#switch to dictwriter for writing data
	writer = csv.DictWriter(open(out_path, 'a', newline=''), fields, extrasaction='ignore')
	for d in data:
		writer.writerow(d)




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

		self.ids = [int(d['id']) for d in self.data]


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

	def get_free_id(self):
		for i in self.ids:
			if not i+1 in self.ids:
				return i+1

	def get_n_free_ids(self, n):
		ids = []
		for i in range(len(self.ids)+n):
			if not i+1 in self.ids and not i+1 in ids:
				ids.append(i+1)
			if len(ids) == n:
				return ids






if __name__ == '__main__':
	create_database_from_calculations()
	# with DatabaseManager(paths.results_table) as dbm:
	# 	print(dbm.get_n_free_ids(10))
	# 	print(dbm.get_by_id(3))

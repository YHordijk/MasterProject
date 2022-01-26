import paths
import csv
import os
import summary, job_results, utility


def create_database_from_calculations(calc_path=paths.calculations, out_path=paths.results_table):
	results = job_results.get_all_results(calc_path)
	fields = results[0]['database_fields']
	writer = csv.writer(open(out_path, 'w+', newline=''))
	writer.writerow(['sep=,'])
	writer.writerow([f.upper() for f in fields])
	for res in results:
		writer.writerow(res['database_data'])


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
		self.results = {d['ID']: job_results.get_results(os.path.join(paths.master,d['CALC_DIRECTORY'])) for d in self.data}
		self.data_dicts = {}
		for i, f in enumerate(self.fieldnames):
			self.data_dicts[f] = [d[f] for d in self.data]
		self.ids = [int(d['ID']) for d in self.data]
		self.hashes = [d['HASH'] for d in self.data]
		all_status = [d['STATUS'] for d in self.data]
		self.Nrunning = all_status.count('R')
		self.Nsucces = all_status.count('S')
		self.Ncancel = all_status.count('C')
		self.Nfail = all_status.count('F')
		self.Nwarn = all_status.count('W')
		self.Ncalcs = len(self.ids)
		self.running_status = self.get_running_status()


	def get_by_id(self, id):
		id = str(id)
		if self.database_reader is None:
			raise Exception('Database not loaded, use "with"-statement syntax or call __enter__ and __exit__ manually')
		if not id in self.data_dicts['ID']:
			print(Exception(f'Calculation with id {id} not found in table {self.master_table_path}'))
			return None

		return self.results[id]


	def get_id_by_hash(self, hash):
		for d in self.data:
			if hash == d['HASH']:
				return d['ID']
		return None


	def get_running_status(self):
		statuses = {}
		running_results = [job_results.get_results(os.path.join(paths.master,d['DIRECTORY'])) for d in self.data if d['STATUS'] == 'R']
		outfiles = [os.path.join(paths.master, r['files']['logfile']) for r in running_results]
		for res, outfile in zip(running_results, outfiles):
			statuses[res['id']] = {}
			#check geom
			if res['task'] in ('GO', 'TSRC'):
				geo_done = False
				freq_progress = 0
				with open(outfile) as f:
					lines = [l.strip() for l in f.readlines()]

					for l in lines:
						if 'constrained gradient max' in l:
							if l.split()[-1] == 'T': 
								geo_done = True
						if '== NUCLEUS' in l:
							freq_progress = int(l.split()[-1])

				# print(res['reaction']+'.'+res['stationary_point'], geo_done, freq_progress)
				statuses[res['id']]['geo_done'] = geo_done
				statuses[res['id']]['freq_progress'] = freq_progress
				statuses[res['id']]['freq_eta'] = 0
		return statuses


def get_hash_collision(hash):
	with DatabaseManager() as dbm:
		return hash in dbm.hashes


def get_n_free_ids(n, write_to_list=True, list_path=paths.id_list):
	old_ids = get_occupied_ids()
	new_ids = []
	for i in range(len(old_ids)+n):
		if not i+1 in old_ids and not i+1 in new_ids:
			new_ids.append(i+1)
		if len(new_ids) == n:
			break
	with open(list_path, 'a') as id_list:
		w = ',' + ','.join([str(i) for i in new_ids])
		id_list.write(w)
	return new_ids


def get_occupied_ids(path=paths.id_list):
	with open(path, 'r') as file:
		return [int(i) for i in file.read().strip().split(',') if i.strip() != '']


def update_id_list(path=paths.id_list):
	create_database_from_calculations()
	with DatabaseManager() as dbm:
		ids = dbm.ids
	with open(path, 'w+') as file:
		file.write(','.join([str(i) for i in ids]))
	return ids


if __name__ == '__main__':
	print(f'Updating database')
	update_id_list()
	print(f'Database updated...')
	with DatabaseManager() as dbm:
		print(f'Found {dbm.Ncalcs} calculations:')
		print(f'\t{dbm.Ncalcs-dbm.Nrunning} finished')
		print(f'\t\t{dbm.Nsucces} success')
		print(f'\t\t{dbm.Nfail} failed')
		print(f'\t\t{dbm.Nwarn} with warnings')
		print(f'\t\t{dbm.Ncancel} canceled')
		print(f'\t{dbm.Nrunning} still running')
		for id, status in dbm.running_status.items():
			res = dbm.get_by_id(id)
			if status['geo_done']:
				natoms = res['natoms']
				s = f'FREQ ({status["freq_progress"]}/{natoms})'
			else:
				s = 'GEO'
			print(f'\t\tCALC: {id :<4} {res["reaction"]}.{res["substituents"]}.{res["stationary_point"]}: {s}')		

	try:
		import excel
	except Exception as e:
		print(e)
		raise
		pass
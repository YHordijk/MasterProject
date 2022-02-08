import paths
import csv
import os
import summary, utility
import job_results2 as job_results



def _result2database(results):
    data = results.data
    d = {}
    d['id'] = data.id
    d['status'] = data.status
    d['task'] = data.task
    d['reaction'] = data.reaction
    d['enantiomer'] = data.enantiomer
    d['R1'] = data.substituents.get('R1', None)
    d['R2'] = data.substituents.get('R2', None)
    d['Rcat'] = data.substituents.get('Rcat', None)
    d['Rch'] = data.substituents.get('Rch', None)
    d['directory'] = os.path.relpath(results.path, paths.results)
    d['stationary_point'] = data.stationary_point
    d['flags'] = ' '.join(list(sorted(data.flags)))
    d['inxyz'] = data.files.get('inxyz', None)
    d['outxyz'] = data.files.get('outxyz', None)
    t = int(data.timings.get('elapsed', 0))
    d['runtime'] = f'{t//3600:0>2}:{(t//60)%60:0>2}:{t%60:0>2}'
    d['hash'] = data.hash
    p = data.progress
    d['progress'] = None
    if p['freq'] > 0:
        d['progress'] = f'FREQ {p["freq"]}'
    elif p['geo'] > 0:
        d['progress'] = f'GEO {p["geo"]}'
    return d





def create_database_table(calc_path=paths.calculations, out_file=paths.results_table):
    fields = ['id',
      'status',
      'task',
      'reaction',
      'enantiomer',
      'R1', 'R2', 'Rcat', 'Rch',
      'stationary_point',
      'directory',
      'flags',
      'inxyz',
      'outxyz',
      'runtime',
      'progress',
      'hash',
      ]

    writer = csv.writer(open(out_file, 'w+', newline=''))
    writer.writerow(['sep=,'])
    writer.writerow([f.upper() for f in fields])

    writer = csv.DictWriter(open(out_file, 'a', newline=''), fields)
    results = job_results.generate_all_results(calc_path)
    if len(results) > 0:
        for res in results:
            data = _result2database(res)
            writer.writerow(data)


class DatabaseManager:
    def __init__(self, master_table_path=paths.results_table):
        self.master_table_path = master_table_path
        if os.path.exists(master_table_path):
            self.database_reader = csv.reader(open(self.master_table_path, newline=''))
            self.read_data()
        else:
            self.data = []

        #get all ids
        self.ids = [d['ID'] for d in self.data]
        self.hashes = [d['HASH'] for d in self.data]
        self.Ncalcs = len(self.data)
        self.statuses = [d['STATUS'] for d in self.data]
        self.progresses = [d['PROGRESS'] for d in self.data]
        self.dirs = [os.path.join(paths.results, d['DIRECTORY']) for d in self.data]


    def __enter__(self):
        return self


    def __exit__(self, *args):
        return None


    #database management code
    def read_data(self):
        lines = [l for l in self.database_reader]
        self.data = []
        if len(lines) > 3:
            fields = lines[1]
            if len(lines) > 2:
                for line in lines[2:]:
                    self.data.append({f:x for f, x in zip(fields, line)})


    def get_by_id(self, id):
        id = str(id)
        i = self.ids.index(id)
        return job_results.Result(self.dirs[i])


    def get_id_by_hash(self, hash):
        for d in self.data:
            if hash == d['HASH']:
                return d['ID']
        return None


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
    create_database_table()
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
        Nrunning = dbm.statuses.count('Running')
        Nsucces = dbm.statuses.count('Success')
        Nfail = dbm.statuses.count('Failed')
        Nwarn = dbm.statuses.count('Warning')
        Nqueued = dbm.statuses.count('Queued')
        print(f'Found {dbm.Ncalcs} calculations:')
        print(f'\t{dbm.Ncalcs-Nrunning-Nqueued} finished')
        print(f'\t\t{Nsucces} success')
        print(f'\t\t{Nwarn} with warnings')
        print(f'\t\t{Nfail} failed')
        print(f'\t{Nqueued} waiting')
        print(f'\t{Nrunning} still running')
        for i in range(dbm.Ncalcs):
            if not dbm.statuses[i] == 'Running': continue
            res = dbm.get_by_id(dbm.ids[i])
            sub = utility.get_sorted_dict_values(res.data.substituents)
            calc_name = f'{res.data.reaction}.{"_".join(sub)}.{res.data.stationary_point}'
            print(f'\t\t{dbm.ids[i] :<4} {calc_name:46} {dbm.progresses[i]:>8}')     

    try:
        import excel
    except Exception as e:
        print('Could not find openpyxl. Did not create nicely formatted excel file ...')
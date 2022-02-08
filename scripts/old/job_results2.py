import scm.plams as plams
import paths, os, struct_generator, utility
import datetime, time, json, shutil
from os.path import join, relpath
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import pickle

np.seterr(all='raise')


# def load_all_results(res_path=paths.results):
#     dirs = []
#     for reaction_dir in os.listdir(res_path):
#       p = join(res_path, reaction_dir)
#       if not os.path.isdir(p): continue
#       for calc_dir in os.listdir(p):
#           calc_dir = join(p, calc_dir)
#           if os.path.isdir(calc_dir):
#               dirs.append(relpath(calc_dir, paths.master))
#     results = []
#     for d in dirs:
#         results.append(Result(join(paths.master, d)))

#     return results


# def get_all_results(res_path=paths.results):
#     dirs = []
#     for reaction_dir in os.listdir(res_path):
#       p = join(res_path, reaction_dir)
#       if not os.path.isdir(p): continue
#       for calc_dir in os.listdir(p):
#           calc_dir = join(p, calc_dir)
#           if os.path.isdir(calc_dir):
#               dirs.append(relpath(calc_dir, paths.master))

#     results = []
#     for d in dirs:
#         if os.path.exists(join(paths.master, d, 'results.dill')):
#             with open(join(paths.master, d, 'results.dill'), 'b') as dill:
#                 results.append(pickle.load(dill))
#         else:
#             results.append(Result(join(paths.master, d)))

    # return results


def generate_all_results(calc_path=paths.calculations):
    dirs = []
    for reaction_dir in os.listdir(calc_path):
      p = join(calc_path, reaction_dir)
      if not os.path.isdir(p): continue
      for calc_dir in os.listdir(p):
          calc_dir = join(p, calc_dir)
          if os.path.isdir(calc_dir):
              dirs.append(relpath(calc_dir, paths.master))

    results = []
    for d in dirs:
        res = generate_result(join(paths.master, d))
        # with open(os.path.join(res.path, 'results.dill'), 'wb+') as dill:
        #     dumps = pickle.dumps(res)
        #     pickle.dump(dumps, dill)
        res.pickle()
        results.append(res)

    return results


def generate_result(path):
    '''
    Function that generates a Result object from a calculation
    path is the directory the calculation was run in
    '''
    def get_files(p):
        files = {}

        metainf = join(os.path.dirname(p), 'meta.info')
        shutil.copy2(metainf, join(p, 'meta.info'))

        for file in os.listdir(path):
            if file.endswith('.out') and not file.startswith('CreateAtoms'):
                files['out'] = file
            if file.endswith('.log'):
                files['log'] = file
            if file.endswith('.in'):
                files['in'] = file
            if file.endswith('.err'):
                files['err'] = file
            if file == 'meta.info':
                files['meta'] = file
            if file == 'adf.rkf':
                files['adfrkf'] = file
            if file == 'ams.rkf':
                files['amsrkf'] = file
            if file == 'output.xyz':
                files['outxyz'] = file

        return files

    def copy_files(files, origin, destination):
        for filet, file in files.items():
            shutil.copy2(join(origin, file), join(destination, file))

    assert os.path.exists(path)

    res_path = join(paths.results, relpath(path, paths.calculations))
    os.makedirs(res_path, exist_ok=True)
    files = get_files(path)
    copy_files(files, path, res_path)
    res = Result(res_path)
    
    return res


class Data:
    pass


class Result:
    '''
    Container that holds results from a calculation and allows for easy access to results
    '''
    def __init__(self, path, regenerate_results=True):
        self.path = path
        #if there is a pickled object availlable unpickle it instead of generating new results object
        if not regenerate_results:
            if self.check_pickle():
                self.data = self.unpickle()
            else:
                self.data = self.get_data()

            #if the job is not yet completed regenerate the object
            if self.data.status in ['Queued', 'Running']:
                self.data = self.get_data()
        else:
            self.data = self.get_data()

        self.pickle()


    def check_pickle(self):
        return os.path.exists(join(self.path, 'results.dill'))


    def pickle(self):
        with open(os.path.join(self.path, 'results.dill'), 'wb+') as dill:
            pickle.dump(self.data, dill, protocol=pickle.HIGHEST_PROTOCOL)


    def unpickle(self):
        with open(os.path.join(self.path, 'results.dill'), 'rb') as dill:
            data = pickle.load(dill)
            return data


    def get_data(self):
        data = Data()
        self._defaults(data)
        self._general(data)
        self._get_files(data)
        self._read_adfrkf(data)
        self._read_log(data)
        self._read_amsrkf(data)
        self._read_meta(data)
        data.hash = struct_generator.get_hash_for_sp(data.reaction, data.substituents, data.stationary_point)
        return data


    def _defaults(self, data):
        data.status = 'Queued'
        data.energies = {}
        data.vibrations = {}
        data.timings = {}
        data.files = {}
        data.Natoms = 0
        data.progress = {'freq': 0, 'geo': 0}
        data.enantiomer = None


    def _general(self, data):
        #first dissect the path itself which is in format
        # {id}.{reaction}.{stationary_point}
        s = os.path.basename(self.path).split('.')
        data.id, data.reaction, data.stationary_point = s[0], s[1], s[-1]
        data.flags = struct_generator.get_flags_for_sp(data.reaction, data.stationary_point)
        for task in ['GO', 'SP', 'TSRC', 'LT']:
            if task in data.flags: 
                data.task = task
        

    def _get_files(self, data):
        for file in os.listdir(self.path):
            if file.endswith('.out') and not file.startswith('CreateAtoms'):
                data.files['out'] = file
            if file.endswith('.log'):
                data.files['log'] = file
            if file.endswith('.in'):
                data.files['in'] = file
            if file.endswith('.err'):
                data.files['err'] = file
            if file == 'adf.rkf':
                data.files['adfrkf'] = file
            if file == 'ams.rkf':
                data.files['amsrkf'] = file
            if file == 'output.xyz':
                data.files['outxyz'] = file
            if file == 'meta.info':
                data.files['meta'] = file


    def _read_adfrkf(self, data):
        if not 'adfrkf' in data.files:
            return
        
        rkf = plams.KFFile(join(self.path, data.files['adfrkf']))
        sections = rkf.sections()

        #read energies
        if 'Energy' in sections:
            data.energies['bond_energy'] = float(rkf.read('Energy', 'Bond Energy'))
            data.energies['eda_elstat'] = float(rkf.read('Energy', 'elstat'))
            data.energies['eda_pauli'] = float(rkf.read('Energy', 'Pauli Total'))
            data.energies['eda_orbit_int'] = float(rkf.read('Energy', 'Orb.Int. Total'))
        if 'Thermodynamics' in sections:
            data.energies['gibbs_energy'] = float(rkf.read('Thermodynamics', 'Gibbs free Energy'))

        #read vibrational info
        if 'vibrations' in sections:
            freqs = rkf.read('Vibrations', 'Frequencies[cm-1]') #freqs is array except for diatoms
            if type(freqs) is float: freqs = [freqs]
            ints = rkf.read('Vibrations', 'Intensities[km/mol]')
            if type(ints) is float: ints = [ints] 
            data.vibrations['frequencies'] = freqs
            data.vibrations['intensities'] = ints
            data.vibrations['freq_nimag'] = len([f for f in freqs if f<0])
            if data.vibrations['freq_nimag'] > 0:
                data.vibrations['freq_imag_displ'] = [float(x) for x in rkf.read('Vibrations', 'NoWeightNormalMode(1)')]
            else:
                data.vibrations['freq_imag_displ'] = []


    def _read_amsrkf(self, data):
        if not 'amsrkf' in data.files:
            return
        try:
            rkf = plams.KFFile(join(self.path, data.files['amsrkf']))

            #get status 
            term = rkf.read('General', 'termination status')
            if 'NORMAL TERMINATION' in term:
                if len(data.warnings) == 0:
                    data.status = 'Success'
                else:
                    data.status = 'Warning'
            elif 'IN PROGRESS' in term:
                data.status = 'Running'
            else:
                data.status = 'Failed'

            data.title = rkf.read('General', 'title')
        except:
            print('error')



    def _read_log(self, data):
        if not 'log' in data.files:
            return

        def get_time(line):
            return datetime.datetime.strptime(' '.join(line.split()[0:2]), '<%b%d-%Y> <%H:%M:%S>')

        with open(join(self.path, data.files['log'])) as log:
            lines = log.readlines()
            start_time = get_time(lines[0])
            end_time = get_time(lines[-1])
            data.timings['elapsed'] = (end_time - start_time).total_seconds()

            #get geometry history:
            GOStep_idx = []
            GOStep_times = []

            #get indices at which *** GOStepx *** is printed
            for i, line in enumerate(lines):
                if '*** GOStep' in line:
                    GOStep_times.append(get_time(line))
                    GOStep_idx.append(i)

            geometries = []
            #starting from the indices add 4 to get to the coordinates and go untill we hit >>>> CORORT
            for i in GOStep_idx:
                g = []
                j = i + 4
                #stop when splits is shorter than 6
                while not 'CORORT' in lines[j+1]:
                    s = lines[j].split()
                    g.append([s[2].split('.')[1], float(s[3]), float(s[4]), float(s[5])])
                    j += 1
                geometries.append(g)
            data.geometries = geometries

            data.Natoms = len(geometries[0])
            data.progress['geo'] = len(geometries)

            utility.write_mol_list(geometries[0], join(self.path, 'input.xyz')) #input molecule
            data.files['inxyz'] = join(self.path, 'input.xyz')
            utility.write_mol_list(geometries[-1], join(self.path, 'output.xyz')) #output molecule
            data.files['outxyz'] = join(self.path, 'output.xyz')

            #get frequency timings
            freq_progress = 0
            freq_times = []
            for line in lines:
                if '=== NUCLEUS:' in line:
                    freq_times.append(get_time(line))
                    freq_progress = int(line.strip().split()[-1])

            data.progress['freq'] = freq_progress

            #Parse timmings
            GOStep_times = [(t - GOStep_times[0]).total_seconds() for t in GOStep_times]
            freq_times = [(t - freq_times[0]).total_seconds() for t in freq_times]
            data.nGOsteps = len(GOStep_times)
            if len(GOStep_times) > 1:
                data.timings['GO']   = [np.mean(np.diff(GOStep_times)), np.std(np.diff(GOStep_times))]
            if len(freq_times) > 1:
                data.timings['freq'] = [np.mean(np.diff(freq_times)), np.std(np.diff(freq_times))]


            #get warnings
            warnings = []
            for line in lines:
                if 'WARNING' in line:
                    w = ' '.join(line.split()[3:])
                    if not w == 'total elapsed time is much higher than the (CPU+system) time.':
                        warnings.append(w)
            data.warnings = warnings


    def _read_meta(self, data):
        if not 'meta' in data.files:
            return

        with open(join(self.path, data.files['meta'])) as meta:
            subs = struct_generator.get_subs_used_for_sp(data.reaction, data.stationary_point)
            data.substituents = {s: None for s in subs}
            lines = meta.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith('R'):
                    sub, n = line.split('=')
                    if sub in data.substituents:
                        data.substituents[sub] = n

            data.flags = data.flags + [f'{R}={n}' for R, n in data.substituents.items()]







if __name__ == '__main__':
    # res = generate_all_results()
    # [print(r.data.hash) for r in res]
    res = generate_result(r"D:\Users\Yuman\Desktop\MasterProject\calculations\achiral_catalyst_H_tBu_ZnCl2\50.achiral_catalyst.H_tBu_ZnCl2.Rsub_cat_complex")
    print(res.data.status)


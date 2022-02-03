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


def get_all_results(res_path=paths.results):
    dirs = []
    for reaction_dir in os.listdir(res_path):
      p = join(res_path, reaction_dir)
      if not os.path.isdir(p): continue
      for calc_dir in os.listdir(p):
          calc_dir = join(p, calc_dir)
          if os.path.isdir(calc_dir):
              dirs.append(relpath(calc_dir, paths.master))

    results = []
    for d in dirs:
        if os.path.exists(join(paths.master, d, 'results.dill')):
            with open(join(paths.master, d, 'results.dill'), 'b') as dill:
                results.append(pickle.load(dill))
        else:
            results.append(Result(join(paths.master, d)))

    return results


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
        with open(os.path.join(res.path, 'results.dill'), 'wb+') as dill:
            dumps = pickle.dumps(res)
            pickle.dump(dumps, dill)
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
            if not os.path.exists(join(destination, file)):
                shutil.copy2(join(origin, file), join(destination, file))

    assert os.path.exists(path)

    res_path = join(paths.results, relpath(path, paths.calculations))
    os.makedirs(res_path, exist_ok=True)
    files = get_files(path)
    copy_files(files, path, res_path)
    res = Result(res_path)
    
    return res


class Result:
    '''
    Container that holds results from a calculation and allows for easy access
    '''
    def __init__(self, path):
        self.path = path
        self.get_data()


    def get_data(self):
        #first dissect the path itself which is in format
        # {id}.{reaction}.{stationary_point}
        self._defaults()
        self._general()
        self._get_files()
        self._read_adfrkf()
        self._read_log()
        self._read_amsrkf()
        self._read_meta()
        self.hash = struct_generator.get_hash_for_sp(self.reaction, self.substituents, self.stationary_point)

    def _defaults(self):
        self.status = 'Queued'
        self.energies = {}
        self.vibrations = {}
        self.timings = {}
        self.files = {}
        self.Natoms = 0
        self.progress = {'freq': 0, 'geo': 0}


    def _general(self):
        s = os.path.basename(self.path).split('.')
        self.id, self.reaction, self.stationary_point = s[0], s[1], s[-1]
        self.flags = struct_generator.get_flags_for_sp(self.reaction, self.stationary_point)
        for task in ['GO', 'SP', 'TSRC', 'LT']:
            if task in self.flags: 
                self.task = task
        

    def _get_files(self):
        for file in os.listdir(self.path):
            if file.endswith('.out') and not file.startswith('CreateAtoms'):
                self.files['out'] = file
            if file.endswith('.log'):
                self.files['log'] = file
            if file.endswith('.in'):
                self.files['in'] = file
            if file.endswith('.err'):
                self.files['err'] = file
            if file == 'adf.rkf':
                self.files['adfrkf'] = file
            if file == 'ams.rkf':
                self.files['amsrkf'] = file
            if file == 'output.xyz':
                self.files['outxyz'] = file
            if file == 'meta.info':
                self.files['meta'] = file


    def _read_adfrkf(self):
        if not 'adfrkf' in self.files:
            return
        
        rkf = plams.KFFile(join(self.path, self.files['adfrkf']))
        sections = rkf.sections()

        #read energies
        if 'Energy' in sections:
            self.energies['bond_energy'] = float(rkf.read('Energy', 'Bond Energy'))
            self.energies['eda_elstat'] = float(rkf.read('Energy', 'elstat'))
            self.energies['eda_pauli'] = float(rkf.read('Energy', 'Pauli Total'))
            self.energies['eda_orbit_int'] = float(rkf.read('Energy', 'Orb.Int. Total'))
        if 'Thermodynamics' in sections:
            self.energies['gibbs_energy'] = float(rkf.read('Thermodynamics', 'Gibbs free Energy'))

        #read vibrational info
        if 'vibrations' in sections:
            freqs = rkf.read('Vibrations', 'Frequencies[cm-1]') #freqs is array except for diatoms
            if type(freqs) is float: freqs = [freqs]
            ints = rkf.read('Vibrations', 'Intensities[km/mol]')
            if type(ints) is float: ints = [ints] 
            self.vibrations['frequencies'] = freqs
            self.vibrations['intensities'] = ints
            self.vibrations['freq_nimag'] = len([f for f in freqs if f<0])
            if self.vibrations['freq_nimag'] > 0:
                self.vibrations['freq_imag_displ'] = [float(x) for x in rkf.read('Vibrations', 'NoWeightNormalMode(1)')]
            else:
                self.vibrations['freq_imag_displ'] = []


    def _read_amsrkf(self):
        if not 'amsrkf' in self.files:
            return
        rkf = plams.KFFile(join(self.path, self.files['amsrkf']))
        sections = rkf.sections()

        #get status 
        term = rkf.read('General', 'termination status')
        if 'NORMAL TERMINATION' in term:
            if len(self.warnings) == 0:
                self.status = 'Success'
            else:
                self.status = 'Warning'
        elif 'IN PROGRESS' in term:
            self.status = 'Running'
        else:
            self.status = 'Failed'

        self.title = rkf.read('General', 'title')



    def _read_log(self):
        if not 'log' in self.files:
            return

        def get_time(line):
            return datetime.datetime.strptime(' '.join(line.split()[0:2]), '<%b%d-%Y> <%H:%M:%S>')

        with open(join(self.path, self.files['log'])) as log:
            lines = log.readlines()
            start_time = get_time(lines[0])
            end_time = get_time(lines[-1])
            self.timings['elapsed'] = (end_time - start_time).total_seconds()

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
                while len(lines[j].split()) == 6:
                    s = lines[j].split()
                    g.append([s[2].split('.')[1], float(s[3]), float(s[4]), float(s[5])])
                    j += 1
                geometries.append(g)
            self.geometries = geometries

            self.Natoms = len(geometries[0])
            self.progress['geo'] = len(geometries)

            utility.write_mol_list(geometries[0], join(self.path, 'input.xyz')) #input molecule
            self.files['inxyz'] = join(self.path, 'input.xyz')
            utility.write_mol_list(geometries[-1], join(self.path, 'output.xyz')) #output molecule
            self.files['outxyz'] = join(self.path, 'output.xyz')

            #get frequency timings
            freq_progress = 0
            freq_times = []
            for line in lines:
                if '=== NUCLEUS:' in line:
                    freq_times.append(get_time(line))
                    freq_progress = int(line.strip().split()[-1])

            self.progress['freq'] = freq_progress

            #Parse timmings
            GOStep_times = [(t - GOStep_times[0]).total_seconds() for t in GOStep_times]
            freq_times = [(t - freq_times[0]).total_seconds() for t in freq_times]
            self.nGOsteps = len(GOStep_times)
            if len(GOStep_times) > 1:
                self.timings['GO']   = [np.mean(np.diff(GOStep_times)), np.std(np.diff(GOStep_times))]
            if len(freq_times) > 1:
                self.timings['freq'] = [np.mean(np.diff(freq_times)), np.std(np.diff(freq_times))]


            #get warnings
            warnings = []
            for line in lines:
                if 'WARNING' in line:
                    w = ' '.join(line.split()[3:])
                    if not w == 'total elapsed time is much higher than the (CPU+system) time.':
                        warnings.append(w)
            self.warnings = warnings


    def _read_meta(self):
        if not 'meta' in self.files:
            return

        with open(join(self.path, self.files['meta'])) as meta:
            subs = struct_generator.get_subs_used_for_sp(self.reaction, self.stationary_point)
            self.substituents = {s: None for s in subs}
            lines = meta.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith('R'):
                    sub, n = line.split('=')
                    if sub in self.substituents:
                        self.substituents[sub] = n

            self.flags = self.flags + [f'{R}={n}' for R, n in self.substituents.items()]









res = generate_all_results()


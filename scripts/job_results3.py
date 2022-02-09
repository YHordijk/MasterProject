import scm.plams as plams
import paths, os, utility
import datetime, time, json, shutil
from os.path import join, relpath
import numpy as np
import matplotlib.pyplot as plt
import json


np.seterr(all='raise')


def get_all_run_dirs(calc_dir=paths.calculations):
    dirs = []
    for system in os.listdir(calc_dir):
        if not os.path.isdir(join(calc_dir, system)): continue
        for dir in os.listdir(join(calc_dir, system)):
            p = join(calc_dir,system,dir)
            dirs.append(p)
    return dirs


def get_all_result_dirs(res_dir=paths.results):
    dirs = []
    for system in os.listdir(res_dir):
        if not os.path.isdir(join(res_dir, system)): continue

        for dir in os.listdir(join(res_dir, system)):
            p = join(res_dir,system,dir)
            dirs.append(p)
    return dirs


def generate_result(calc_path, calc_dir=paths.calculations, res_dir=paths.results):
    assert os.path.exists(calc_path), f'Provided calculation path {calc_path} does not exist'

    res_path = join(res_dir, os.path.relpath(calc_path, calc_dir))
    os.makedirs(res_path, exist_ok=True)

    def get_files():
        files = {}
        job_name = os.path.basename(calc_path)

        #preopt file endings:
        preopt_adf_kf = f'{job_name}_preopt.adf.rkf'
        preopt_ams_kf = f'{job_name}_preopt.ams.rkf'
        preopt_out = f'{job_name}_preopt.out'
        preopt_log = f'{job_name}_preopt.log'
        #opt file endings
        opt_adf_kf = f'{job_name}.adf.rkf'
        opt_ams_kf = f'{job_name}.ams.rkf'
        opt_out = f'{job_name}.out'
        opt_log = f'{job_name}.log'

        files['calc_path'] = calc_path
        files['res_path'] = res_path
        for file in os.listdir(calc_path):
            #input files
            if file == 'input.xyz':
                files['input xyz'] = join(calc_path, file)
            elif file == 'run.info':
                files['run info'] = join(calc_path, file)
            elif file.endswith(job_name):
                files['run script'] = join(calc_path, file)

            #preopt result files
            elif file.endswith(preopt_adf_kf):
                files['preopt adfrkf'] = join(calc_path, file)
            elif file.endswith(preopt_ams_kf):
                files['preopt amsrkf'] = join(calc_path, file)
            elif file.endswith(preopt_out):
                files['preopt out'] = join(calc_path, file)
            elif file.endswith(preopt_log):
                files['preopt log'] = join(calc_path, file)
            #opt results files
            elif file.endswith(opt_adf_kf):
                files['opt adfrkf'] = join(calc_path, file)
            elif file.endswith(opt_ams_kf):
                files['opt amsrkf'] = join(calc_path, file)
            elif file.endswith(opt_out):
                files['opt out'] = join(calc_path, file)
            elif file.endswith(opt_log):
                files['opt log'] = join(calc_path, file)

            if os.path.isdir(join(calc_path, file)) and file == 'ams.results':
                files['ams.results'] = join(calc_path, file)

        if 'ams.results' in files:
            if not 'preopt amsrkf' in files:
                files['preopt amsrkf'] = join(files['ams.results'], 'ams.rkf')
                files['preopt log'] = join(files['ams.results'], 'ams.log')
            elif not 'opt amsrkf' in files:
                files['opt amsrkf'] = join(files['ams.results'], 'ams.rkf')
                files['opt log'] = join(files['ams.results'], 'ams.log')

        return files

    def get_info():
        general = {}
        general['substituents'] = []
        files = data['files']

        if 'run info' in files:
            with open(files['run info']) as info:
                for line in info.readlines():
                    arg, val = line.strip().split('=')

                    if arg == 'reaction':
                        general['reaction'] = val
                    elif arg == 'task':
                        general['task'] = val
                    elif arg == 'phase':
                        general['phase'] = val
                    elif arg == 'radical':
                        general['radical'] = val
                    elif arg == 'stationary_point':
                        general['stationary point'] = val
                    elif arg == 'enantiomer':
                        general['enantiomer'] = val
                    elif arg == 'TSRC':
                        general['TSRC'] = val
                    elif arg.startswith('R'):
                        general['substituents'].append((arg, val))
                    elif arg == 'unique':
                        general['unique'] = val

            if 'unique' in general:
                general['hash'] = utility.hash2(general['unique'])

        if 'opt log' in files:
            logfile = files['opt log']
        elif 'preopt log' in files:
            logfile = files['preopt log']
        else:
            return general


        with open(logfile) as log:
            lines = log.readlines()
            last = lines[-1]
            last = ' '.join(last.split()[0:2])

            now = datetime.datetime.now()
            then = datetime.datetime.strptime(last, '<%b%d-%Y> <%H:%M:%S>')
            general['time silent'] = (now-then).seconds

        return general

    def get_GO_data(preopt=False):
        def write_xyz(file, elements, coords, comment=''):
            with open(file, 'w+') as xyz:
                xyz.write(str(len(elements)))
                xyz.write(f'\n{comment}\n')
                for e, c in zip(elements, coords):
                    xyz.write(f'{e:2}\t{c[0]:11.9f}\t{c[1]:11.9f}\t{c[2]:11.9f}\n')

        def get_runtime(file):
            with open(file) as log:
                lines = log.readlines()
                first = datetime.datetime.strptime(' '.join(lines[0].split()[0:2]),  '<%b%d-%Y> <%H:%M:%S>')
                last  = datetime.datetime.strptime(' '.join(lines[-1].split()[0:2]), '<%b%d-%Y> <%H:%M:%S>')
            return (last-first).seconds


        opt = {}
        files = data['files'] 
        amsrkf = 'opt amsrkf'
        input_xyz = join(res_path, 'input.xyz')
        output_xyz = join(res_path, 'output.xyz')
        log = 'opt log'
        # xyz_comment =
        if preopt:
            amsrkf = 'preopt amsrkf'
            input_xyz = join(res_path, 'input_preopt.xyz')
            output_xyz = join(res_path, 'output_preopt.xyz')
            log = 'preopt log'

        if amsrkf in files:
            ams = plams.KFFile(files[amsrkf])
            opt['runtime'] = get_runtime(files[log])

            opt['natoms'] = ams.read('InputMolecule', 'nAtoms')
            opt['elements'] = ams.read('InputMolecule', 'AtomSymbols').split()
            c = ams.read('InputMolecule', 'Coords')
            opt['input coords'] = [c[i:i+3] for i in range(0, len(c), 3)]
            write_xyz(input_xyz, opt['elements'], opt['input coords'])
            opt['input xyz'] = input_xyz


            try:
                opt['steps'] = ams.read('History', 'nEntries')
                c = ams.read('History', f'Coords({opt["steps"]})')
                opt['output coords'] = [c[i:i+3] for i in range(0, len(c), 3)]
                opt['output energy'] = ams.read('History', f'Energy({opt["steps"]})')
                write_xyz(output_xyz, opt['elements'], opt['output coords'])
                opt['output xyz'] = output_xyz
            except: 
                opt['steps'] = 0
            
            opt['status'] = ams.read('General', 'termination status')

        return opt

    def get_freq_data():
        freq = {}
        files = data['files']
        adfrkf = 'adfrkf'

        if adfrkf in files:
            adf = plams.KFFile(files[adfrkf])
            freq['nmodes'] = adf.read('Vibrations', 'nNormalModes')
            freq['frequencies'] = adf.read('Vibrations', 'Frequencies[cm-1]')
            freq['intensities'] = adf.read('Vibrations', 'Intensities[km/mol]')
            freq['nimag'] = len(f for f in freq['frequencies'] if f < 0)
            if freq['nimag'] > 0:
                freq['imag mode'] = adf.read('Vibrations', 'NoWeightNormalMode(1)')
            freq['ZPE'] = adf.read('Vibrations', 'ZeroPointEnergy')

        return freq

    def get_thermo_data():
        thermo = {}
        files = data['files']
        adfrkf = 'opt adfrkf'

        if adfrkf in files:
            adf = plams.KFFile(files[adfrkf])
            thermo['gibbs'] = adf.read('Thermodynamics', 'Gibbs free Energy')
            thermo['enthalpy'] = adf.read('Thermodynamics', 'Enthalpy')

        return thermo

    def get_misc_data():
        misc = {}
        files = data['files']
        adfrkf = 'opt adfrkf'

        # if adfrkf in files:
        #     adf = plams.KFFile(files[adfrkf])

        return misc


    data = {}
    data['files']   = get_files()
    data['info']    = get_info()
    data['pre opt'] = get_GO_data(True)
    data['opt']     = get_GO_data()
    data['freq']    = get_freq_data()
    data['thermo']  = get_thermo_data()
    data['misc']    = get_misc_data()

    with open(join(res_path, 'results.json'), 'w+') as res:
        res.write(json.dumps(data, indent=2))

    return Result(join(res_path))


def get_all_results(calc_dir=paths.calculations, res_dir=paths.results, regenerate_all=True):
    calc_paths = get_all_run_dirs(calc_dir) #get all calc directories
    res = []
    for calc_path in calc_paths:

        if regenerate_all:
            r = generate_result(calc_path, calc_dir=calc_dir, res_dir=res_dir)
        else:
            #get res_path
            res_path = join(res_dir, os.path.relpath(calc_path, calc_dir))
            #check if it exists first
            if not os.path.exists(res_path):
                #if not, generate the results
                r = generate_result(calc_path, calc_dir=calc_dir, res_dir=res_dir)
            else:
                r = Result(res_path)
                #check status
                if r.status in ['Running', 'Queued']:
                    if os.path.exists(calc_path):
                        r = generate_result(calc_path, calc_dir=calc_dir, res_dir=res_dir)

        res.append(r)
    return res


class Result:
    def __init__(self, path):
        self.path = path
        self.data = self.read_data()
        self._set_status()

    def read_data(self):
        with open(join(self.path, 'results.json'), 'r') as res:
            return json.loads(res.read())

    # some handy functions
    @property
    def gibbs(self):
        return self.data['thermo'].get('gibbs', None)

    @property
    def energy(self):
        if 'opt' in self.data:
            return self.data['opt'].get('output energy')
        return self
    

    @property
    def hash(self):
        return self.data['info'].get('hash', None)

    @property
    def reaction(self):
        return self.data['info'].get('reaction', None)

    @property
    def substituents(self):
        return self.data['info'].get('substituents', [])

    @property
    def radical(self):
        return self.data['info'].get('radical', False)

    @property
    def task(self):
        return self.data['info'].get('task', None)

    @property
    def enantiomer(self):
        return self.data['info'].get('enantiomer', None)

    def get_substituent(self, R):
        s = self.substituents
        for sub in s:
            if sub[0] == R:
                return sub[1]

    @property
    def stationary_point(self):
        return self.data['info'].get('stationary point', None)

    @property
    def runtime(self):
        return self.data['pre opt'].get('runtime', 0) + self.data['opt'].get('runtime', 0)

    def _set_status(self):
        preopt_status = self.data['pre opt'].get('status', None)
        opt_status    = self.data['opt'].get('status', None)

        self.status = 'Queued'
        self.step = ''

        if preopt_status in ['IN PROGRESS', 'NORMAL TERMINATION']:
            self.status = 'Running'
            self.step = 'PreOpt'

        if opt_status is None:
            return 

        self.step = 'GO'
        if opt_status == 'IN PROGRESS':
            self.status = 'Running'
        elif opt_status == 'NORMAL TERMINATION':
            self.status = 'Success'
            self.step = ''
        elif 'WARNING' in opt_status:
            self.status = 'Warning'
            self.step = ''
        else:
            self.status = 'Failed'
            self.step = ''

        if self.status in ['Success', 'Warning']:
            self.step = ''

        #if the job is still running check the time since the last message
        #if the job has not sent a message for 2 hours we consider it cancelled
        if self.status == 'Running':
            time = self.data['info']['time silent']
            if time/3600 >= 2:
                self.status = 'Canceled'

        with open(self.data['files']['opt log']) as log:
            if '=== NUCLEUS' in log.read():
                self.step = 'FREQ'


def summarize_calculations(res, tabs=0):
    statuses = [r.status for r in res]
    njobs = len(statuses)
    nSuccess = statuses.count('Success')
    nFailed = statuses.count('Failed')
    nWarn = statuses.count('Warning')
    nQueued = statuses.count('Queued')
    nRunning = statuses.count('Running')
    nCanceled = statuses.count('Canceled')

    print('\t'*tabs + f'Found {njobs} jobs!')
    print('\t'*tabs + f'\tSuccessful     = {nSuccess} ({nWarn} warnings)')
    print('\t'*tabs + f'\tFailed         = {nFailed}')
    print('\t'*tabs + f'\tCanceled       = {nCanceled}')
    print('\t'*tabs + f'\tQueued         = {nQueued}')
    print('\t'*tabs + f'\tRunning        = {nRunning}')

    reaction_len = max(len(r.reaction) for r in res)
    point_len = max(len(r.stationary_point) for r in res)
    sub_names = list(sorted(set( s[0] for r in res for s in r.substituents )))
    subs = [[r.get_substituent(R) for R in sub_names] for r in res]
    subs = [[(s, '')[s is None] for s in sub] for sub in subs]
    sub_len = max([max(len(s) for s in sub) for sub in subs])
    step_len = max(len(r.step) for r in res)

    print()
    header = '\t'*tabs +  f'\t{"Reaction".ljust(reaction_len)} | {"Point".ljust(point_len)} | '
    header += (' | ').join([s.ljust(sub_len) for s in sub_names])
    header += f' | {"Step".ljust(step_len)}'
    print(header)
    print('\t'*(tabs+1) + '-'*(len(header)))
    for r in res:
        if r.status != 'Running': continue
        l = '\t'*tabs + f'\t{r.reaction:{reaction_len}} | {r.stationary_point:{point_len}} | '
        l += (' | ').join([r.get_substituent(R).ljust(sub_len) for R in sub_names])
        l += f' | {r.step.ljust(step_len)}'
        print(l)
    print()


if __name__ == '__main__':
    res = get_all_results(calc_dir=join(paths.master, 'calculations_test'), regenerate_all=True)
    summarize_calculations(res)

 
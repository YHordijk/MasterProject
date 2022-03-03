import scm.plams as plams
import paths, os, utility, struct_generator2
os.environ['PYGAME_HIDE_SUPPORT_PROMPT'] = "hide"
import datetime, time, json, shutil
from os.path import join, relpath
import geometry
import numpy as np
import matplotlib.pyplot as plt
try:
    import molviewer2.molecule as molecule
except:
    pass
import json


b2a = utility.bohr2angstrom

np.seterr(all='raise')


def get_all_run_dirs(calc_dir=paths.calculations):
    dirs = []
    for system in os.listdir(calc_dir):
        if not os.path.isdir(join(calc_dir, system)): continue
        for dir in os.listdir(join(calc_dir, system)):
            p = join(calc_dir,system,dir)
            dirs.append(p)
    return dirs

all_run_dirs = get_all_run_dirs()


def get_all_result_dirs(res_dir=paths.results):
    dirs = []
    for system in os.listdir(res_dir):
        if not os.path.isdir(join(res_dir, system)): continue

        for dir in os.listdir(join(res_dir, system)):
            p = join(res_dir,system,dir)
            dirs.append(p)
    return dirs

all_result_dirs = get_all_result_dirs()


def generate_result(calc_path, calc_dir=paths.calculations, res_dir=paths.results):
    assert os.path.exists(calc_path), f'Provided calculation path {calc_path} does not exist'

    res_path = join(res_dir, os.path.relpath(calc_path, calc_dir))
    os.makedirs(res_path, exist_ok=True)

    def get_files():
        files = {}
        job_name = os.path.basename(calc_path).split('.')[0]

        #opt file endings
        GO_adf_kf = f'{job_name}.adf.rkf'
        GO_ams_kf = f'{job_name}.ams.rkf'
        GO_out = f'{job_name}.out'
        GO_log = f'{job_name}.log'

        rel_calc_path = os.path.relpath(calc_path, calc_dir)
        files['calc_path'] = rel_calc_path
        files['res_path'] = os.path.relpath(res_path, res_dir)
        for file in os.listdir(calc_path):
            #input files
            if file == 'input.xyz':
                files['input xyz'] = file
            elif file == 'run.info':
                files['run info'] = file
            elif file.endswith(job_name):
                files['GO run script'] = file

            #fragment files 
            elif file == 'Fragment_1.log':
                files['frag1 log'] = file
            elif file == 'Fragment_1.out':
                files['frag1 out'] = file
            elif file == 'Fragment_1.rkf':
                files['frag1 rkf'] = file
            elif file == 'Fragment_1.adf.rkf':
                files['frag1 adfrkf'] = file

            elif file == 'Fragment_2.log':
                files['frag2 log'] = file
            elif file == 'Fragment_2.out':
                files['frag2 out'] = file
            elif file == 'Fragment_2.rkf':
                files['frag2 rkf'] = file
            elif file == 'Fragment_2.adf.rkf':
                files['frag2 adfrkf'] = file

            #EDA files
            elif file == 'EDA.log':
                files['EDA log'] = file
            elif file == 'EDA.ams.rkf':
                files['EDA amsrkf'] = file
            elif file == 'EDA.adf.rkf':
                files['EDA adfrkf'] = file
            elif file == 'frag.run.out':
                files['EDA out'] = file
            elif file == 'frag.run':
                files['EDA run script'] = file

            #GO results files
            elif file.endswith(GO_adf_kf):
                files['GO adfrkf'] = file
            elif file.endswith(GO_ams_kf):
                files['GO amsrkf'] = file
            elif file.endswith(GO_out):
                files['GO out'] = file
            elif file.endswith(GO_log):
                files['GO log'] = file

        return files


    def get_info():
        general = {}
        general['substituents'] = []
        files = data['files']

        if 'run info' in files:
            fileinfo = join(calc_dir, calc_path, files['run info'])
            with open(fileinfo) as info:
                for line in info.readlines():
                    try:
                        arg, val = line.strip().split('=')
                    except:
                        print(files['run info'])

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
                    elif arg == 'functional':
                        general['functional'] = val
                    elif arg == 'basis':
                        general['basis'] = val
                    elif arg == 'numerical_quality':
                        general['numerical_quality'] = val

            if 'unique' in general:
                general['hash'] = utility.hash2(general['unique'])

            try:
                subdict = {x[0]:x[1] for x in general['substituents']}
                mol = struct_generator2.generate_stationary_points(template=general['reaction'], substituents=subdict)[general['stationary point']]
                if hasattr(mol, 'active_atom_idx'):
                    general['active atom index'] = mol.active_atom_idx
            except:
                pass

        return general


    def get_GO_data():
        def write_xyz(file, elements, coords, comment=''):
            with open(file, 'w+') as xyz:
                xyz.write(str(len(elements)))
                xyz.write(f'\n{comment}\n')
                for e, c in zip(elements, coords):
                    xyz.write(f'{e:2}\t{b2a(c[0]):11.9f}\t{b2a(c[1]):11.9f}\t{b2a(c[2]):11.9f}\n')

        GO = {}
        files = data['files'] 
        amsrkf = join(calc_dir, calc_path, files['GO amsrkf'])
        try:
            adfrkf = join(calc_dir, calc_path, files['GO adfrkf'])
        except:
            print(files)
        input_xyz = join(res_path, 'input.xyz')
        output_xyz = join(res_path, 'output.xyz')
        aligned_xyz = join(res_path, 'aligned.xyz')
        log = join(calc_dir, calc_path, files['GO log'])
        GO['runtime'] = 0
        GO['done'] = False

        if 'GO amsrkf' in files:
            ams = plams.KFFile(amsrkf)
            try:
                GO['natoms'] = ams.read('InputMolecule', 'nAtoms')
                GO['elements'] = ams.read('InputMolecule', 'AtomSymbols').split()
                c = ams.read('InputMolecule', 'Coords')
                GO['input coords'] = [c[i:i+3] for i in range(0, len(c), 3)]
                write_xyz(input_xyz, GO['elements'], GO['input coords'])
                GO['input xyz'] = os.path.relpath(input_xyz, res_path)
            except:
                raise
            
            try:
                GO['steps'] = ams.read('History', 'nEntries')
                c = ams.read('History', f'Coords({GO["steps"]})')
                GO['output coords'] = [c[i:i+3] for i in range(0, len(c), 3)]
                GO['output energy'] = ams.read('History', f'Energy({GO["steps"]})')
                write_xyz(output_xyz, GO['elements'], GO['output coords'])
                GO['output xyz'] = os.path.relpath(output_xyz, res_path)
            except: 
                GO['steps'] = 0

            try:
                xyz = np.array(GO['output coords'])
                el = GO['elements']
                m = plams.Molecule()
                for e, x in zip(el, xyz):
                    m.add_atom(plams.Atom(symbol=e, coords=x))
                subdict = {x[0]:x[1] for x in data['info']['substituents']}
                tmol = struct_generator2.get_mol(data['info']['reaction'], subdict, data['info']['stationary point'])
                geometry.align_molecule_to_plane(m, tmol.plane_idx)
                geometry.rotate_molecule_in_plane(m, tmol.align_idx)
                geometry.center_molecule(m, tmol.center_idx)
                GO['aligned coords'] = [a.coords for a in m.atoms]
                print(GO['aligned coords'])
                write_xyz(aligned_xyz, GO['elements'], GO['aligned coords'])
            except:
                raise

            GO['status'] = ams.read('General', 'termination status')

        if 'GO adfrkf' in files:
            adf = plams.KFFile(adfrkf)
            GO['bond energy'] = adf.read('Energy', 'Bond Energy')
            GO['Mulliken charge'] = adf.read('Properties', 'AtomCharge Mulliken')
            GO['electron dens at nuclei'] = adf.read('Properties', 'Electron Density at Nuclei')
            GO['Hirshfeld charge'] = adf.read('Properties', 'FragmentCharge Hirshfeld')
            GO['Voronoi charge'] = adf.read('Properties', 'AtomCharge_SCF Voronoi')

        if 'GO log' in files:
            with open(log) as log:
                lines = log.readlines()
                first = datetime.datetime.strptime(' '.join(lines[0].split()[0:2]),  '<%b%d-%Y> <%H:%M:%S>')
                last  = datetime.datetime.strptime(' '.join(lines[-1].split()[0:2]), '<%b%d-%Y> <%H:%M:%S>')
                GO['runtime'] = int((last-first).total_seconds())

                if any('NORMAL TERMINATION' in line for line in lines):
                    GO['done'] = True

                warning_lines = [line for line in lines if 'WARNING:' in line]
                warning_lines = [l for l in warning_lines if 'total elapsed time is much higher than the (CPU+system) time' not in l]
                GO['warnings'] = warning_lines

                error_lines = [line for line in lines if 'ERROR:' in line]
                GO['errors'] = error_lines

        return GO


    def get_freq_data():
        freq = {}
        files = data['files']
        adfrkf = join(calc_dir, calc_path, files['GO adfrkf'])

        if 'GO adfrkf' in files:
            adf = plams.KFFile(adfrkf)
            try:
                freq['nmodes'] = adf.read('Vibrations', 'nNormalModes')
                freq['frequencies'] = adf.read('Vibrations', 'Frequencies[cm-1]')
                freq['intensities'] = adf.read('Vibrations', 'Intensities[km/mol]')
                if type(freq['frequencies']) is float: freq['frequencies'] = [freq['frequencies']]
                if type(freq['intensities']) is float: freq['intensities'] = [freq['intensities']]
                freq['nimag'] = len([f for f in freq['frequencies'] if f < 0])
                if freq['nimag'] > 0:
                    freq['imag mode'] = adf.read('Vibrations', 'NoWeightNormalMode(1)')
                freq['ZPE'] = adf.read('Vibrations', 'ZeroPointEnergy')
            except:
                pass

        return freq

    def get_thermo_data():
        thermo = {}
        files = data['files']
        adfrkf = join(calc_dir, calc_path, files['GO adfrkf'])
        if 'GO adfrkf' in files:
            adf = plams.KFFile(adfrkf)
            if 'Thermodynamics' in adf.sections():
                thermo['gibbs'] = adf.read('Thermodynamics', 'Gibbs free Energy')
                thermo['enthalpy'] = adf.read('Thermodynamics', 'Enthalpy')

        return thermo


    def get_EDA_data():
        files = data['files']
        EDA = {}
        if not 'EDA adfrkf' in files:
            return
        try:
            p = join(calc_dir, calc_path, files['EDA adfrkf'])
            EDAadfrkf = plams.KFFile(p)
            EDA['pauli'] = EDAadfrkf.read('Energy', 'Pauli Total')
            EDA['bonding'] = EDAadfrkf.read('Energy', 'Bond Energy')
            EDA['elstat'] = EDAadfrkf.read('Energy', 'Electrostatic Interaction')
            EDA['orbital interaction'] = EDAadfrkf.read('Energy', 'Orb.Int. Total')
        except:
            print(files['calc_path'])
            raise

        return EDA



    data = {}
    data['files']   = get_files()
    data['info']    = get_info()
    if data['info']['reaction'] == 'achiral_catalyst':
        if data['info']['stationary point'] not in ['sub_cat_complex', 'TS', 'P1_cat_complex']: return
    data['GO']      = get_GO_data()
    data['freq']    = get_freq_data()
    data['thermo']  = get_thermo_data()
    data['EDA']     = get_EDA_data()
    with open(join(res_path, 'results.json'), 'w+') as res:
        res.write(json.dumps(data, indent=2))

    return Result(join(res_path))


def get_all_results(calc_dir=paths.calculations, res_dir=paths.results, regenerate_all=False):
    def check_regenerate(calc_path, res_path):
        check = False
        reason = 'None'
        if regenerate_all: 
            check = True
            reason = 'regenerate_all'
        #if the path does not exist then in all cases do not regenerate
        if not os.path.exists(calc_path):
            check = False
            reason = 'calc_path not found'

        #if it does exist check the status
        try:
            r = Result(res_path, calc_path=calc_path)
            if r.status in ['Running', 'Queued', 'Canceled']:
                reason = 'calc still running'
                check = True
        except:
            pass
        return check

    calc_path_dict = {}
    for calc_path in get_all_run_dirs(calc_dir):
        res_path = join(res_dir, os.path.relpath(calc_path, calc_dir))
        calc_path_dict[res_path] = calc_path
        if check_regenerate(calc_path, res_path):
            generate_result(calc_path, calc_dir=calc_dir, res_dir=res_dir)
    
    res = []
    for res_path in get_all_result_dirs(res_dir):
        print(res_path)
        r = Result(res_path, calc_path=calc_path_dict.get(res_path, None))
        res.append(r)

    return res



class Result:
    def __init__(self, path, calc_path=None):
        self.path = path
        self._calc_path = calc_path
        self.data = self.read_data()
        self._set_status()
        # try:
        #     self.write_aligned_xyz()
        # except:
        #     pass

    def read_data(self):
        with open(join(self.path, 'results.json'), 'r') as res:
            return json.loads(res.read())

    # some handy functions
    @property
    def gibbs(self):
        return self.data['thermo'].get('gibbs', None)

    @property
    def energy(self):
        if 'bond energy' in self.data['GO']:
            return self.data['GO'].get('bond energy')

        return self.data['GO'].get('output energy')

    @property
    def phase(self):
        return self.data['info'].get('phase', 'vacuum')
    
    @property
    def hash(self):
        return self.data['info'].get('hash', None)

    @property
    def functional(self):
        return self.data['info'].get('functional', None)

    @property
    def basis(self):
        return self.data['info'].get('basis', None)

    @property
    def numerical_quality(self):
        return self.data['info'].get('numerical_quality', None)
    
    @property
    def reaction(self):
        return self.data['info'].get('reaction', None)

    @property
    def substituents(self):
        s = self.data['info'].get('substituents', [])
        return {x[0]:x[1] for x in s}

    @property
    def sorted_substituents(self):
        return utility.get_sorted_dict_values(self.substituents)

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
        return ''

    @property
    def stationary_point(self):
        return self.data['info'].get('stationary point', None)

    @property
    def runtime(self):
        return self.data['GO']['runtime']

    @property
    def formatted_runtime(self):
        t = self.runtime
        return f'{t//3600:0>2}:{(t//60)%60:0>2}:{t%60:0>2}'

    @property
    def calc_path(self):
        return self._calc_path

    @property
    def natoms(self):
        return self.data['GO'].get('natoms')

    @property
    def get_imaginary_mode(self):
        return self.data['freq'].get('imag mode', None)

    def get_mol(self):
        xyz = join(self.path, self.data['GO'].get('output xyz'))
        mol = plams.Molecule(xyz)
        mol.name = self.stationary_point
        mol.reaction = self.reaction
        mol.normalmode = self.get_imaginary_mode
        mol.substituents = self.substituents
        mol.template_mol = self.get_template_mol()
        return mol

    def get_template_mol(self):
        mol = struct_generator2.get_mol(self.reaction, self.substituents, self.stationary_point)
        return mol

    def get_aligned_mol(self):
        xyz = join(self.path, self.data['GO'].get('aligned xyz'))
        mol = plams.Molecule(xyz)
        mol.name = self.stationary_point
        mol.reaction = self.reaction
        mol.normalmode = self.get_imaginary_mode
        mol.substituents = self.substituents
        mol.template_mol = self.get_template_mol()
        return mol

    def _set_status(self):
        opt_status    = self.data['GO'].get('status', None)

        self.status = 'Queued'

        if opt_status is None: 
            return 
        if opt_status == 'IN PROGRESS':
            self.status = 'Running'
        elif opt_status == 'NORMAL TERMINATION':
            self.status = 'Success'
        elif 'warnings' in opt_status:
            self.status = 'Warning'
        elif 'errors' in opt_status:
            self.status = 'Error'
        else:
            self.status = 'Failed'

        #if the job is still running check the time since the last message
        #if the job has not sent a message for 2 hours we consider it cancelled
        if self.status == 'Running':
            time = self.data['info']['time silent']
            if time/3600 >= 2:
                self.status = 'Canceled'

        if self.status == 'Warning':
            warnings = self.data['GO']['warnings']
            if len(warnings) == 0:
                self.status = 'Success'
            else:
                if all('CPKS failed to converge' in warning for warning in warnings):
                    if self.natoms == 1:
                        self.status = 'Success'


    
def summarize_calculations(res, tabs=0):
    statuses = [r.status for r in res]
    njobs = len(statuses)
    nSuccess = statuses.count('Success')
    nFailed = statuses.count('Failed')
    nWarn = statuses.count('Warning')
    nError = statuses.count('Error')
    nQueued = statuses.count('Queued')
    nRunning = statuses.count('Running')
    nCanceled = statuses.count('Canceled')

    print('\t'*tabs + f'Found {njobs} jobs!')
    print('\t'*tabs + f'\tSuccessful     = {nSuccess+nWarn+nError} ({nWarn} warnings, {nError} errors)')
    print('\t'*tabs + f'\tFailed         = {nFailed}')
    print('\t'*tabs + f'\tCanceled       = {nCanceled}')
    print('\t'*tabs + f'\tQueued         = {nQueued}')
    print('\t'*tabs + f'\tRunning        = {nRunning}')

    reaction_len = max(max(len(r.reaction) for r in res), len('Reaction'))
    point_len = max(max(len(r.stationary_point) for r in res), len('Point'))
    sub_names = list(sorted(set( s for r in res for s in r.substituents if r.status == 'Running' )))
    subs = [[r.substituents.get(R, '') for r in res if r.status == 'Running'] for R in sub_names]
    sub_lens = [max(max(len(s) for s in sub), len(R)) for sub, R in zip(subs, sub_names)]
    step_len = max(max(len(r.step) for r in res), len('Step'))

    header = '\t'*tabs +  f'\t{"Reaction".ljust(reaction_len)} | {"Point".ljust(point_len)} | '
    header += (' | ').join([R.ljust(l) for R, l in zip(sub_names, sub_lens)])
    header += f' | {"Step".ljust(step_len)} | Progress | Runtime  '
    print()
    print(header)
    print('\t'*(tabs+1) + '-'*(len(header)))
    for r in res:
        if r.status != 'Running': continue
        if r.step == 'FREQ':
            progress = f"{r.data['opt'].get('freq progress','???')}/{r.data['opt'].get('natoms', '???')}"
        else:
            progress = ''
        l = '\t'*tabs + f'\t{r.reaction:{reaction_len}} | {r.stationary_point:{point_len}} | '
        l += (' | ').join([r.substituents.get(R,'').ljust(l) for R, l in zip(sub_names, sub_lens)])
        l += f' | {r.step.ljust(step_len)} | {progress.rjust(8)} | {r.formatted_runtime}'
        print(l)
    print()

def print_failed(res, tabs=0):
    failed = [r for r in res if r.status == 'Failed']
    print(f'Found {len(failed)} failed jobs')
    for r in failed:
        print(r.calc_path)

def print_canceled(res, tabs=0):
    canceled = [r for r in res if r.status == 'Canceled']
    print(f'Found {len(canceled)} canceled jobs')
    for r in canceled:
        print(r.calc_path)


def get_result(template, substituents, stationary_point):
    results = all_results
    results = [r for r in results if r.reaction == template]
    results = [r for r in results if all(sub == substituents[R] for R, sub in r.substituents.items())]
    results = [r for r in results if all((r.functional == functional, r.basis == basis, r.numerical_quality == numerical_quality))]
    results = [r for r in results if r.phase == phase]
    results = [r for r in results if r.stationary_point == stationary_point]
    if len(results) == 0:
        return
    else:
        return results[0]


all_results = get_all_results()


if __name__ == '__main__':
    calc_dir = join(paths.master, 'calculations_test')
    calc_dir = paths.calculations
    res = get_all_results(calc_dir=calc_dir, regenerate_all=True)
    # summarize_calculations(res)
    print_failed(res)
    print_canceled(res)

 
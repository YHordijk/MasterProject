import os, paths, datetime, json
import scm.plams as plams
join = os.path.join


match_dict = lambda a, b: all(a[k] == b.get(k, None) for k in a.keys())
match_dict2 = lambda a, b: match_dict(a,b) or match_dict(b,a)


class Result:
    def __init__(self, data):
        self.data = data

    def get_main_molecule(self):
        return plams.Molecule(self._path('output_xyz'))

    def get_solute_path(self):
        return self._path('coskf')

    def _path(self, p, calc_dir=paths.calculations2):
        return join(calc_dir, self.data['files'][p])

    def match(self, **settings):
        base = all(self[key] == val for key, val in settings.items() if key != 'substituents')

        if 'substituents' in settings:
            return base and match_dict2(settings['substituents'], self['substituents'])
        return base

    def __eq__(self, other):
        keys = ['reaction', 
                'functional', 
                'basis', 
                'numerical_quality', 
                'frozen_core', 
                'phase', 
                'stationary_point']

        base = all(self[key] == other[key] for key in keys)
        return base and match_dict2(other['substituents'], self['substituents'])

    def __getitem__(self, key):
        if key in ['reaction', 'functional', 'basis', 'numerical_quality', 'frozen_core', 'phase', 'task', 'stationary_point']:
            return self.data['general'][key]
            
        elif key == 'substituents':
            return {R: self.data['general'][R] for R in self.data['general'].keys() if R.startswith('R')}

        elif key == 'gibbs':
            if self['task'] == 'COSMORS':
                return self.data['COSMORS']['gibbs']
            elif self['task'] == 'FREQ':
                return self.data['FREQ']['gibbs']
            return



def load_result(run_dir):
    with open(join(run_dir, 'job.results'), 'r') as res:
        data = json.loads(res.read())
    return Result(data)


def generate_result(run_dir, out_file=None, calc_dir=paths.calculations2):
    path = lambda p: join(calc_dir, data['files'][p])

    if out_file is None:
        out_file = join(run_dir, 'job.results')

    def get_files():
        files = {}

        rel_calc_path = os.path.relpath(run_dir, calc_dir)
        files['calc_path']     = rel_calc_path
        files['adfrkf']        = join(rel_calc_path, 'ADF','adf.rkf')
        files['amsrkf']        = join(rel_calc_path, 'ADF','ams.rkf')
        files['out']           = join(rel_calc_path, 'ADF','ADF.out')
        files['log']           = join(rel_calc_path, 'ADF','ams.log')
        files['err']           = join(rel_calc_path, 'ADF','ADF.err')
        files['coskf']         = join(rel_calc_path, 'ADF','solute.coskf')
        files['COSMORS_out']   = join(rel_calc_path, 'COSMORS', 'COSMORS.out')
        files['COSMORS_in']    = join(rel_calc_path, 'COSMORS', 'COSMORS.in')
        files['COSMORS_crskf'] = join(rel_calc_path, 'COSMORS', 'COSMORS.crskf')
        files['in']            = join(rel_calc_path, 'job.input')
        files['hash']          = join(rel_calc_path, 'job.hash')
        files['info']          = join(rel_calc_path, 'job.info')
        files['input_xyz']     = join(rel_calc_path, 'input.xyz')
        files['output_xyz']    = join(rel_calc_path, 'ADF', 'output.xyz')

        #check the files
        missing_files = []
        for name, file in files.items():
            if not os.path.exists(join(calc_dir, file)):
                # print('Could not find file', name, join(calc_dir, file))
                missing_files.append(name)
        [files.pop(name) for name in missing_files]

        return files

    def read_info():
        files = data['files']
        general = {}

        with open(path('info')) as info:
            for line in info.readlines():
                arg, val = line.strip().split('=')
                general[arg] = val

        if 'plane_idx' in general:
            general['plane_idx'] = [int(i) for i in general['plane_idx'].split('_')]
        if 'align_idx' in general:
            general['align_idx'] = [int(i) for i in general['align_idx'].split('_')]
        if 'TSRC_idx' in general:
            if general['TSRC_idx'] == '':
                general['TSRC_idx'] = None
            else:
                general['TSRC_idx'] = [int(i) for i in general['TSRC_idx'].split('_')]
        if 'center_idx' in general:
            general['center_idx'] = int(general['center_idx'])

        return general


    def read_log():
        first = datetime.datetime.strptime(' '.join(loglines[0].split()[0:2]),  '<%b%d-%Y> <%H:%M:%S>')
        last  = datetime.datetime.strptime(' '.join(loglines[-1].split()[0:2]), '<%b%d-%Y> <%H:%M:%S>')
        data['general']['runtime'] = int((last-first).total_seconds())
        data['general']['status'] = ams.read('General', 'termination status')


    def read_GO():
        GO = {}

        GO['natoms']   = ams.read('InputMolecule', 'nAtoms')
        GO['elements'] = ams.read('InputMolecule', 'AtomSymbols').split()
        GO['steps']    = ams.read('History', 'nEntries')
        GO['energy']   = ams.read('History', f'Energy({GO["steps"]})')

        c = ams.read('InputMolecule', 'Coords')
        GO['input coords'] = [c[i:i+3] for i in range(0, len(c), 3)]
        c = ams.read('History', f'Coords({GO["steps"]})')
        GO['output coords'] = [c[i:i+3] for i in range(0, len(c), 3)]
        
        GO['bond energy']             = adf.read('Energy', 'Bond Energy')
        GO['Mulliken charge']         = adf.read('Properties', 'AtomCharge Mulliken')
        GO['electron dens at nuclei'] = adf.read('Properties', 'Electron Density at Nuclei')
        GO['Hirshfeld charge']        = adf.read('Properties', 'FragmentCharge Hirshfeld')
        GO['Voronoi charge']          = adf.read('Properties', 'AtomCharge_SCF Voronoi')

        warning_lines = [line for line in loglines if 'WARNING:' in line]
        warning_lines = [l for l in warning_lines if 'total elapsed time is much higher than the (CPU+system) time' not in l]
        GO['warnings'] = warning_lines

        error_lines = [line for line in loglines if 'ERROR:' in line]
        GO['errors'] = error_lines
        return GO

    def read_FREQ():
        FREQ = {}

        #vibration data
        FREQ['nmodes'] = adf.read('Vibrations', 'nNormalModes')
        FREQ['frequencies'] = adf.read('Vibrations', 'Frequencies[cm-1]')
        FREQ['intensities'] = adf.read('Vibrations', 'Intensities[km/mol]')
        #correct for case when there is only one vibration (will be float instead of list of floats)
        if type(FREQ['frequencies']) is float: FREQ['frequencies'] = [FREQ['frequencies']]
        if type(FREQ['intensities']) is float: FREQ['intensities'] = [FREQ['intensities']]
        #look for imaginary frequencies
        FREQ['nimag'] = len([f for f in FREQ['frequencies'] if f < 0])
        if FREQ['nimag'] > 0:
            FREQ['imag mode'] = adf.read('Vibrations', 'NoWeightNormalMode(1)')
        FREQ['ZPE'] = adf.read('Vibrations', 'ZeroPointEnergy')

        #thermodynamics data
        FREQ['gibbs'] = adf.read('Thermodynamics', 'Gibbs free Energy')
        FREQ['enthalpy'] = adf.read('Thermodynamics', 'Enthalpy')
        return FREQ

    def read_COSMORS():
        COSMORS = {}
        crskf = plams.KFFile(path('COSMORS_crskf'))
        G = crskf.read('ACTIVITYCOEF', 'deltag')[1]
        COSMORS['gibbs'] = plams.Units.convert(G, 'kcal/mol', 'hartree')
        return COSMORS

    def write_results():
        with open(out_file, 'w+') as res:
            res.write(json.dumps(data, indent=2))

    data = {}
    data['files'] = get_files()

    data['general'] = read_info()

    if data['general']['task'] in ['GO', 'constrained_GO', 'FREQ']:
        ams = plams.KFFile(path('amsrkf'))

        if 'log' in data['files']:
            with open(path('log')) as log: 
                loglines = log.readlines()
        read_log()

        if data['general']['status'] in ['IN PROGRESS']:   
            write_results()
            return

        #adf.rkf is created once the job is done, so we must check status beforehand
        adf = plams.KFFile(path('adfrkf'))

        if data['general']['task'] in ['GO', 'constrained_GO']:
            data['GO'] = read_GO()
        elif data['general']['task'] in ['FREQ']:
            data['FREQ'] = read_FREQ()


    elif data['general']['task'] in ['COSMORS']:
        data['COSMORS'] = read_COSMORS()

    write_results()

    return Result(data)



def get_all_results(calc_dir=paths.calculations2):
    results = []
    for a in os.listdir(calc_dir):
        pa = join(calc_dir, a)
        for b in os.listdir(pa):
            pb = join(pa, b)
            if not os.path.isdir(pb):
                continue
            for c in os.listdir(pb):
                pc = join(pb, c)
                if not os.path.isdir(pc):
                    continue
                if os.path.exists(join(pc, 'SUCCESS')):
                    results.append(generate_result(pc))

    return results



results = get_all_results()
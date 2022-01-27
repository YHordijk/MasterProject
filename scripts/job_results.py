import scm.plams as plams
import paths, os, struct_generator, utility
import datetime, time, json


def get_all_results(calc_path=paths.calculations):
    dirs = []
    for reaction_dir in os.listdir(calc_path):
      p = os.path.join(calc_path, reaction_dir)
      if not os.path.isdir(p): continue
      for calc_dir in os.listdir(p):
          calc_dir = os.path.join(p, calc_dir)
          if os.path.isdir(calc_dir):
              dirs.append(os.path.relpath(calc_dir, paths.master))
    results = []
    for d in dirs:
        results.append(get_results(os.path.join(paths.master, d)))

    return results


def get_results(path):
    results = {}
    results['path'] = path
    results['directory'] = os.path.relpath(path, paths.master)
    results['results_dir'] = os.path.join(paths.results, os.path.relpath(path, paths.calculations))
    results['results_path'] = os.path.join(results['results_dir'], 'results.json')
    os.makedirs(results['results_dir'], exist_ok=True)
    name = os.path.basename(path)
    results['id'] = name.split('.')[0]
    results['stationary_point'] = name.split('.')[3]

    _get_job_files(results)
    _get_runtime(results)
    _read_metainfo(results)
    _get_input_file(results)
    _get_status(results)
    _get_warning(results)
    _get_flags(results)
    _get_natoms(results)
    _get_energies(results)
    _get_normal_modes(results)
    # _write_mols(results)

    mv2s = _make_molview2_start(results, 'out')
    if mv2s is not None: results['files']['molview2_start_out'] = mv2s

    mv2s = _make_molview2_start(results, 'in')
    if mv2s is not None: results['files']['molview2_start_in'] = mv2s

    results['hash'] = utility.hash(results['reaction'], results['stationary_point'], results['flags'])
    results['flags_list'] = results['flags']
    results['flags'] = ' '.join(results['flags'])

    _prepare_for_database(results)
    _write_to_json(results)

    return results


def _get_job_files(results):
    files = {}
    for file in os.listdir(results['path']):
        if file == 'ams.rkf': files['amsrkf'] = os.path.join(results['path'], file)
        if file == 'adf.rkf': files['adfrkf'] = os.path.join(results['path'], file)
        if file == 'output.xyz': files['outxyz'] = os.path.join(results['path'], file)
        if file.endswith('.log'): files['logfile'] = os.path.join(results['path'], file)
        if file.endswith('.err'): files['errfile'] = os.path.join(results['path'], file)
        if file.endswith('.in'): files['infile'] = os.path.join(results['path'], file)
        if file.endswith('.dill'): files['dillfile'] = os.path.join(results['path'], file)
        if file.endswith('.out'): files['outfile'] = os.path.join(results['path'], file)

    files = {n:os.path.relpath(f, paths.master) for n, f in files.items()}
    # for n, f in files.items():
    #     results[n] = f
    results['files'] = files


def _get_status(results):
    if 'amsrkf' in results['files']:
        try:
            rkf = plams.KFFile(os.path.join(paths.master, results['files']['amsrkf']))
            term = rkf.read('General', 'termination status')
            if term == 'IN PROGRESS':
                with open(os.path.join(paths.master, results['files']['errfile'])) as err:
                    err_cont = err.read().strip()
                    if err_cont != '':
                        results['status'] = 'C'; return
                results['status'] = 'R'; return
            elif term == 'NORMAL TERMINATION':
                results['status'] = 'S'; return
            elif 'NORMAL TERMINATION' in term:
                results['status'] = 'W'; return

            results['status'] = 'F'; return
        except:
            results['status'] = 'C'; return
    else:
        results['status'] = 'C'; return


def _get_natoms(results):
    if 'amsrkf' in results['files']:
        rkf = plams.KFFile(os.path.join(paths.master, results['files']['amsrkf']))
        results['natoms'] = rkf.read('Molecule', 'nAtoms')
    else:
        results['natoms'] = 0


def _get_flags(results):
    sorted_Rnames = list(sorted(results['all_substituents'].keys()))
    sorted_R = [results['all_substituents'][R] for R in sorted_Rnames]
    input_file = struct_generator.get_mol_path(paths.input_xyz, results['reaction']+'_'+'_'.join(sorted_R), results['stationary_point'])
    if not os.path.exists(input_file):
        results['flags'] = []
    with open(input_file) as file:
        lines = file.readlines()
        results['flags'] = [f.strip() for f in lines[1].split(',')]

    results['substituents'] = {}
    for f in results['flags']:
        if f.startswith('R'):
            results['substituents'][f.split('=')[0]] = f.split('=')[1].strip()
            results[f.split('=')[0]] = f.split('=')[1].strip()

    results['radical'] = 'radical' in results['flags']
    results['COSMO'] = 'COSMO' in results['flags']

    for task in ['GO', 'SP', 'TSRC', 'LT']:
        if task in results['flags']: 
            results['task'] = task


def _get_input_file(results):
    sorted_Rnames = list(sorted(results['all_substituents'].keys()))
    sorted_R = [results['all_substituents'][R] for R in sorted_Rnames]
    input_file = struct_generator.get_mol_path(paths.input_xyz, results['reaction']+'_'+'_'.join(sorted_R), results['stationary_point'])
    results['files']['inxyz'] = os.path.relpath(input_file, paths.master)
    

def _make_molview2_start(results, stage='out'):
    files = results['files']
    f = stage + 'xyz'
    if f'{stage}xyz' in files:
        with open(os.path.join(paths.master, files[f'{stage}xyz'])) as xyz:
            results[f'{stage}_geometry'] = xyz.readlines()
        molview2_start = os.path.join(results['results_dir'], f'molview2_start_{stage}.bat')
        with open(molview2_start, 'w+') as mvs:
            mvs.write('ECHO ON\n')
            mvs.write(f'cd {paths.scripts}\n')
            mvs.write(f'{paths.driveletter}\n')
            xyz = os.path.join(paths.master, files[f"{stage}xyz"])
            # line = f'python mol_viewer2.py {xyz} n={results["stationary_point"]} reaction={results["reaction"]}'
            # for R, s in results['substituents'].items():
            #     line = line + f' {R}={s}'

            # if 'vibrations' in results and results['vibrations']['freq_nimag'] > 0:
            #     line = line + f' normalmode={",".join(str(d) for d in results["vibrations"]["freq_imag_displ"])}'
            mvs.write(f'python mol_viewer2.py {xyz} stage={stage} res={results["results_path"]}\n')

        return molview2_start


def _get_warning(results):
    warns = []
    if not 'logfile' in results['files']: return []
    logfile = os.path.join(paths.master, results['files']['logfile'])
    with open(logfile) as log:
        for l in log.readlines():
            if ' '.join(l.split()[2:]).startswith('WARNING:'):
                warns.append(' '.join(l.split()[2:]))

    # some warnings can be ignored:
    ignorable_warnings = ['WARNING: total elapsed time is much higher than the (CPU+system) time.']
    results['warnings'] = [w for w in warns if not w in ignorable_warnings]
    if results['status'] == 'W' and len(results['warnings']) == 0:
        results['status'] = 'S'

    results['warnings'] = warns


def _get_runtime(results):
    if not 'logfile' in results['files']: return 0
    logfile = os.path.join(paths.master, results['files']['logfile'])
    with open(logfile) as log:
        lines = log.readlines()
        first = ' '.join(lines[0].split()[0:2])
        last = ' '.join(lines[-1].split()[0:2])
        first = datetime.datetime.strptime(first, '<%b%d-%Y> <%H:%M:%S>')
        last = datetime.datetime.strptime(last, '<%b%d-%Y> <%H:%M:%S>')
        runtime = str(last-first)
        if 'day' in runtime:
            hours = int(runtime.split()[0])*24 + int(runtime.split()[2].split(':')[0])
            mins = int(runtime.split()[2].split(':')[1])
            seconds = int(runtime.split()[2].split(':')[2])
            runtime = f'{hours}:{mins}:{seconds}'
        results['runtime'] = runtime


def _get_energies(results):
    energies = {}
    if 'adfrkf' in results['files']:
        kfreader = plams.KFFile(os.path.join(paths.master, results['files']['adfrkf']))
        energies['bond_energy'] = float(kfreader.read('Energy', 'Bond Energy'))
        energies['eda_elstat'] = float(kfreader.read('Energy', 'elstat'))
        energies['eda_pauli'] = float(kfreader.read('Energy', 'Pauli Total'))
        energies['eda_orbit_int'] = float(kfreader.read('Energy', 'Orb.Int. Total'))
        energies['gibbs_energy'] = float(kfreader.read('Thermodynamics', 'Gibbs free Energy'))
    results['energies'] = energies


def _get_normal_modes(results):
    if 'adfrkf' in results['files']:
        results['vibrations'] = {}
        kfreader = plams.KFFile(os.path.join(paths.master, results['files']['adfrkf']))
        freqs = kfreader.read('Vibrations', 'Frequencies[cm-1]')
        if type(freqs) is float: freqs = [freqs]
        ints = kfreader.read('Vibrations', 'Intensities[km/mol]')
        if type(ints) is float: ints = [ints]
        results['vibrations']['frequencies'] = freqs
        results['vibrations']['intensities'] = ints
        results['vibrations']['freq_nimag'] = len([f for f in freqs if f<0])
        if results['vibrations']['freq_nimag'] > 0:
            results['vibrations']['freq_imag_displ'] = [float(x) for x in kfreader.read('Vibrations', 'NoWeightNormalMode(1)')]
        else:
            results['vibrations']['freq_imag_displ'] = []


def _read_metainfo(results):
    results['all_substituents'] = {}
    results['files']['metainfo'] = os.path.join(os.path.dirname(results['path']), 'meta.info')
    with open(results['files']['metainfo']) as meta:
        lines = meta.readlines()
        for line in lines:
            if line.startswith('reaction='):
                results['reaction'] = line.split('=')[1].strip()
            if line.startswith('R'):
                results['all_substituents'][line.split('=')[0]] = line.split('=')[1].strip()


def _write_to_json(results):
    p = results['results_path']
    with open(p, 'w+') as resfile:
        json.dump(results, resfile, indent='\t', sort_keys=True)






def _prepare_for_database(results):
    fields = ['id',
          'status',
          'task',
          'reaction',
          'R1', 'R2', 'Rcat',
          'stationary_point',
          'calc_directory',
          'result_directory',
          'flags',
          'inxyz',
          'outxyz',
          'runtime',
          'hash',
          ]
    results['database_fields'] = fields

    d = [
    results['id'],
    results['status'],
    results['task'],
    results['reaction'],
    results['substituents'].get('R1', None),
    results['substituents'].get('R2', None),
    results['substituents'].get('Rcat', None),
    results['stationary_point'],
    results['directory'],
    results['results_dir'],
    results['flags'],
    results['files'].get('inxyz', None),
    results['files'].get('outxyz', None),
    results.get('runtime', None),
    results['hash'],
        ]
    results['database_data'] = d


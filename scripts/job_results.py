import scm.plams as plams
import paths, os, struct_generator, utility
import datetime, time


def get_results(path):
    results = {}
    files = _get_job_files(path)
    files = {n:os.path.relpath(f, paths.master) for n, f in files.items()}
    results['files'] = files
    name = os.path.basename(path)
    results['directory'] = os.path.relpath(path, paths.master)

    results['id'] = name.split('.')[0]
    results['stationary_point'] = name.split('.')[3]

    results['metainfo'] = os.path.join(os.path.dirname(path), 'meta.info')

    results['runtime'] = _get_runtime(results)

    results['all_substituents'] = {}
    with open(results['metainfo']) as meta:
        lines = meta.readlines()
        for line in lines:
            if line.startswith('reaction='):
                results['reaction'] = line.split('=')[1].strip()
            if line.startswith('R'):
                results['all_substituents'][line.split('=')[0]] = line.split('=')[1].strip()

    results['files']['inxyz'] = _get_input_file(results)
    
    for n, f in files.items():
        results[n] = f
    results['status'] = _get_status(files)
    warn = _get_warning(results)
    # some warnings can be ignored:
    ignorable_warnings = ['WARNING: total elapsed time is much higher than the (CPU+system) time.']
    results['warnings'] = [w for w in warn if not w in ignorable_warnings]
    if results['status'] == 'W' and len(results['warnings']) == 0:
        results['status'] = 'S'


    flags = _get_flags(results)
    results['flags'] = ' '.join(flags)

    results['substituents'] = {}
    for f in results['flags'].split():
        if f.startswith('R'):
            results['substituents'][f.split('=')[0]] = f.split('=')[1].strip()
            results[f.split('=')[0]] = f.split('=')[1].strip()

    results['radical'] = 'radical' in flags
    results['COSMO'] = 'COSMO' in flags
    results['natoms'] = _get_natoms(files)

    mv2s = _make_molview2_start(results, 'out')
    if mv2s is not None: results['molview2_start_out'] = mv2s

    mv2s = _make_molview2_start(results, 'in')
    if mv2s is not None: results['molview2_start_in'] = mv2s

    for task in ['GO', 'SP', 'TSRC', 'LT']:
        if task in flags: 
            results['task'] = task

    results['hash'] = utility.hash(results['reaction'], results['stationary_point'], results['flags'].split())
    return results


def _get_job_files(path):
    files = {}
    for file in os.listdir(path):
        if file == 'ams.rkf': files['amsrkf'] = os.path.join(path, file)
        if file == 'adf.rkf': files['adfrkf'] = os.path.join(path, file)
        if file == 'output.xyz': files['outxyz'] = os.path.join(path, file)
        if file.endswith('.log'): files['logfile'] = os.path.join(path, file)
        if file.endswith('.err'): files['errfile'] = os.path.join(path, file)
        if file.endswith('.in'): files['infile'] = os.path.join(path, file)
        if file.endswith('.dill'): files['dillfile'] = os.path.join(path, file)
        if file.endswith('.out'): files['outfile'] = os.path.join(path, file)
    return files


def _get_status(files):
    if 'amsrkf' in files:
        rkf = plams.KFFile(os.path.join(paths.master, files['amsrkf']))
        term = rkf.read('General', 'termination status')
        if term == 'IN PROGRESS':
            with open(os.path.join(paths.master, files['errfile'])) as err:
                err_cont = err.read().strip()
                if err_cont != '':
                    return 'C'
            return 'R'
        elif term == 'NORMAL TERMINATION':
            return 'S'
        elif 'NORMAL TERMINATION' in term:
            return 'W'

        return 'F'
    else:
        return 'C'

def _get_natoms(files):
    if 'amsrkf' in files:
        rkf = plams.KFFile(os.path.join(paths.master, files['amsrkf']))
        return rkf.read('Molecule', 'nAtoms')
    else:
        return 0


def _get_flags(results):
    sorted_Rnames = list(sorted(results['all_substituents'].keys()))
    sorted_R = [results['all_substituents'][R] for R in sorted_Rnames]
    input_file = struct_generator.get_mol_path(paths.input_xyz, results['reaction']+'_'+'_'.join(sorted_R), results['stationary_point'])
    if not os.path.exists(input_file): 
        print(input_file)
        return []
    with open(input_file) as file:
        lines = file.readlines()
        return [f.strip() for f in lines[1].split(',')]


def _get_input_file(results):
    sorted_Rnames = list(sorted(results['all_substituents'].keys()))
    sorted_R = [results['all_substituents'][R] for R in sorted_Rnames]
    input_file = struct_generator.get_mol_path(paths.input_xyz, results['reaction']+'_'+'_'.join(sorted_R), results['stationary_point'])
    return input_file
    

def _make_molview2_start(results, stage='out'):
    files = results['files']
    f = stage + 'xyz'
    if f'{stage}xyz' in files:
        results['geometry'] = os.path.relpath(files[f'{stage}xyz'], results['directory'])
        molview2_start = os.path.join(paths.master, results['directory'], f'molview2_start_{stage}.bat')
        with open(os.path.join(paths.master, results['directory'], f'molview2_start_{stage}.bat'), 'w+') as mvs:
            mvs.write('ECHO ON\n')
            mvs.write(f'cd {paths.scripts}\n')
            mvs.write(f'{paths.driveletter}\n')
            xyz = os.path.join(paths.master, files[f"{stage}xyz"])
            line = f'python mol_viewer2.py {xyz} n={results["stationary_point"]} reaction={results["reaction"]}'
            for R, s in results['substituents'].items():
                line = line + f' {R}={s}'
            mvs.write(line)

        return molview2_start

def _get_warning(results):
    warns = []
    if not 'logfile' in results['files']: return []
    logfile = os.path.join(paths.master, results['files']['logfile'])
    with open(logfile) as log:
        for l in log.readlines():
            if ' '.join(l.split()[2:]).startswith('WARNING:'):
                warns.append(' '.join(l.split()[2:]))
    return warns

def _get_runtime(results):
    if not 'logfile' in results['files']: return 0
    logfile = os.path.join(paths.master, results['files']['logfile'])
    with open(logfile) as log:
        lines = log.readlines()
        first = ' '.join(lines[0].split()[0:2])
        last = ' '.join(lines[-1].split()[0:2])
        first = datetime.datetime.strptime(first, '<%b%d-%Y> <%H:%M:%S>')
        last = datetime.datetime.strptime(last, '<%b%d-%Y> <%H:%M:%S>')
        return str(last-first)


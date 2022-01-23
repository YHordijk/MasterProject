import scm.plams as plams
import paths, os, struct_generator, utility


def get_results(path):
    results = {}
    files = _get_job_files(path)
    files = {n:os.path.relpath(f, paths.master) for n, f in files.items()}
    name = os.path.basename(path)
    results['id'] = name.split('.')[0]
    results['metainfo'] = os.path.join(os.path.dirname(path), 'meta.info')
    results['substituents'] = {}
    with open(results['metainfo']) as meta:
        lines = meta.readlines()
        for line in lines:
            if line.startswith('reaction='):
                results['reaction'] = line.split('=')[1].strip()
            if line.startswith('R'):
                results['substituents'][line.split('=')[0]] = line.split('=')[1].strip()

    results['stationary_point'] = name.split('.')[3]
    results['directory'] = os.path.relpath(path, paths.master)
    results['files'] = files
    for n, f in files.items():
        results[n] = f
    results['status'] = _get_status(files)
    flags = _get_flags(results)
    results['flags'] = ' '.join(flags)

    results['radical'] = 'radical' in flags
    results['COSMO'] = 'COSMO' in flags
    results['natoms'] = _get_natoms(files)

    if 'outxyz' in files:
        results['geometry'] = os.path.relpath(files['outxyz'], results['directory'])
        results['molview2_start'] = os.path.join(paths.master, path, 'molview2_start.bat')
        with open(os.path.join(paths.master, path, 'molview2_start.bat'), 'w+') as mvs:
            mvs.write('ECHO ON\n')
            mvs.write(f'cd {paths.scripts}\n')
            mvs.write(f'{paths.driveletter}\n')
            mvs.write(f'python mol_viewer2.py {os.path.join(paths.master, files["outxyz"])}')

    for task in ['GO', 'SP', 'TSRC', 'LT']:
        if task in flags: results['task'] = task
    for f in flags:
        if f.startswith('R'):
            name, sub = f.split('=')
            results[name] = sub



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
    sorted_Rnames = list(sorted(results['substituents'].keys()))
    sorted_R = [results['substituents'][R] for R in sorted_Rnames]
    input_file = struct_generator.get_mol_path(paths.input_xyz, results['reaction']+'_'+'_'.join(sorted_R), results['stationary_point'])
    if not os.path.exists(input_file): return []
    with open(input_file) as file:
        lines = file.readlines()
        return [f.strip() for f in lines[1].split(',')]
    





# get_results(r"D:\Users\Yuman\Desktop\MasterProject\calculations\no_catalyst_H_H\3.no_catalyst_H_H.TS")
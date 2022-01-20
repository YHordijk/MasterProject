import scm.plams as plams
import paths, os, struct_generator

def get_results(path):
    results = {}
    files = _get_job_files(path)
    files = {n:os.path.relpath(f, paths.master) for n, f in files.items()}
    name = os.path.basename(path)
    results['id'] = name.split('.')[0]
    results['reaction'] = name.split('.')[1]
    results['stationary_point'] = name.split('.')[2]
    results['directory'] = os.path.relpath(path, paths.master)
    results['files'] = files
    for n, f in files.items():
        results[n] = f
    results['status'] = _get_status(files)
    flags = _get_flags(results)
    results['flags'] = flags
    results['radical'] = 'radical' in flags
    for task in ['GO', 'SP', 'TSRC', 'LT']:
        if task in flags: results['task'] = task
    for f in flags:
        if f.startswith('R'):
            name, sub = f.split('=')
            results[name] = sub
    print(results)
    return results


def _get_job_files(path):
    files = {}
    for file in os.listdir(path):
        if file == 'ams.rkf': files['amsrkf'] = os.path.join(path, file)
        if file == 'adf.rkf': files['adfrkf'] = os.path.join(path, file)
        if file.endswith('.log'): files['logfile'] = os.path.join(path, file)
        if file.endswith('.err'): files['errfile'] = os.path.join(path, file)
        if file.endswith('.in'): files['infile'] = os.path.join(path, file)
        if file.endswith('.dill'): files['dillfile'] = os.path.join(path, file)
        if file.endswith('.out'): files['outfile'] = os.path.join(path, file)
    return files


def _get_status(files):
    rkf = plams.KFFile(os.path.join(paths.master, files['amsrkf']))
    term = rkf.read('General', 'termination status')
    if term == 'IN PROGRESS':
        return 'R'
    elif term == 'NORMAL TERMINATION':
        return 'S'
    elif 'NORMAL TERMINATION' in term:
        return 'W'

    return 'F'


def _get_flags(results):
    input_file = struct_generator.get_mol_path(paths.input_xyz, results['reaction'], results['stationary_point'])
    with open(input_file) as file:
        lines = file.readlines()
        return [f.strip() for f in lines[1].split(',')]
    





# get_results(r"D:\Users\Yuman\Desktop\MasterProject\calculations\no_catalyst_H_H\3.no_catalyst_H_H.TS")
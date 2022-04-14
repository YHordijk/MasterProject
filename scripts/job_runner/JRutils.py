import utility, os, paths
join = os.path.join


def info_file_parser(file):
    data = {}
    with open(file) as info:
        for line in info.readlines():
            arg, val = line.strip().split('=')
            data[arg] = val

    return data


def get_run_dir(calc_path, job_name, mol_info, ignore_status=True):
    existing_job_dirs = [join(calc_path, d) for d in os.listdir(calc_path) if os.path.isdir(join(calc_path, d))]
    hash = utility.hash2(get_unique(mol_info))
    collisions = [check_collision(hash, d) for d in existing_job_dirs]
    successes = ['SUCCESS' in os.listdir(d) for d in existing_job_dirs]

    for c, s in zip(collisions, successes):
        if c:
            if ignore_status:
                return
            elif s:
                return

    i = 1
    base = join(calc_path, job_name)
    while os.path.exists(base + f'.{str(i).zfill(3)}'):
        i += 1
    return base + f'.{str(i).zfill(3)}'


def check_collision(h1, path):
    with open(join(path, 'run.info')) as info:
        for line in info.readlines():
            arg, val = line.strip().split('=')
            if arg == 'unique':
                h2 = utility.hash2(val)

    return h1 == h2


def get_unique(run_info):
    features = [
        run_info['reaction'],
        run_info['phase'],
        run_info['task'],
        run_info['stationary_point'],
        run_info["sorted_substituents"],
        run_info['TSRC_idx'],
        run_info['radical'],
        run_info['enantiomer'],
        run_info['functional'],
        run_info['basis'],
        run_info['numerical_quality']
        ]
    return '.'.join(features)



def get_job_runner_template(task):
    file = join(paths.JR_templates, task)
    assert os.path.exists(file)
    return file



def read_mol_info(calc_path):
    data = info_file_parser(join(calc_path, 'mol.info'))
    for arg in ['plane_idx', 'align_idx']:
        if arg in data:
            data[arg] = [int(i) for i in data[arg].split("_")]
    for arg in ['center_idx']:
        if arg in data:
            data[arg] = int(data[arg])
    return data
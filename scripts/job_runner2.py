import scm.plams as plams
import struct_generator2, paths, os, utility


join = os.path.join



def run_jobs(template, substituents={}, calc_dir=paths.calculations, phase='vacuum', test_mode=False):
    def get_mols():
        print(f'Generating molecules:')
        print(f'\tTemplate: {template}')
        for name, sub in substituents.items():
            print(f'\t{name:<8}: {sub}')
        mols = struct_generator2.generate_stationary_points(template, substituents)
        print(f'\nGenerated {len(mols)} molecules:')
        struct_generator2.print_mols(mols, tabs=1)
        return mols


    def get_job_runner_template(mol):
        file = join(paths.JR_templates, mol.task)
        assert os.path.exists(file)
        return file


    def sort_substituents(subs):
        sorted_Rnames = list(sorted(subs))
        sorted_R = [substituents[R] for R in sorted_Rnames]
        return sorted_R


    def get_job_dir(mol):
        sorted_R = sort_substituents(substituents)
        return join(calc_dir, template + '.' + '_'.join(sorted_R) + '.' + mol.phase, mol.name)


    def write_runscript(mol):
        RADICAL = ''
        if mol.radical:
            RADICAL = '\n  SpinPolarization 1\n  Unrestricted Yes\n'
        ATOMS = ''
        for atom in mol.atoms:
            ATOMS += f'    {atom.symbol:<2}\t{atom.coords[0]:=11.8f}\t{atom.coords[1]:=11.8f}\t{atom.coords[2]:=11.8f}\n'
        TSRC = ''
        if mol.task == 'TSRC':
            TSRC = f'    Distance {mol.TSRC_idx[0]} {mol.TSRC_idx[1]} -1\n'
        SLURM_SUBMIT_DIR = runscript_path(mol)

        CORES = '64'
        NODES = '1'
        if mol.radical:
            CORES = '64'
            NODES = '2'

        print(mol.name, NODES, mol.radical)

        blocks = {
                '[RADICAL]':          RADICAL, 
                '[ATOMS]':            ATOMS, 
                '[TSRC]':             TSRC, 
                '[SLURM_SUBMIT_DIR]': SLURM_SUBMIT_DIR,
                '[CORES]':            CORES,
                '[NODES]':            NODES,
                }

        with open(JR_template) as temp:
            temp_lines = temp.readlines()
            with open(runscript_path(mol), 'w+') as run:
                for line in temp_lines:
                    for block in blocks:
                        if block in line:
                            line = line.replace(block, blocks[block])

                    run.write(line)


    def get_unique(mol):
        return f'{mol.reaction}.{mol.phase}.{mol.task}.{mol.name}.{sort_substituents(mol.substituents)}.{sorted(mol.TSRC_idx)}.{mol.radical}.{mol.enantiomer}'


    def hash(mol):
        return utility.hash2(get_unique(mol))


    def write_info(mol):
        #this function writes information used by crawlers later on
        info_path = join(os.path.dirname(runscript_path(mol)), 'run.info')
        with open(info_path, 'w+') as info:
            info.write(f'unique={get_unique(mol)}\n')
            info.write(f'reaction={template}\n')
            info.write(f'stationary_point={mol.name}\n')
            info.write(f'task={mol.task}\n')
            info.write(f'TSRC_idx={"_".join(str(i) for i in mol.TSRC_idx)}\n' if len(mol.TSRC_idx) > 0 else 'TSRC=N\\A\n')
            info.write(f'radical={mol.radical}\n')
            info.write(f'phase={mol.phase}\n')
            info.write(f'sorted_substituents={"_".join(sort_substituents(mol.substituents))}\n')
            for sub in mol.substituents:
                info.write(f'{sub}={substituents[sub]}\n')
            info.write(f'enantiomer={mol.enantiomer}\n')
            

    def runscript_path(mol):
        return join(job_dir, f'{mol.task}.{template}.' + '_'.join(sort_substituents(substituents)) + '.' + mol.name)


    def run_job(file):
        start_script = f"""#!/bin/bash
cd {os.path.dirname(file)}
sbatch {os.path.basename(file)}
"""     
        os.system(start_script)

    mols = get_mols()
    i = 0
    for mol_name, mol in mols.items():
        i += 1
        # if mol_name != 'P2': continue
        mol.phase = phase

        if utility.hash_collision(hash(mol), calc_dir):
            print(f'Hash collision detected {i} {mol_name}')
            continue

        JR_template = get_job_runner_template(mol)
        job_dir = get_job_dir(mol)
        if not test_mode:
            os.makedirs(job_dir, exist_ok=True)

            #copy input molecule to new dir
            mol.write(join(job_dir, 'input.xyz'))
            #make job_run script
            write_runscript(mol)
            write_info(mol)
        #run the job we just made
        print(f'Running job {os.path.relpath(runscript_path(mol), calc_dir)}')
        if not test_mode:
            run_job(runscript_path(mol))

        
if __name__ == '__main__':
    # for R2 in ['Ph', 'tBu']:
    #     for cat in ['I2', 'ZnCl2', 'TiCl4', 'BF3', 'AlF3', 'SnCl4']:
    #         run_jobs('achiral_catalyst', {'R1':'H', 'R2':R2, 'Rcat':cat}, phase='vacuum', calc_dir=join(paths.master, 'calculations_test'), test_mode=False)

    #     run_jobs('no_catalyst', {'R1':'H', 'R2':R2}, phase='vacuum', calc_dir=join(paths.master, 'calculations_test'), test_mode=False)
    run_jobs('squaramide', {'R1':'H', 'R2':'Ph', 'Rc1':'H', 'Rc2':'Ph', 'Rch':'O', 'Rch2':'O'}, phase='vacuum', calc_dir=join(paths.master, 'calculations_test'), test_mode=False)
    run_jobs('achiral_catalyst', {'R1':'H', 'R2':'H', 'Rcat':'SnCl4'}, phase='vacuum', calc_dir=join(paths.master, 'calculations_test'), test_mode=False)
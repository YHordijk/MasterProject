import scm.plams as plams
import struct_generator2, paths, os, utility, job_results3, time, sys, geometry
from utility import bohr2angstrom
b2a = bohr2angstrom(1)

join = os.path.join



def run_frag_job(calc_path):
    def run_job(file):
        start_script = f"""#!/bin/bash
cd {os.path.dirname(file)}
sbatch {os.path.basename(file)}
"""     
        os.system(start_script)


    #read information about run
    substituents = {}
    plane_idx = None
    align_idx = None
    center_idx = None
    frag1idx = None
    frag2idx = None
    with open(join(calc_path, 'run.info')) as info:
        for line in info.readlines():
            line = line.strip()
            key, value = line.split('=')
            if key.startswith('R'):        substituents[key] = value
            if key == 'reaction':          template            = value
            if key == 'functional':        functional          = value
            if key == 'basis':             basis               = value
            if key == 'phase':             phase               = value
            if key == 'stationary_point':  stationary_point    = value
            if key == 'radical':           radical             = value
            if key == 'numerical_quality': numerical_quality   = value
            if key == 'plane_idx':         plane_idx           = [int(i) for i in value.split('_')]
            if key == 'align_idx':         align_idx           = [int(i) for i in value.split('_')]
            if key == 'center_idx':        center_idx          = int(value)
            if key == 'frag1idx':          frag1idx            = [int(i) for i in value.split('_')]
            if key == 'frag2idx':          frag2idx            = [int(i) for i in value.split('_')]

    if functional == 'OLYP': return
    # if frag1idx is None or frag2idx is None: return
    molfile = join(calc_path, 'output.xyz')
    align_mol = plams.Molecule(molfile)
    geometry.align_mol(align_mol, plane_idx=plane_idx, align_idx=align_idx, center_idx=center_idx)

    ATOMS = ''
    for i, atom in enumerate(align_mol.atoms, 1):
        if i in frag1idx:
            ATOMS += f'    {atom.symbol:<2}\t{atom.coords[0]*b2a:=11.8f}\t{atom.coords[1]*b2a:=11.8f}\t{atom.coords[2]*b2a:=11.8f}\tf=f1\n'
        elif i in frag2idx:
            ATOMS += f'    {atom.symbol:<2}\t{atom.coords[0]*b2a:=11.8f}\t{atom.coords[1]*b2a:=11.8f}\t{atom.coords[2]*b2a:=11.8f}\tf=f2\n'

    ATOMS1 = ''
    for i, atom in enumerate(align_mol.atoms, 1):
        if i in frag1idx:
            ATOMS1 += f'    {atom.symbol:<2}\t{atom.coords[0]*b2a:=11.8f}\t{atom.coords[1]*b2a:=11.8f}\t{atom.coords[2]*b2a:=11.8f}\tregion=Region_1\n'

    ATOMS2 = ''
    for i, atom in enumerate(align_mol.atoms, 1):
        if i in frag2idx:
            ATOMS2 += f'    {atom.symbol:<2}\t{atom.coords[0]*b2a:=11.8f}\t{atom.coords[1]*b2a:=11.8f}\t{atom.coords[2]*b2a:=11.8f}\tregion=Region_2\n'

    FUNCTIONAL = {
            'OLYP': 'GGA OLYP',
            'BLYP-D3(BJ)': 'GGA BLYP\n    DISPERSION GRIMME3 BJDAMP',
            'M06L': 'MetaGGA M06L'
        }[functional]

    RADICAL = ''
    if radical == 'True':
        RADICAL = '\n  SpinPolarization 1\n  Unrestricted Yes\n'
    CORES = '64'
    NODES = '1'
    
    blocks = {
        '[RADICAL]':          RADICAL, 
        '[ATOMS]':            ATOMS, 
        '[ATOMS1]':           ATOMS1, 
        '[ATOMS2]':           ATOMS2, 
        '[CORES]':            CORES,
        '[NODES]':            NODES,
        '[BASIS]':            basis,
        '[FUNCTIONAL]':       FUNCTIONAL,
        '[NUMERICALQUALITY]': numerical_quality,
        }

    JR_template = join(paths.JR_templates, 'FRAG')
    runscript = join(calc_path, 'frag.run')
    with open(JR_template, 'r') as temp:
        temp_lines = temp.readlines()
        with open(runscript, 'w+') as run:
            for line in temp_lines:
                for block in blocks:
                    if block in line:
                        line = line.replace(block, blocks[block])

                run.write(line)

    run_job(runscript)
    # print('running?', runscript)


# with open(JR_template) as temp:
#             temp_lines = temp.readlines()
#             with open(runscript_path(mol), 'w+') as run:
#                 for line in temp_lines:
#                     for block in blocks:
#                         if block in line:
#                             line = line.replace(block, blocks[block])

#                     run.write(line)


def run_SP_job(calc_path):
    def run_job(file):
        start_script = f"""#!/bin/bash
cd {os.path.dirname(file)}
sbatch {os.path.basename(file)}
"""     
        os.system(start_script)

    #read information about run
    substituents = {}
    plane_idx = None
    align_idx = None
    center_idx = None
    with open(join(calc_path, 'run.info')) as info:
        for line in info.readlines():
            line = line.strip()
            key, value = line.split('=')
            if key.startswith('R'):        substituents[key] = value
            if key == 'reaction':          template            = value
            if key == 'functional':        functional          = value
            if key == 'basis':             basis               = value
            if key == 'phase':             phase               = value
            if key == 'stationary_point':  stationary_point    = value
            if key == 'radical':           radical             = value
            if key == 'numerical_quality': numerical_quality   = value
            if key == 'plane_idx':         plane_idx           = [int(i) for i in value.split('_')]
            if key == 'align_idx':         align_idx           = [int(i) for i in value.split('_')]
            if key == 'center_idx':        center_idx          = int(value)

    if functional == 'OLYP': return
    molfile = join(calc_path, 'output.xyz')
    align_mol = plams.Molecule(molfile)
    geometry.align_mol(align_mol, plane_idx=plane_idx, align_idx=align_idx, center_idx=center_idx)
    ATOMS = ''
    for i, atom in enumerate(align_mol.atoms, 1):
        ATOMS += f'    {atom.symbol:<2}\t{atom.coords[0]*b2a:=11.8f}\t{atom.coords[1]*b2a:=11.8f}\t{atom.coords[2]*b2a:=11.8f}\n'

    FUNCTIONAL = {
            'OLYP': 'GGA OLYP',
            'BLYP-D3(BJ)': 'GGA BLYP\n    DISPERSION GRIMME3 BJDAMP',
            'M06L': 'MetaGGA M06L'
        }[functional]

    RADICAL = ''
    if radical == 'True':
        RADICAL = '\n  SpinPolarization 1\n  Unrestricted Yes\n'
    CORES = '64'
    NODES = '1'
    
    blocks = {
        '[RADICAL]':          RADICAL, 
        '[ATOMS]':            ATOMS, 
        '[CORES]':            CORES,
        '[NODES]':            NODES,
        '[BASIS]':            basis,
        '[FUNCTIONAL]':       FUNCTIONAL,
        '[NUMERICALQUALITY]': numerical_quality,
        }

    JR_template = join(paths.JR_templates, 'SP')
    runscript = join(calc_path, 'SP.run')
    with open(JR_template, 'r') as temp:
        temp_lines = temp.readlines()
        with open(runscript, 'w+') as run:
            for line in temp_lines:
                for block in blocks:
                    if block in line:
                        line = line.replace(block, blocks[block])

                run.write(line)

    run_job(runscript)
    # print('running?', runscript)



def run_jobs(template, substituents={}, calc_dir=paths.calculations, phase='vacuum', 
                test_mode=False, basis='DZP', functional='BLYP-D3(BJ)', numerical_quality='Good'):

    assert basis in ['SZ', 'DZ', 'DZP', 'TZP', 'TZ2P', 'QZ4P']
    assert functional in ['OLYP', 'BLYP-D3(BJ)', 'M06L']
    assert numerical_quality in ['Basic', 'Normal', 'Good', 'VeryGood', 'Excellent']

    def get_mols():
        print(f'Generating molecules:')
        print(f'\tTemplate  : {template}')
        for name, sub in substituents.items():
            print(f'\t{name:<8}  : {sub}')
        print(f'\tBasis     : {basis}')
        print(f'\tFunctional: {functional}')
        print(f'\tQuality   : {numerical_quality}')
        mols = struct_generator2.generate_stationary_points(template, substituents)
        print(f'\nGenerated {len(mols)} molecules:')
        struct_generator2.print_mols(mols, tabs=1)
        print()
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
        if mol.radical and mol.reaction in ['urea_tBu_Ph', 'squaramide']:
            CORES = '64'
            NODES = '2'

        FUNCTIONAL = {
            'OLYP': 'GGA OLYP',
            'BLYP-D3(BJ)': 'GGA BLYP\n    DISPERSION GRIMME3 BJDAMP',
            'M06L': 'MetaGGA M06L'
        }[functional]


        blocks = {
                '[RADICAL]':          RADICAL, 
                '[ATOMS]':            ATOMS, 
                '[TSRC]':             TSRC, 
                '[SLURM_SUBMIT_DIR]': SLURM_SUBMIT_DIR,
                '[CORES]':            CORES,
                '[NODES]':            NODES,
                '[BASIS]':            basis,
                '[FUNCTIONAL]':       FUNCTIONAL,
                '[NUMERICALQUALITY]': numerical_quality,
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
        return f'{mol.reaction}.{mol.phase}.{mol.task}.{mol.name}.{sort_substituents(mol.substituents)}.{sorted(mol.TSRC_idx)}.{mol.radical}.{mol.enantiomer}.{functional}.{basis}.{numerical_quality}'


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
            info.write(f'functional={functional}\n')
            info.write(f'basis={basis}\n')
            info.write(f'numerical_quality={numerical_quality}\n')
            if hasattr(mol, 'active_atom_idx'):
                info.write(f'active_atom_idx={mol.active_atom_idx}\n')
            if hasattr(mol, 'plane_idx'):
                info.write(f'plane_idx={"_".join([str(i) for i in mol.plane_idx])}\n')
            if hasattr(mol, 'align_idx'):
                info.write(f'align_idx={"_".join([str(i) for i in mol.align_idx])}\n')
            if hasattr(mol, 'center_idx'):
                info.write(f'center_idx={mol.center_idx}\n')
            if hasattr(mol, 'frag1idx'):
                info.write(f'frag1idx={"_".join([str(i) for i in mol.frag1idx])}\n')
            if hasattr(mol, 'frag2idx'):
                info.write(f'frag2idx={"_".join([str(i) for i in mol.frag2idx])}\n')
            

    def runscript_path(mol):
        return join(job_dir, f'{mol.task}.{template}.' + '_'.join(sort_substituents(substituents)) + '.' + mol.name)


    def run_job(file):
        start_script = f"""#!/bin/bash
cd {os.path.dirname(file)}
sbatch {os.path.basename(file)}
"""     
        os.system(start_script)

    mols = get_mols()
    Njobs = 0
    i = 0
    for mol_name, mol in mols.items():
        i += 1
        # if mol_name != 'P2': continue
        mol.phase = phase

        if utility.hash_collision(hash(mol), calc_dir):
            print(f'Hash collision detected {i} {mol_name}')
            continue
        Njobs += 1

        JR_template = get_job_runner_template(mol)
        job_dir = get_job_dir(mol)
        if not test_mode:
            #check if the dir already exists
            try:
                os.makedirs(job_dir, exist_ok=False)
            except:
                #generate new dir in PLAMS style
                i = 2
                job_dir_next = job_dir + '.' + str(i).zfill(3)
                while os.path.exists(job_dir_next):
                    i += 1
                    job_dir_next = job_dir + '.' + str(i).zfill(3)
                job_dir = job_dir_next
                os.makedirs(job_dir, exist_ok=False)

            #copy input molecule to new dir
            mol.write(join(job_dir, 'input.xyz'))
            #make job_run script
            write_runscript(mol)
            write_info(mol)
        #run the job we just made
        print(f'Running job {os.path.relpath(runscript_path(mol), calc_dir)}')
        if not test_mode:
            run_job(runscript_path(mol))

    return Njobs


    
if __name__ == '__main__':
    mode = 'GO'
    calc_dir = paths.calculations
    test_mode = True


    n = 0
    if mode == 'GO':
        basis = 'TZ2P'
        functional = 'BLYP-D3(BJ)'
        numerical_quality = 'Good'

        for R1 in ['Et', 'NMe2']:
            for R2 in ['Bz', 'NHMe', 'OEt', 'p-HOPh', 'm-HOPh', 'o-HOPh']:
                for cat in ['I2', 'ZnCl2', 'TiCl4', 'BF3', 'AlF3', 'SnCl4']:
                    n += run_jobs('achiral_catalyst', {'R1':R1, 'R2':R2, 'Rcat':cat}, phase='vacuum', calc_dir=calc_dir, test_mode=test_mode, basis=basis, functional=functional, numerical_quality=numerical_quality)
                    # time.sleep(1)

        # for R1  in ['NH2', 'Me', 'OMe']:
        #     for R2 in ['Et', 'NMe2']:
        #         for Rch in ['S']:
        #             n += run_jobs('urea_tBu_Ph', {'R1':R1, 'R2':R2, 'Rch':Rch}, phase='vacuum', calc_dir=calc_dir, test_mode=test_mode, basis=basis, functional=functional, numerical_quality=numerical_quality)
        #             time.sleep(2)

    if mode == 'SP':
        filtered_dirs = []
        dirs = job_results3.get_all_run_dirs() 
        #filter sub_cat_complex stationary points for achiral_catalysts
        filtered_dirs += [d for d in dirs if 'achiral_catalyst' in d and os.path.basename(d) in ['sub_cat_complex.002', 'sub_cat_complex']]
        n = 0
        for d in filtered_dirs:
            # if not os.path.exists(join(d, 'SP.run.out')):
                print(d)
                if not test_mode:
                    run_SP_job(d)
                    time.sleep(1)
                n += 1
                

    if mode == 'fragment':
        filtered_dirs = []
        dirs = job_results3.get_all_run_dirs() 
        #filter sub_cat_complex stationary points for achiral_catalysts
        filtered_dirs += [d for d in dirs if 'achiral_catalyst' in d and os.path.basename(d) in ['sub_cat_complex.002', 'sub_cat_complex']]
        filtered_dirs = ['/scistor/tc/yhk800/MasterProject/calculations/achiral_catalyst.H_Ph_AlF3.vacuum/sub_cat_complex',
                         '/scistor/tc/yhk800/MasterProject/calculations/achiral_catalyst.H_Ph_I2.vacuum/sub_cat_complex',
                         '/scistor/tc/yhk800/MasterProject/calculations/achiral_catalyst.H_tBu_BF3.vacuum/sub_cat_complex',
                         '/scistor/tc/yhk800/MasterProject/calculations/achiral_catalyst.H_Ph_BF3.vacuum/sub_cat_complex',
                         '/scistor/tc/yhk800/MasterProject/calculations/achiral_catalyst.H_Ph_SnCl4.vacuum/sub_cat_complex']
        n = 0
        for d in filtered_dirs:
            # if not os.path.exists(join(d, 'frag.run.out')):
                print(d)
                if not test_mode:
                    run_frag_job(d)
                    time.sleep(1)
                n += 1
                
            
    print(n)
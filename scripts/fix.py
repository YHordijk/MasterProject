import struct_generator2, job_results3
import os, shutil

for d in job_results3.get_all_run_dirs():
    if not 'no_catalyst' in d: continue
    runinfo = os.path.join(d, 'run.info')
    runinfoout = os.path.join(d, 'run.info.test')
    try:
        mol = struct_generator2.get_mol_from_runinfo(runinfo)
    except: continue
    with open(runinfo, 'r') as info:
        info_lines = info.readlines()
        keys = [l.strip().split('=')[0] for l in info_lines]
        values = [l.strip().split('=')[1] for l in info_lines]
        data = {k:v for k, v in zip(keys, values)}

        Rs = {R:sub for R, sub in zip(keys, values) if R.startswith('R')}
        mol = struct_generator2.get_mol('no_catalyst', Rs, data['stationary_point'])

        with open(runinfoout, 'w+') as out:
            for key, value in zip(keys, values):
                out.write(f'{key}={value}\n')

            if not 'plane_idx' in keys:
                out.write(f'plane_idx={"_".join([str(i) for i in mol.plane_idx])}\n')
            if not 'align_idx' in keys:
                out.write(f'align_idx={"_".join([str(i) for i in mol.align_idx])}\n')
            if not 'center_idx' in keys:
                out.write(f'center_idx={mol.center_idx}\n')
            try:
                if not 'active_atom_idx' in keys:
                    out.write(f'active_atom_idx={mol.active_atom_idx}\n')
            except: pass


    os.replace(runinfoout, runinfo)
import struct_generator2, job_results3
import os, shutil

for d in job_results3.get_all_run_dirs():
    runinfo = os.path.join(d, 'run.info')
    runinfoout = os.path.join(d, 'run.info.test')
    try:
        mol = struct_generator2.get_mol_from_runinfo(runinfo)
    except: continue
    with open(runinfo, 'r') as info:
        info_lines = info.readlines()
        if any(l.startswith('frag1idx') for l in info_lines): continue
        with open(runinfoout, 'w+') as infoout:
            [infoout.write(l) for l in info_lines]
            # if hasattr(mol, 'active_atom_idx'):
            #     infoout.write(f'active_atom_idx={mol.active_atom_idx}\n')
            # if hasattr(mol, 'plane_idx'):
            #     infoout.write(f'plane_idx={"_".join([str(i) for i in mol.plane_idx])}\n')
            # if hasattr(mol, 'align_idx'):
            #     infoout.write(f'align_idx={"_".join([str(i) for i in mol.align_idx])}\n')
            # if hasattr(mol, 'center_idx'):
            #     infoout.write(f'center_idx={mol.center_idx}\n')
            if hasattr(mol, 'frag1idx'):
                infoout.write(f'frag1idx={"_".join([str(i) for i in mol.frag1idx])}\n')
            if hasattr(mol, 'frag2idx'):
                infoout.write(f'frag2idx={"_".join([str(i) for i in mol.frag2idx])}\n')

    os.replace(runinfoout, runinfo)
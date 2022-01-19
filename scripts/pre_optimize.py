import scm.plams as plams
import os 


def pre_optimize(mol, path):
    '''
    Uses DFTB to preoptimize a molecule in xyz file
    '''
    # mol = plams.Molecule(file)
    # name = os.path.basename(file.strip('.xyz'))

    name = mol.name

    settings = plams.Settings()
    settings.input.ams.task = 'GeometryOptimization'
    settings.input.DFTB
    settings.input.DFTB.model = 'GFN1-xTB'

    job = plams.AMSJob(molecule=plams.Molecule(mol.path), settings=settings, name=name + '_preopt')
    res = job.run()

    if check_termination_succes(res):
        molecule = res.get_main_molecule()
        molecule.write(path)

        return molecule
    return mol


def check_termination_succes(result :Results) -> bool:
    term = result.readrkf('General', 'termination status', 'ams')
    if term == 'NORMAL TERMINATION':
        return True
    elif 'NORMAL TERMINATION' in term:
        return 'WARNING'
    return False
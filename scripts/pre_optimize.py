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

    job = plams.AMSJob(molecule=mol, settings=settings, name=name)
    res = job.run()

    if self.check_termination_succes(res):
    molecule = res.get_main_molecule()
    molecule.write(path)

    return molecule
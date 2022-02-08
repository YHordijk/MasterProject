import scm.plams as plams
import struct_generator, paths, pre_optimize, os, multiprocessing, utility
import results_database2 as results_database


class ReactionRunner:
    ''' 
    This class holds functions and methods to generate and run AMSJobs for a given reaction and R-groups
    - First generates coordinates using struct_generator.generate_stationary_points
    - Pre-optimizes them using DFTB optionally (not sure if spinpolarization is possible with dftb)
    - 
    '''
    def __init__(self, template, substituents, pre_optimize=False):
        self.template = template
        self.substituents = substituents

        results_database.update_id_list()
        self.begin_calculations(pre_optimize=pre_optimize)

    def begin_calculations(self, pre_optimize):
        mols = struct_generator.generate_stationary_points(self.template, self.substituents)


        sorted_Rnames = list(sorted(self.substituents.keys()))
        sorted_R = [self.substituents[R] for R in sorted_Rnames]
        reaction_name = f'{self.template}_{"_".join(sorted_R)}'

        plams.init(path=paths.calculations, folder=reaction_name)
        plams.config.default_jobrunner = plams.JobRunner(parallel=True, maxjobs=2)

        workdir = plams.config.default_jobmanager.workdir
        with open(os.path.join(workdir, 'meta.info'), 'w+') as meta:
            meta.write(f'reaction={self.template}\n')
            for R, s in self.substituents.items():
                meta.write(f'{R}={s}\n')

        #here we make sure PLAMS does parallelization
        

        print(f'=== Beginning Reaction {reaction_name}')
        print(f'== Found {len(mols)} molecules to calculate:')
        for m in mols:
            print(f'\t{m.name}: \tpath={os.path.relpath(m.path, paths.master)}, \tflags={m.flags}')

        print('== Checking for hash collisions ...')
        collisions = []
        for m in mols:
            hash = utility.hash(self.template, m.name, m.flags)
            if results_database.get_hash_collision(hash):
                print(f'\t{m.name}\t{hash}')
                collisions.append(m.name)

        
        if pre_optimize: 
            print('== Pre-optimizing molecules ...')
            mols = [pre_optimize.pre_optimize(m, m.path) for m in mols]
        else:
            print('== Skipping pre-optimization step')

        ids = results_database.get_n_free_ids(len(mols))
        dft_jobs = []
        for i, mol in zip(ids, mols):
            if not mol.name in collisions:
                m = utility.load_mol(mol.path)
                s = self.define_settings(m)

                dft_jobs.append(plams.AMSJob(name=f'{i}.{self.template}.{"_".join(sorted_R)}.{mol.name}', settings=s, molecule=m))

        
        #starts runs and wait for finish
        dft_results = [j.run() for j in dft_jobs]
        dft_results = [r for r in dft_results if r.ok()]
        handled_res = []
        while not all([r.ok() for r in dft_results]):
            for res in dft_results:
                if res.ok() and not res in handled_res:
                    res.dill()
                    opt_mol = res.get_main_molecule()
                    opt_mol.flags = res.job.molecule.flags
                    utility.write_mol(opt_mol, )

        print('== Calculations finished')

        plams.finish()

    def define_settings(self, mol):
        s = plams.Settings()
        if 'GO' in mol.flags:
            s.input.ams.task = 'GeometryOptimization'
        elif 'TSRC' in mol.flags:
            s.input.ams.task = 'TransitionStateSearch'
            tsrc = [f for f in mol.flags if f.startswith('TSRC=')][0].split('=')[1].split('_')
            s.input.ams.TransitionStateSearch.ReactionCoordinate.Distance = f'{tsrc[0]} {tsrc[1]} -1.0'

        if 'radical' in mol.flags:
            s.input.adf.Unrestricted = 'Yes'
            s.input.adf.SpinPolarization = '1.0'

        s.input.adf.basis.core = 'None'
        s.input.adf.XC.GGA = 'OLYP'
        s.input.adf.basis.type = 'TZ2P'
        s.input.adf.NumericalQuality = 'VeryGood'
        s.input.ams.Properties.NormalModes = 'Yes'
        s.input.ams.NormalModes.ReScanFreqRange = '-10000000.0 20.0'
        s.input.adf.SYMMETRY = 'NOSYM'
        s.input.adf.SCF.Iterations = '99'
        s.input.adf.SCF.Converge = '1.0e-6'
        return s


if __name__ == '__main__':
    r = ReactionRunner('urea_tBu_Ph', {'R1':'H', 'R2':'Ph', 'Rch':'O'})

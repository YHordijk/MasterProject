import scm.plams as plams
import struct_generator, paths, results_database, pre_optimize, os, multiprocessing, utility



class ReactionRunner:
    ''' 
    This class holds functions and methods to generate and run AMSJobs for a given reaction and R-groups
    - First generates coordinates
    - Pre-optimizes them
    - 
    '''
    def __init__(self, template, substituents):
        self.template = template
        self.substituents = substituents
        self.begin_calculations()

    def begin_calculations(self):
        plams.init()

        #here we make sure PLAMS does parallelization
        plams.config.default_jobrunner = plams.JobRunner(parallel=True, maxjobs=multiprocessing.cpu_count())
        plams.config.job.runscript.nproc = 1

        sorted_Rnames = list(sorted(self.substituents.keys()))
        sorted_R = [self.substituents[R] for R in sorted_Rnames]
        reaction_name = f'{self.template}_{"_".join(sorted_R)}'

        print(f'Beginning Reaction {reaction_name}')
        mols = struct_generator.generate_stationary_points(self.template, self.substituents)
        print(f'Found {len(mols)} molecules to calculate:')
        for m in mols:
            print(f'\t{m.name}: \tpath={os.path.relpath(m.path, paths.master)}, \tflags={m.flags}')

        print('Pre-optimizing molecules ...')
        # mols = [pre_optimize.pre_optimize(m, m.path) for m in mols]

        dft_jobs = []
        for mol in mols:
            m = utility.load_mol(mol.path)
            s = self.define_settings(m)
            dft_jobs.append(plams.AMSJob(name=mol.name, settings=s, molecule=m))

        #starts runs and wait for finish
        dft_results = [j.run() for j in dft_jobs]
        dft_results = [r for r in dft_results if r.ok()]

        print(dft_results)



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
        print(s)
        return s


if __name__ == '__main__':
    r = ReactionRunner('no_catalyst', {'R1':'H', 'R2':'H'})

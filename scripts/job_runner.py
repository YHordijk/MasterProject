import scm.plams as plams
import struct_generator, paths, results_database, pre_optimize, os



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

        sorted_Rnames = list(sorted(self.substituents.keys()))
        sorted_R = [self.substituents[R] for R in sorted_Rnames]
        reaction_name = f'{self.template}_{"_".join(sorted_R)}'

        print(f'Beginning Reaction {reaction_name}')
        self.init_mols = struct_generator.generate_stationary_points(self.template, self.substituents)
        print(f'Found {len(self.init_mols)} molecules to calculate:')
        for m in self.init_mols:
            print(f'\t{m.name}: path={os.path.relpath(m.path, paths.master)}, flags=({m.flags})')

        print('Pre-optimizing molecules ...')
        self.preopt_mols = [pre_optimize.pre_optimize(m, m.path) for m in self.init_mols]




        plams.finish()

    def define_settings(self, mol):
        s = plams.Settings()
        if 'GO' in mol.flags:
            s.input.ams.task = 'GeometryOptimization'
        elif 'TSRC' in mol.flags:
            s.input.ams.task = 'TransitionStateSearch'

        if 'radical' in mol.flags:
            s.input.adf.Unrestricted = 'Yes'
            s.input.adf.SpinPolarization = '1.0'

        s.input.adf.basis.core = 'None'
        s.input.adf.xc = 'OLYP'
        s.input.adf.basis.type = 'TZ2P'
        s.input.adf.NumericalQuality = 'Very Good'


if __name__ == '__main__':
    r = ReactionRunner('no_catalyst', {'R1':'H', 'R2':'H'})

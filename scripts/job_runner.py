import scm.plams as plams
import struct_generator, paths, results_database, pre_optimize



def ReactionRunner:
    ''' 
    This class holds functions and methods to generate and run AMSJobs for a given reaction and R-groups
    - First generates coordinates
    - Pre-optimizes them
    - 
    '''
    def __init__(self, template, substituents):
        sorted_Rnames = list(sorted(substituents.keys()))
        sorted_R = [substituents[R] for R in sorted_Rnames]
        reaction_name = f'{template}_{"_".join(sorted_R)}'

        print(f'Beginning Reaction {reaction_name}')
        plams.init()
        self.template = template
        self.substituents = substituents

        self.init_mols = struct_generator.generate_stationary_points(template, substituents)
        print(f'Found {len(self.init_mols)} molecules to calculate:')
        for m in self.init_mols:
            print(f'\t{m.name}: {m.path}')

        print('Pre-optimizing molecules ...')
        self.preopt_mols = [pre_optimize.pre_optimize(m, m.path) for m in self.init_mols]
        plams.finish()


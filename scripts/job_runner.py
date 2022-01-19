import scm.plams as plams
import struct_generator, paths, results_database, pre_optimize



class ReactionRunner:
    ''' 
    This class holds functions and methods to generate and run AMSJobs for a given reaction and R-groups
    - First generates coordinates
    - Pre-optimizes them
    - 
    '''
    def __init__(self, template, R):
        plams.init()
        self.template = template
        self.R = R

        self.init_mols = struct_generator.generate_stationary_points(template, R)
        print(self.init_mols)
        self.preopt_mols = [pre_optimize.pre_optimize(m, m.path) for m in self.init_mols]
        plams.finish()


if __name__ == '__main__':
    r = ReactionRunner('no_catalyst', {'R1':'H', 'R2':'H'})

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
		self.template = template
		self.R = R

		self.init_mols = struct_generator.generate_stationary_points(template, R)
		self.preopt_mols = [pre_optimize.pre_optimize(m) for m in self.init_mols]


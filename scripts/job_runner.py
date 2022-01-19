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
		


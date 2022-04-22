import scm.plams as plams
import workflow_module.results as wfr
import paths, struct_generator2, os, utility
join = os.path.join



def sort_substituents(subs):
	sorted_Rnames = list(sorted(subs))
	sorted_R = [subs[R] for R in sorted_Rnames]
	return sorted_R



class Task:
	def __init__(self, mol, settings, task_name, **kwargs):
		self.mol = mol
		self.settings = settings
		self.task_name = task_name
		self.kwargs = kwargs

		print(f'Job (task={task_name}, mol={kwargs["stationary point"]})')


	def write_info(self):
		# WRITE ALSO SOME INFO
		with open(os.path.join(self.workdir, 'job.info'), 'w+') as ji:
			ji.write(f'reaction={self.kwargs["reaction"]}\n')
			ji.write(f'stationary_point={self.kwargs["stationary point"]}\n')
			ji.write(f'task={self.task_name}\n')
			ji.write(f'TSRC_idx={"_".join([str(i) for i in self.kwargs.get("TSRC idx")])}\n')
			ji.write(f'radical={self.kwargs["radical"]}\n')
			ji.write(f'phase={self.kwargs["solvent"]}\n')
			ji.write(f'enantiomer={self.kwargs["enantiomer"]}\n')
			sorted_substituents = "_".join(sort_substituents(self.kwargs['substituents']))
			ji.write(f'sorted_substituents={sorted_substituents}\n')
			for R, sub in self.kwargs['substituents'].items():
				ji.write(f'{R}={sub}\n')
			ji.write(f'functional={self.kwargs["functional"]}\n')
			ji.write(f'basis={self.kwargs["basis"]}\n')
			ji.write(f'numerical_quality={self.kwargs["quality"]}\n')
			ji.write(f'frozen_core={self.kwargs["frozen core"]}\n')
			ji.write(f'plane_idx={"_".join([str(i) for i in self.kwargs["plane idx"]])}\n')
			ji.write(f'align_idx={"_".join([str(i) for i in self.kwargs["align idx"]])}\n')
			ji.write(f'center_idx={self.kwargs["center idx"]}\n')
			ji.write(f'active_atom_idx={self.kwargs["active atom idx"]}\n')


	def postrun(self):
		result = wfr.generate_result(self.workdir)
		with open(os.path.join(self.workdir, 'SUCCESS'), 'w+') as _: pass
		return result


	def write_hash(self, job):
		job_hash = job.hash_input()
		with open(os.path.join(self.workdir, 'job.hash'), 'w+') as jh:
			jh.write(job_hash)


	def write_mol(self, mol):
		mol.write(os.path.join(self.workdir, 'input.xyz'))


	def write_input(self, job):
		with open(os.path.join(self.workdir, 'job.input'), 'w+') as ji:
			ji.write(job.get_input())


	def skip_job(self, job):
		job_hash = job.hash_input()
		p = self.kwargs['calc dir']
		for a in os.listdir(p):
			pa = os.path.join(p, a)
			# print('pa',pa)
			for b in os.listdir(pa):
				pb = os.path.join(pa, b)
				# print('pb',pb)
				if os.path.isdir(pb):
					for c in os.listdir(pb):
						pc = os.path.join(pb, c)
						# print('pc',pc)
						# print(os.path.exists(os.path.join(pc, 'job.hash')), os.path.join(pc, 'job.hash'))
						# print(os.path.exists(os.path.join(pc, 'SUCCESS')), os.path.join(pc, 'SUCCESS'))
						if os.path.exists(os.path.join(pc, 'job.hash')):
							with open(os.path.join(pc, 'job.hash')) as pd:
								 if pd.read().strip() == job_hash:
								 	if os.path.exists(os.path.join(pc, 'SUCCESS')):
									 	return pc
		return False



class Optimization(Task):
	def run(self):
		### CREATE JOB AND CHECK IF IT HAS ALREADY BEEN RUN BEFORE
		job = plams.AMSJob(molecule=self.mol, settings=self.settings, name='ADF')

		skip_dir = self.skip_job(job)
		if not self.kwargs['rerun calculations'] and skip_dir != False:
			print(f'Skipping, already found the job in {skip_dir}')
			return wfr.load_result(skip_dir)

		print('Running ...')

		### MAKE SOME PATHS
		job_dir = os.path.join(self.kwargs['run dir'], self.kwargs['stationary point'])
		os.makedirs(job_dir, exist_ok=True)

		plams.init(folder=self.task_name, path=job_dir)
		self.workdir = plams.config.default_jobmanager.workdir

		self.write_info()
		self.write_hash(job)
		self.write_mol(self.mol)
		self.write_input(job)


		### RUN THE JOB
		res = job.run()
		result = self.postrun(job, res)

		plams.finish()

		return result
		# return self.mol


	def postrun(self, job, res):
	    resfile = plams.KFFile(res['adf.rkf'])
	    cosmo_data = resfile.read_section('COSMO')
	    coskf = plams.KFFile(join(job.path, 'solute.coskf'))
	    for k, v in cosmo_data.items():
	        coskf.write('COSMO', k, v)
	    res.collect()

	    result = wfr.generate_result(self.workdir)
	    with open(os.path.join(self.workdir, 'SUCCESS'), 'w+') as _: pass
	    return result




class COSMORS(Task):
	def __init__(self, settings, task_name, **kwargs):
		self.settings = settings
		self.task_name = task_name
		self.kwargs = kwargs

		print(f'Job (task={task_name}, mol={kwargs["stationary point"]})')

	def run(self):
		job_dir = os.path.join(self.kwargs['run dir'], self.kwargs['stationary point'])
		os.makedirs(job_dir, exist_ok=True)

		print('Running ...')

		plams.init(folder=self.task_name, path=job_dir)
		self.workdir = plams.config.default_jobmanager.workdir

		job = plams.CRSJob(settings=self.settings, name=self.task_name)

		self.write_info()
		self.write_input(job)

		res = job.run().get_results()
		result = self.postrun()

		plams.finish()

		return result

	

class Frequencies(Task):
	def run(self):
		### CREATE JOB AND CHECK IF IT HAS ALREADY BEEN RUN BEFORE
		job = plams.AMSJob(molecule=self.mol, settings=self.settings, name='ADF')

		skip_dir = self.skip_job(job)
		if not self.kwargs['rerun calculations'] and skip_dir != False:
			print(f'Skipping, already found the job in {skip_dir}')
			return wfr.load_result(skip_dir)

		print('Running ...')
		### MAKE SOME PATHS
		job_dir = os.path.join(self.kwargs['run dir'], self.kwargs['stationary point'])
		os.makedirs(job_dir, exist_ok=True)

		plams.init(folder=self.task_name, path=job_dir)
		self.workdir = plams.config.default_jobmanager.workdir

		### BEFORE RUNNING THE JOB WRITE SOME FILES TO RUNDIR
		# WRITE HASH SO WE CAN READ IT IN LATER
		self.write_info()
		self.write_hash(job)
		self.write_input(job)

		### RUN THE JOB
		res = job.run()
		result = self.postrun()

		plams.finish()

		return result
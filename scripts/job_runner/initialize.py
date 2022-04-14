import struct_generator2, paths, os, utility

join = os.path.join


def initialize(template, substituents={}, calc_dir=paths.calculations, test_mode=False):
	def get_mols():
		print(f'Generating molecules:')
		print(f'\tTemplate  : {template}')
		for name, sub in substituents.items():
			print(f'\t{name:<8}  : {sub}')
		mols = struct_generator2.generate_stationary_points(template, substituents)
		print(f'\nGenerated {len(mols)} molecules:')
		struct_generator2.print_mols(mols, tabs=1)
		print()
		return mols

	def sort_substituents(subs):
		sorted_Rnames = list(sorted(subs))
		sorted_R = [substituents[R] for R in sorted_Rnames]
		return sorted_R

	def get_job_dir(mol):
		sorted_R = sort_substituents(substituents)
		return join(calc_dir, template + '.' + '_'.join(sorted_R), mol.name)

	def meta_info():
		sorted_R = sort_substituents(substituents)
		with open(join(calc_dir, template + '.' + '_'.join(sorted_R), 'reaction.info'), 'w+') as meta:
			meta.write(f'reaction={template}\n')
			meta.write(f'sorted_substituents={"_".join(sort_substituents(substituents))}\n')
			for sub in substituents:
				meta.write(f'{sub}={substituents[sub]}\n')

	def write_info(mol):
		with open(join(get_job_dir(mol), 'mol.info'), 'w+') as info:
			info.write(f'reaction={template}\n')
			info.write(f'stationary_point={mol.name}\n')
			info.write(f'TSRC_idx={"_".join(str(i) for i in mol.TSRC_idx)}\n' if len(mol.TSRC_idx) > 0 else 'TSRC_idx=N/A\n')
			info.write(f'radical={mol.radical}\n')
			info.write(f'sorted_substituents={"_".join(sort_substituents(mol.substituents))}\n')
			for sub in mol.substituents:
				info.write(f'{sub}={substituents[sub]}\n')
			info.write(f'enantiomer={mol.enantiomer}\n')
			if hasattr(mol, 'active_atom_idx'):
				info.write(f'active_atom_idx={mol.active_atom_idx}\n')
			if hasattr(mol, 'plane_idx'):
				info.write(f'plane_idx={"_".join([str(i) for i in mol.plane_idx])}\n')
			if hasattr(mol, 'align_idx'):
				info.write(f'align_idx={"_".join([str(i) for i in mol.align_idx])}\n')
			if hasattr(mol, 'center_idx'):
				info.write(f'center_idx={mol.center_idx}\n')
			if hasattr(mol, 'frag1idx'):
				info.write(f'frag1idx={"_".join([str(i) for i in mol.frag1idx])}\n')
			if hasattr(mol, 'frag2idx'):
				info.write(f'frag2idx={"_".join([str(i) for i in mol.frag2idx])}\n')

	mols = get_mols()
	for name, mol in mols.items():
		mol_job_dir = get_job_dir(mol)
		os.makedirs(mol_job_dir, exist_ok=True)
		write_info(mol)
		utility.write_mol(mol, join(mol_job_dir, 'init.xyz'))

	meta_info()

	return join(calc_dir, template + '.' + '_'.join(sort_substituents(substituents)))

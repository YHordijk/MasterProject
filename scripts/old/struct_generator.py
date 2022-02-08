import scm.plams as plams
import os, paths, utility


def sub_path(name):
    return os.path.join(paths.SGT_substituents, name + '.sub')


def reaction_path(name):
    return os.path.join(paths.SGT, name + '.tmplt')


def generate_stationary_points(template_name, substituents={}, keep_dummy=False):
    '''
    This function generates stationary points according to a given template
    It will split the file up into the different molecules
    It will then find substituents at each molecule and try to substitute them

    template_name: name of template file
    substituent: dictionary with substituent group as key and substituent name as value
    '''

    def parse_contents(lines, issub=False):
        lines = [l.strip() for l in lines]

        #read from file
        name = lines[0]
        flags = lines[1].split()
        struct = [l.split() for l in lines[2:]]
        elements = [s[0] for s in struct]
        coords = [[float(x) for x in s[1:4]] for s in struct]
        tags = [s[4:] for s in struct]
        #create molecule
        mol = plams.Molecule()
        mol.name = name
        mol.reaction = template_name
        mol.flags = flags
        atoms = [plams.Atom(symbol=s, coords=c) for s, c in zip(elements, coords)]
        
        if issub: mol.connector = [None, None]
        else: 
            mol.connector = {}
            mol.atoms_to_delete = []
            TSRC = []

        for i, a, t in zip(range(len(atoms)), atoms, tags):
            if not issub:
                for f in t:
                    if f.startswith('R'):
                        name = f[0:-1]
                        if not name in mol.connector: mol.connector[name] = [None, None]
                        if f[-1] == 'a': mol.connector[name][0] = a
                        elif f[-1] == 'b': mol.connector[name][1] = a

                    elif f == 'TSRC':
                        TSRC.append(a)

                    elif f == 'delete':
                        mol.atoms_to_delete.append(a)

                    elif f.startswith('dist='):
                        #find Rgroup
                        for x in t:
                            if f.startswith('R'):
                                name = f[0:-1]
                        mol.connector_distance = (name, float(f.split('=')[-1]))

            if issub:
                if   'a' in t: mol.connector[0] = i
                elif 'b' in t: mol.connector[1] = i
        
        [mol.add_atom(a) for a in atoms]

        if issub: mol.connector = tuple(mol.connector)
        else: 
            mol.connector = {R: tuple(c) for R, c in mol.connector.items()}
            if len(TSRC) > 0:
                TSRCidx = (mol.atoms.index(TSRC[0])+1, mol.atoms.index(TSRC[1])+1)
                mol.flags.append(f'TSRC={TSRCidx[0]}_{TSRCidx[1]}')

        return mol

    #First parse substituents
    sub_mols = {}
    for R, p in substituents.items():
        with open(sub_path(p), 'r') as file:
            sub_mols[R] = parse_contents(file.readlines(), issub=True)
    with open(sub_path('H'), 'r') as file:
        default_sub = parse_contents(file.readlines(), issub=True)

    all_substituent_names = {}
    # if R 
    # substituents[R] = sub_mols[R].name if R in substituents else default_sub.name
    with open(reaction_path(template_name), 'r') as reaction:
        content = reaction.read().split('\n\n')
        mols = [parse_contents(c.split('\n')) for c in content]
        for m in mols:
            m.substituents = {}
            for R in m.connector.keys():
                if R in substituents:
                    m.substituents[R] = sub_mols[R].name
                elif not keep_dummy: 
                    default_sub.name

                if not R in all_substituent_names:
                    all_substituent_names[R] = sub_mols[R].name if R in substituents else default_sub.name

            for R in m.connector:
                use_default = not R in sub_mols
                if use_default:   s = default_sub.copy()
                else:             s = sub_mols[R].copy()

                if not (use_default and keep_dummy):
                    lconn = (s.atoms[s.connector[0]], s.atoms[s.connector[1]])

                    bond_length = None
                    if hasattr(m, 'connector_distance'):
                        if m.connector_distance[0] == R:
                            bond_length = m.connector_distance[1]
                        
                    m.substitute(m.connector[R], s, lconn, bond_length=bond_length)
                    for dela in m.atoms_to_delete:
                        m.delete_atom(dela)

        sorted_Rnames = list(sorted(all_substituent_names.keys()))
        sorted_R = [all_substituent_names[R] for R in sorted_Rnames]
        for m in mols:
            mpath = get_mol_path(paths.input_xyz, template_name + '_' + '_'.join(sorted_R), f'{m.name}')
            m.path = mpath
            os.makedirs(os.path.dirname(mpath), exist_ok=True)
            [m.flags.append(r) for r in [R+"="+n for R, n in m.substituents.items()]]
            comment = ", ".join(m.flags)
            utility.write_mol(m, m.path, comment=comment)

    return mols


def get_mol_path(base, dir, name):
    return os.path.join(base, dir, name + '.xyz')


def get_hashes_for_calc(template, substituents):
    mols = generate_stationary_points(template, substituents)
    hashes = {}
    for mol in mols:
        hashes[mol.name] = utility.hash(template, mol.name, mol.flags)
    return hashes


def get_subs_used_for_sp(template, sp):
    with open(reaction_path(template)) as tmplt:
        text = tmplt.read()
        mols = text.split('\n\n')
        sp_names = [s.split()[0] for s in mols]
        mol = mols[sp_names.index(sp)]
        substituents = set()
        for line in mol.split('\n')[2:]:
            tags = line.split()[4:]
            for tag in tags:
                if tag.startswith('R'):
                    substituents.add(tag[:-1])
    return substituents

def get_flags_for_sp(template, sp):
    subs = get_subs_used_for_sp(template, sp)
    with open(reaction_path(template)) as tmplt:
        text = tmplt.read()
        mols = text.split('\n\n')
        sp_names = [s.split()[0] for s in mols]
        mol = mols[sp_names.index(sp)]
        flags = mol.split('\n')[1].split()
        return flags


def get_hash_for_sp(template, substituents, sp):
    return get_hashes_for_calc(template, substituents)[sp]

# def _get_flags(results):
#     sorted_Rnames = list(sorted(results['all_substituents'].keys()))
#     sorted_R = [results['all_substituents'][R] for R in sorted_Rnames]
#     input_file = struct_generator.get_mol_path(paths.input_xyz, results['reaction']+'_'+'_'.join(sorted_R), results['stationary_point'])
#     if not os.path.exists(input_file):
#         results['flags'] = []
#     with open(input_file) as file:
#         lines = file.readlines()
#         results['flags'] = [f.strip() for f in lines[1].split(',')]

#     results['substituents'] = {}
#     for f in results['flags']:
#         if f.startswith('R'):
#             results['substituents'][f.split('=')[0]] = f.split('=')[1].strip()
#             results[f.split('=')[0]] = f.split('=')[1].strip()

#     results['radical'] = 'radical' in results['flags']
#     results['COSMO'] = 'COSMO' in results['flags']

#     for task in ['GO', 'SP', 'TSRC', 'LT']:
#         if task in results['flags']: 
#             results['task'] = task




if __name__ == '__main__':
    # for Rcat in ['TiCl4', 'SnCl4', 'I2', 'ZnCl2', 'AlF3', 'BF3']:
    #     generate_stationary_points('achiral_catalyst', {'R2':'Ph', 'R1':'H', 'Rcat':Rcat})

    generate_stationary_points('urea_tBu_Ph', {'R2':'H', 'R1':'H'})
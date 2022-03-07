import scm.plams as plams
import os, paths, utility
import matplotlib.pyplot as plt
try:
    import mol_viewer2
except: pass
join = os.path.join


def get_mol(template, substituents, stationary_point):
    mols = generate_stationary_points(template, substituents)
    if stationary_point in mols:
        return mols[stationary_point]
    else:
        return None


def get_all_substituents(template):
    template_dir = join(paths.SGT, template)
    template_meta = join(template_dir, 'meta.info')

    assert os.path.isdir(template_dir), f'Directory {template_dir} does not exist'
    assert os.path.exists(template_meta), f'No meta.info file found in {template_dir}'

    with open(template_meta) as meta:
        subs = []
        for line in meta.readlines():
            if line.startswith('subdefault'):
                subs.append(line.split()[1])
    return subs


def get_all_stationary_points(template):
    template_dir = join(paths.SGT, template)
    assert os.path.isdir(template_dir), f'Directory {template_dir} does not exist'

    sps = set()
    for file in os.listdir(template_dir):
        if file.endswith('.xyz'):
            sps.add(file[:-4])
    return sps


def generate_stationary_points(template, substituents=None, keep_dummy=False):
    template_dir = join(paths.SGT, template)
    template_meta = join(template_dir, 'meta.info')

    assert os.path.isdir(template_dir), f'Directory {template_dir} does not exist'
    assert os.path.exists(template_meta), f'No meta.info file found in {template_dir}'

    def load_mol(file):
        name = os.path.basename(file).split('.')[0]
        with open(file, 'r') as f:
            mol = plams.Molecule(file)

            lines = [l.strip() for l in f.readlines()]

            flags = lines[1].split()
            #parse flags
            task = 'GO'
            radical = False
            get_frags = False
            enant = 'N/A'
            TSRC_idx = []
            substituent_idx = {}
            substituent_dist = {}
            delete_idx = []
            active_atom_idx = False
            for flag in flags:
                if flag in ['GO']:
                    task = flag
                if flag == 'radical':
                    radical = True
                if flag.startswith('enant='):
                    enant = flag.split('=')[1]
                if flag.startswith('R'):
                    n, idx = flag.split('=')
                    if len(idx.split('_')) == 2:
                        substituent_idx[n] = [int(i) for i in idx.split('_')]
                        substituent_dist[n] = None
                    elif len(idx.split('_')) == 3:
                        substituent_idx[n] = [int(i) for i in idx.split('_')[0:2]]
                        substituent_dist[n] = float(idx.split('_')[2])
                if flag.startswith('TSRC'):
                    task = 'TSRC'
                    idx = flag.split('=')[1]
                    TSRC_idx = [int(i) for i in idx.split('_')]
                if flag.startswith('delete'):
                    delete_idx = [int(i) for i in flag.split('=')[1].split('_')]
                if flag == 'FRAG':
                    get_frags = True
                if flag.startswith('active_atom='):
                    active_atom_idx = flag.split('=')[1]
                if flag.startswith('plane='):
                    plane = [int(i) for i in flag.split('=')[1].split('_')]
                    mol.plane_idx = [mol.atoms[i-1] for i in plane]
                if flag.startswith('align='):
                    align = [int(i) for i in flag.split('=')[1].split('_')]
                    mol.align_idx = [mol.atoms[i-1] for i in align]
                if flag.startswith('center='):
                    center = int(flag.split('=')[1])
                    mol.center_idx = mol.atoms[center-1]

            #construct our molecule object
            mol.name = name
            mol.reaction = template
            mol.task = task
            mol.radical = radical
            mol.enantiomer = enant
            mol.TSRC_idx = TSRC_idx
            mol.substituent_idx = substituent_idx
            mol.substituents = substituent_idx.keys()
            mol.substituent_dist = substituent_dist
            mol.delete_atoms = [mol.atoms[i-1] for i in delete_idx]
            mol.get_frags = get_frags
            mol.template_mol = mol
            if active_atom_idx:
                mol.active_atom = mol.atoms[int(active_atom_idx)]
            
        return mol


    def get_all_subs(mols):
        all_substituents = set()
        for mol in template_mols:
            for sub in mol.substituent_idx.keys():
                all_substituents.add(sub)
        return list(all_substituents)


    def get_default_subs(all_substituents):
        default = {name: 'H' for name in all_substituents}
        with open(template_meta) as meta:
            for line in meta.readlines():
                if line.startswith('subdefault'):
                    name, sub = line.split()[1:3]
                    default[name] = sub
        return default


    def get_sub_files(subs):
        substituent_files = {}
        for name, sub in subs.items():
            substituent_files[name] = join(paths.SGT_substituents, sub + '.xyz')
            assert os.path.exists(substituent_files[name]), f'File {substituent_files[name]} does not exist'
        return substituent_files


    def substitute_mol(mol):
        def load_substituent(file):
            substituent_idx = []
            with open(file) as f:
                for i, line in enumerate(f.readlines()[2:], start=1):
                    #sub indices are tagged with R
                    tag = line.split()[4:]
                    if 'R' in tag:
                        substituent_idx.append(i)
            mol = plams.Molecule(file)
            mol.substituent_idx = substituent_idx
            return mol

        def get_connector(mol, idx):
            atoms = mol.atoms[idx[0]-1], mol.atoms[idx[1]-1]
            elements = [atom.symbol for atom in atoms] #sort based on element
            #Xx element must always be last
            if elements[1] == 'Xx':
                return atoms
            else:
                return atoms[::-1]


        #save TSRC atoms so we can update them later as new idx
        TSRC_atoms = [mol.atoms[i-1] for i in mol.TSRC_idx]

        #first get sub_atoms as the order of atoms can change during substitutions
        main_conns = {}
        for sub_name, idx in mol.substituent_idx.items():
            main_conns[sub_name] = get_connector(mol, idx) #indices start at 1

        #now do the substitutions
        for sub_name, atoms in main_conns.items():
            if not sub_name in substituent_files:
                continue
            sub_file = substituent_files[sub_name]
            #load the substituent molecule from file
            #also get the connector atoms
            sub_mol = load_substituent(sub_file)
            sub_conn = get_connector(sub_mol, sub_mol.substituent_idx)

            #substitute the molecule
            main_conn = main_conns[sub_name]
            # print(sub_name, template, substituents, mol.name)
            mol.substitute(main_conn, sub_mol, sub_conn, bond_length=mol.substituent_dist[sub_name], steps=1)

        [mol.delete_atom(a) for a in mol.delete_atoms]

        #atom indices have changed so update the TSRC indices
        mol.TSRC_idx = [mol.atoms.index(a) + 1 for a in TSRC_atoms]

    def get_fragments(mol):
        #we must copy the atoms into a new molecule as 
        #the properties we glued on mol are not picklable

        newm = plams.Molecule()
        newm.add_molecule(mol)
        newm.delete_all_bonds()
        newm.guess_bonds()
        frag1, frag2 = newm.separate()
        frag1coord = [a.coords for a in frag1.atoms]
        frag2coord = [a.coords for a in frag2.atoms]
        frag1idx = []
        frag2idx = []
        for i, a in enumerate(newm.atoms, 1):
            if a.coords in frag1coord:
                frag1idx.append(i)
            elif a.coords in frag2coord:
                frag2idx.append(i)
            else: raise
        mol.frag1idx = frag1idx
        mol.frag2idx = frag2idx

    def set_active_atom_idx(mol):
        if hasattr(mol, 'active_atom'):
            i = mol.atoms.index(mol.active_atom)
            mol.active_atom_idx = i

    def set_plane_idx(mol):
        if hasattr(mol, 'plane_idx'):
            idx = [mol.atoms.index(i) for i in mol.plane_idx]
            mol.plane_idx = idx

    def set_align_idx(mol):
        if hasattr(mol, 'align_idx'):
            idx = [mol.atoms.index(i) for i in mol.align_idx]
            mol.align_idx = idx

    def set_center_idx(mol):
        if hasattr(mol, 'center_idx'):
            idx = mol.atoms.index(mol.center_idx)
            mol.center_idx = idx

    template_files = [join(template_dir, f) for f in os.listdir(template_dir) if f.endswith('.xyz')]
    #load the molecules here
    template_mols = []
    for file in template_files:
        template_mols.append(load_mol(file))

    #load the substituents
    #get a list of all substituents first
    all_substituents = get_all_subs(template_mols)
    #read in default substituents from meta.info
    if not keep_dummy:
        _substituents = get_default_subs(all_substituents)
        _substituents.update(substituents) #update with provided subs
    else:
        _substituents = substituents
    #finally get the correct files
    substituent_files = get_sub_files(_substituents)

    #read default substituents:
    for mol in template_mols:
        substitute_mol(mol)

    for mol in template_mols:
        try:
            if mol.get_frags:
                get_fragments(mol)
        except:
            pass
        set_active_atom_idx(mol)
        set_plane_idx(mol)
        set_align_idx(mol)
        set_center_idx(mol)

    return {mol.name: mol for mol in template_mols}



def print_mols(mols, tabs=0):
    name_len = max(len(name) for name in mols)
    header = '\t'*tabs + f'Mol   {"Name".center(name_len)}   Task   Radical   Enantiomer   TSRC indices'
    print(header)
    print('\t'*tabs + '-'*len(header))
    i = 0
    for name, mol in mols.items():
        i += 1
        print('\t'*tabs + f'{i:<3} | {name:<{name_len}} | {mol.task.center(4)} | {str(mol.radical)[0].center(7)} | {mol.enantiomer.center(10)} | {" ".join(str(i).center(3) for i in mol.TSRC_idx)}')

def show_reaction(template, substituents=None, simple=False):
    if substituents is None:
        substituents = {}
    mols = generate_stationary_points(template, substituents, keep_dummy=True)
    mol_viewer2.show(list(mols.values()), simple=simple)

if __name__ == '__main__':
    # mols = generate_stationary_points('achiral_catalyst', {'Rcat':'AlF3'})
    # print_mols(mols)
    show_reaction('achiral_catalyst', {'Rcat':'SnCl4', 'R2':'m-FPh', 'R1':'Me'})

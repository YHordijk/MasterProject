import scm.plams as plams
import os, paths, utility

join = os.path.join




def generate_stationary_points(template, substituents=None):
    template_dir = join(paths.SGT, template)
    template_meta = join(template_dir, 'meta.info')

    
    assert os.path.isdir(template_dir), f'Directory {template_dir} does not exist'
    assert os.path.exists(template_meta), f'No meta.info file found in {template_dir}'


    def load_mol(file):
        name = os.path.basename(file).split('.')[0]
        with open(file, 'r') as f:
            lines = [l.strip() for l in f.readlines()]

            flags = lines[1].split()
            #parse flags
            task = 'GO'
            radical = False
            enant = 'N/A'
            for flag in flags:
                if flag in ['GO', 'TSRC']:
                    task = flag
                if flag == 'radical':
                    radical = True
                if flag.startswith('enant='):
                    enant = flag.split('=')[1]

            #read the rest of the lines
            #get TSRC indices and substituent indices
            TSRC_idx = []
            substituent_idx = {}
            for i, line in enumerate(lines[2:], start=1):
                tags = line.split()[4:]
                for tag in tags:
                    if tag.startswith('R'): #all substituent tags start with R
                        if tag in substituent_idx:
                            substituent_idx[tag].append(i)
                        else:
                            substituent_idx[tag] = [i]

                    if tag == 'TSRC':
                        TSRC_idx.append(i)

            #construct our molecule object
            mol = plams.Molecule(file)
            mol.name = name
            mol.reaction = template
            mol.task = task
            mol.radical = radical
            mol.enantiomer = enant
            mol.TSRC_idx = TSRC_idx
            mol.substituent_idx = substituent_idx
            mol.substituents = substituent_idx.keys()
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
            sub_file = substituent_files[sub_name]
            #load the substituent molecule from file
            #also get the connector atoms
            sub_mol = load_substituent(sub_file)
            sub_conn = get_connector(sub_mol, sub_mol.substituent_idx)

            #substitute the molecule
            main_conn = main_conns[sub_name]
            mol.substitute(main_conn, sub_mol, sub_conn)

        #atom indices have changed so update the TSRC indices
        mol.TSRC_idx = [mol.atoms.index(a) + 1 for a in TSRC_atoms]


    template_files = [join(template_dir, f) for f in os.listdir(template_dir) if f.endswith('.xyz')]
    #load the molecules here
    template_mols = []
    for file in template_files:
        template_mols.append(load_mol(file))

    #load the substituents
    #get a list of all substituents first
    all_substituents = get_all_subs(template_mols)
    #read in default substituents from meta.info
    _substituents = get_default_subs(all_substituents)
    _substituents.update(substituents) #update with provided subs
    #finally get the correct files
    substituent_files = get_sub_files(_substituents)

    #read default substituents:
    for mol in template_mols:
        substitute_mol(mol)

    return {mol.name: mol for mol in template_mols}



def print_mols(mols, tabs=0):
    name_len = max(len(name) for name in mols)
    header = '\t'*tabs + f'Mol   {"Name".center(name_len)}   Task   Radical   Enantiomer   TSRC indices'
    print(header)
    print('\t'*tabs + '_'*len(header))
    i = 0
    for name, mol in mols.items():
        i += 1
        print('\t'*tabs + f'{i:<3} | {name:<{name_len}} | {mol.task.center(4)} | {str(mol.radical)[0].center(7)} | {mol.enantiomer.center(10)} | {" ".join(str(i).center(3) for i in mol.TSRC_idx)}')


if __name__ == '__main__':
    mols = generate_stationary_points('urea_tBu_Ph', {'R1':'H', 'Rch': 'S'})
    print_mols(mols)
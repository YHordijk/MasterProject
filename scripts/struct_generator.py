import scm.plams as plams
import os, paths, utility


def sub_path(name):
    return os.path.join(paths.SGT_substituents, name + '.sub')

def reaction_path(name):
    return os.path.join(paths.SGT, name + '.tmplt')


def generate_stationary_points(template_name, substituents={}):
    '''
    This function generates stationary points according to a given template
    It will split the file up into the different molecules
    It will then find substituents at each molecule and try to substitute them

    template_name: name of template file
    substituent: dictionary with substituent group as key and substituent name as value
    '''

    substituents['default'] = 'H'

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
        mol.flags = flags
        atoms = [plams.Atom(symbol=s, coords=c) for s, c in zip(elements, coords)]
        
        if issub: mol.connector = [None, None]
        else: 
            mol.connector = {}
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


            if issub:
                if   'a' in t: mol.connector[0] = i
                elif 'b' in t: mol.connector[1] = i
        

        [mol.add_atom(a) for a in atoms]

        if issub: mol.connector = tuple(mol.connector)
        else: 
            mol.connector = {R: tuple(c) for R, c in mol.connector.items()}
            if len(TSRC) > 0:
                TSRCidx = (mol.atoms.index(TSRC[0]), mol.atoms.index(TSRC[1]))
                mol.flags.append(f'TSRC={TSRCidx[0]}_{TSRCidx[1]}')

        print(mol.flags)

        
        return mol

    #First parse substituents
    sub_mols = {}
    for R, p in substituents.items():
        with open(sub_path(p), 'r') as file:
            sub_mols[R] = parse_contents(file.readlines(), issub=True)

    with open(reaction_path(template_name), 'r') as reaction:
        content = reaction.read().split('\n\n')
        mols = [parse_contents(c.split('\n')) for c in content]
        for m in mols:
            m.substituents = {R: sub_mols[R].name for R in m.connector.keys()}
            for R in m.connector:
                if not R in sub_mols:   s = sub_mols['default'].copy()
                else:                   s = sub_mols[R].copy()
                lconn = (s.atoms[s.connector[0]], s.atoms[s.connector[1]])
                m.substitute(m.connector[R], s, lconn)

            sorted_Rnames = list(sorted(m.substituents.keys()))
            sorted_R = [m.substituents[R] for R in sorted_Rnames]
            mpath = os.path.join(paths.xyz, template_name, '_'.join(sorted_R), f'{m.name}.xyz')
            m.path = mpath
            os.makedirs(os.path.dirname(mpath), exist_ok=True)
            [m.flags.append(r) for r in [R+"="+n for R, n in m.substituents.items()]]
            print(m.flags)
            comment = ", ".join(m.flags)
            utility.write_mol(m, m.path, comment=comment)

    return mols


if __name__ == '__main__':
    mols = generate_stationary_points('no_catalyst', {'R1':'F', 'R2': 'F'})
    # for mol in mols:
    #     print(type(mol))
    #     print(mol.name)
    #     print(mol)
    #     print()
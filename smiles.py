
from rdkit import Chem
from rdkit.Chem import RWMol, SanitizeMol, rdchem

def is_valence_ok(atom_element, existing_bonds, new_bond_order=1):
    max_valences = {'H': 1, 'O': 2, 'C': 4}
    return (existing_bonds + new_bond_order) <= max_valences.get(atom_element, 4)

def add_atom_and_bond(molecule, atom_index, element):
    atom = molecule.GetAtomWithIdx(atom_index)
    existing_bonds = sum([bond.GetBondTypeAsDouble() for bond in atom.GetBonds()])
    if not is_valence_ok(atom.GetSymbol(), existing_bonds):
        return None

    new_atom_index = molecule.AddAtom(Chem.Atom(element))
    molecule.AddBond(atom_index, new_atom_index, rdchem.BondType.SINGLE)
    try:
        Chem.SanitizeMol(molecule)
        return new_atom_index
    except:
        molecule.RemoveAtom(new_atom_index)  # Rollback if invalid
        return None

def generate_molecules(max_depth, current_molecule, current_depth=0, molecules=None):
    if molecules is None:
        molecules = []

    smiles = Chem.MolToSmiles(current_molecule)
    molecules.append(smiles)

    if current_depth <= max_depth:
        for atom_index in range(current_molecule.GetNumAtoms()):
            for element in ['C', 'O', 'H']:
                existing_bonds = sum([bond.GetBondTypeAsDouble() for bond in current_molecule.GetAtomWithIdx(atom_index).GetBonds()])
                if is_valence_ok(element, existing_bonds):
                    new_molecule = Chem.RWMol(current_molecule)
                    if add_atom_and_bond(new_molecule, atom_index, element) is not None:
                        generate_molecules(max_depth, new_molecule, current_depth + 1, molecules)

    return list(set(molecules))  # Élimine les doublons avant de retourner la liste

if __name__ == '__main__':
    max_depth = 7
    init_molecule = Chem.RWMol()
    init_molecule.AddAtom(Chem.Atom('C'))  # Starting with a carbon atom
    generated_molecules = generate_molecules(max_depth, init_molecule)
    for smiles in generated_molecules:
        print(smiles)
    print(f'Generated {len(generated_molecules)} unique molecules.')
    file_name = f'smiles_depth_{max_depth}.txt'

    with open(file_name, 'w') as file:
        for smiles in generated_molecules:
            file.write(smiles + '\n')

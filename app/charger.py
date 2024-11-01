from rdkit import Chem
from rdkit.Chem import rdChemReactions
import pandas as pd
import os

def handler(signum, frame):
    print("Timeout!")
    raise Exception("Timeout!")

def get_rules():
    local_csv_path = "data/definition_prot_rules.csv"
    try:
        df_rules = pd.read_csv(local_csv_path, index_col=0)
    except FileNotFoundError:
        raise FileNotFoundError(f"Could not find the file {local_csv_path}. Please ensure it exists.")
    
    rules = []
    for idx, row in df_rules.iterrows():
        rules.append([[row["smirks"]], row["type"], row["pka"], idx])
    return rules

def reformat_rules(rules):
    rule_counter = 0
    for rule in rules:
        trans_counter = 0
        for transformation in rule[0]:
            educt, product = transformation.split('>>')
            # Primary amines
            product = product.replace('[NH3+]', '[N+]([H])([H])[H]')
            # Secondary amines
            product = product.replace('[NH2+]', '[N+]([H])([H])')
            # Tertiary amines
            product = product.replace('[NH+]', '[N+]([H])')
            product = product.replace('[nH+]', '[n+]([H])')
            product = product.replace('[nH]', '[n]([H])')
            reaction = educt + '>>' + product
            rules[rule_counter][0][trans_counter] = reaction
            trans_counter += 1
        rule_counter += 1
    return rules

def charge_molecules(
    input_data,
    ph=7.0,
    ph_range=1.0,
    is_string=True,
    code=None,
    debug=False
):
    if not is_string:
        if not os.path.isfile(input_data):
            raise FileNotFoundError(f"File path not valid or file does not exist: {input_data}")
        infile_list = input_data.split('.')
        file_type = infile_list[-1].lower()
        if file_type == 'smi' or file_type == 'csv':
            try:
                mol_list = Chem.SmilesMolSupplier(input_data, delimiter="\t")
            except Exception as e:
                raise IOError(f"Unable to open file: {input_data}") from e
        else:
            raise ValueError('Unknown file type. Supported types are .smi and .csv')
    else:
        if isinstance(input_data, str):
            mol_input = Chem.MolFromSmiles(input_data)
            if mol_input is None:
                raise ValueError("Invalid SMILES string provided.")
            mol_input.SetProp("_Name", "input_molecule")
            mol_list = [mol_input]
        elif isinstance(input_data, list):
            mol_list = []
            for smi in input_data:
                mol = Chem.MolFromSmiles(smi)
                if mol is None:
                    raise ValueError(f"Invalid SMILES string in list: {smi}")
                mol.SetProp("_Name", "input_molecule")
                mol_list.append(mol)
        else:
            raise TypeError("input_data must be a string or a list of strings when is_string is True")

    rules = get_rules()

    reaction, type_idx, pka = 0, 1, 2
    counter = 0
    for rule in rules:
        if rule[type_idx] == 'acid':
            if (rule[pka] < (ph - ph_range)) or (ph_range == 0 and rule[pka] <= (ph - ph_range)):
                direction = 'forward'
            elif rule[pka] > (ph + ph_range):
                direction = 'backwards'
            else:
                direction = 'both'
        else:
            if rule[pka] > (ph + ph_range):
                direction = 'forward'
            elif rule[pka] < (ph - ph_range) or (ph_range == 0 and rule[pka] >= (ph + ph_range)):
                direction = 'backwards'
            else:
                direction = 'both'

        if direction != 'forward':
            start = rule[reaction][0].find('>')
            backward_direction = rule[reaction][0][start+2:] + '>>' + rule[reaction][0][:start]
        if direction == 'backwards':
            rules[counter][reaction][0] = backward_direction
        elif direction == 'both':
            rules[counter][reaction].append(backward_direction)
        counter += 1
    rules = reformat_rules(rules)

    charged_smiles_list = []

    for mol in mol_list:
        try:
            frag_mol = Chem.GetMolFrags(mol, True, True)
        except Exception:
            print("Could not read molecule")
            continue
        if len(frag_mol) > 1:
            length = 0
            keep_frag = frag_mol[0]
            for frag in frag_mol:
                if frag.GetNumAtoms() > length:
                    length = frag.GetNumAtoms()
                    keep_frag = frag
            name = mol.GetProp("_Name")
            mol = keep_frag
            mol.SetProp("_Name", name)
        Chem.AssignStereochemistry(mol)
        smiles_list = [Chem.MolToSmiles(mol)]
        for rule in rules:
            tmp_list = []
            for smiles in smiles_list:
                for transformation in rule[0]:
                    work_mol = Chem.MolFromSmiles(smiles)
                    try:
                        work_mol = Chem.rdmolops.AddHs(work_mol)
                    except Exception:
                        print(f"Could not add hydrogen atoms, skipping {name}")
                        break
                    to_do_list = [work_mol]
                    while len(to_do_list) > 0:
                        new_to_do = []
                        new_to_do_smiles = []
                        for work_mol in to_do_list:
                            umr = rdChemReactions.ReactionFromSmarts(transformation)
                            try:
                                product = umr.RunReactant(work_mol, 0)
                            except Exception:
                                work_mol = Chem.MolFromSmiles(smiles)
                                if smiles not in tmp_list:
                                    tmp_list.append(smiles)
                                continue
                            counter = 0
                            if len(product) == 0 and smiles not in tmp_list:
                                tmp_list.append(smiles)
                            for ps in product:
                                for prod in ps:
                                    counter += 1
                                    try:
                                        Chem.SanitizeMol(prod)
                                    except Exception:
                                        if smiles not in tmp_list:
                                            tmp_list.append(smiles)
                                        continue
                                    if prod.GetNumAtoms() + 3 < work_mol.GetNumAtoms():
                                        if smiles not in tmp_list:
                                            tmp_list.append(smiles)
                                        continue
                                    if umr.IsMoleculeReactant(prod):
                                        prod = Chem.rdmolops.RemoveHs(prod)
                                        tmp_smi = Chem.MolToSmiles(prod)
                                        if tmp_smi not in new_to_do_smiles:
                                            new_to_do_smiles.append(tmp_smi)
                                            new_to_do.append(prod)
                                        if debug:
                                            print(f"Reaction can still proceed on molecule: {tmp_smi}")
                                    else:
                                        prod = Chem.rdmolops.RemoveHs(prod)
                                        tmp_smi = Chem.MolToSmiles(prod)
                                        if tmp_smi not in tmp_list:
                                            tmp_list.append(tmp_smi)
                                            if debug:
                                                print(f"Final product of reaction: {tmp_smi}")
                        to_do_list = new_to_do
            smiles_list = tmp_list.copy()

        if debug:
            print('Final charged forms:')
            print(smiles_list)
        charged_smiles_list.extend(smiles_list)
    return charged_smiles_list
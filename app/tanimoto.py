from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity
import csv
import os

def calculate_tanimoto_similarity(template_smiles, comparison_smiles):
    # Ensure directory exists
    if not os.path.exists("tmp"):
        os.makedirs("tmp")

    template_molecules = []
    for smile in template_smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            print(f"Invalid SMILES string in template: {smile}")
        else:
            template_molecules.append(mol)

    comparison_molecules = []
    for smile in comparison_smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            print(f"Invalid SMILES string in comparison: {smile}")
        else:
            comparison_molecules.append(mol)

    if not template_molecules or not comparison_molecules:
        print("No valid molecules found.")
        return []

    template_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2) for mol in template_molecules]
    comparison_fps = [AllChem.GetMorganFingerprintAsBitVect(mol, 2) for mol in comparison_molecules]

    csv_files = []
    for i, template_fp in enumerate(template_fps):
        file_name = f'tmp/tanimoto_scores_template_{i+1}.csv'  # Save in tmp/ directory
        with open(file_name, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Pair', 'Tanimoto Score'])
            for j, comparison_fp in enumerate(comparison_fps):
                score = TanimotoSimilarity(template_fp, comparison_fp)
                writer.writerow([f'{template_smiles[i]}-{comparison_smiles[j]}', score])
        csv_files.append(file_name)

    return csv_files

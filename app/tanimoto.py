from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.DataStructs.cDataStructs import TanimotoSimilarity

def calculate_tanimoto_similarity(template_smiles, comparison_smiles):
    template_molecules = []
    for smile in template_smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            print(f"Invalid SMILES string in template: {smile}")
        else:
            template_molecules.append((smile, mol))

    comparison_molecules = []
    for smile in comparison_smiles:
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            print(f"Invalid SMILES string in comparison: {smile}")
        else:
            comparison_molecules.append((smile, mol))

    if not template_molecules or not comparison_molecules:
        print("No valid molecules found.")
        return []

    template_fps = [(smile, AllChem.GetMorganFingerprintAsBitVect(mol, 2)) for (smile, mol) in template_molecules]
    comparison_fps = [(smile, AllChem.GetMorganFingerprintAsBitVect(mol, 2)) for (smile, mol) in comparison_molecules]

    scores = []
    for t_smile, t_fp in template_fps:
        for c_smile, c_fp in comparison_fps:
            score = TanimotoSimilarity(t_fp, c_fp)
            scores.append({
                'template_smile': t_smile, 
                'comparison_smile': c_smile, 
                'tanimoto_score': round(score, 4)
            })

    return scores

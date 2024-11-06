from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, QED, rdFreeSASA, AllChem, Lipinski, Descriptors3D, rdMolDescriptors
import inspect
from io import BytesIO
import base64

def get_all_descriptors():
    all_descriptors = {
        'chem': {name: func for name, func in inspect.getmembers(Descriptors, inspect.isfunction) if filter_method(name)},
        'lipinski': {name: func for name, func in inspect.getmembers(Lipinski, inspect.isfunction) if filter_method(name)},
        'crippen': {name: func for name, func in inspect.getmembers(Crippen, inspect.isfunction) if filter_method(name)},
        'qed': {name: func for name, func in inspect.getmembers(QED, inspect.isfunction) if filter_method(name)},
        'rdfreesasa': {name: func for name, func in inspect.getmembers(rdFreeSASA, inspect.isfunction) if filter_method(name)},
        'descriptor3d': {name: func for name, func in inspect.getmembers(Descriptors3D, inspect.isfunction) if filter_method(name)},
        'rdmoldescriptors': {name: func for name, func in inspect.getmembers(rdMolDescriptors, callable) if filter_method(name)}
    }
    return all_descriptors

def filter_method(name):
    """Filter out non-working methods."""
    not_working_methods = [
        "rundoctest", "auto", 'namedtuple', 'setdescriptordersion', '_init', '_readpatts', 
        'ads', 'AtomPairsParameters', 'CalcChiNn', 'CalcChiNv', 'CustomProp_VSA_', 
        'GetAtomFeatures', 'GetAtomPairAtomCode', 'GetAtomPairCode', 
        'GetHashedMorganFingerprint', 'GetMorganFingerprint', 'GetUSR', 
        'MakePropertyRangeQuery', 'NumRotatableBondsOptions', 'Properties', 
        'PropertyFunctor', 'PropertyRangeQuery'
    ]
    if name[0] == "_":
        return False
    filtered_descriptors = not any(nm in name.lower() for nm in not_working_methods)
    return filtered_descriptors

def compute_descriptors(smiles, selected_options):
    molecule = Chem.MolFromSmiles(smiles)
    descriptors = {}
    img_str = None
    all_descriptors = get_all_descriptors()

    if molecule is not None:
        for option in selected_options:
            descriptors['SMILES'] = smiles
            method_name = option
            if method_name == 'CalcMolDescriptors':
                full_descriptors = Descriptors.CalcMolDescriptors(molecule)
                for key, value in full_descriptors.items():
                    descriptors[key] = value
                continue
            elif method_name == 'CalcMolDescriptors3D':
                if molecule.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(molecule, AllChem.ETKDG())
                full_descriptors = Descriptors3D.CalcMolDescriptors3D(molecule)
                for key, value in full_descriptors.items():
                    descriptors[key] = value
                continue
            elif method_name in all_descriptors['chem']:
                descriptors[method_name] = all_descriptors['chem'][method_name](molecule)
            elif method_name in all_descriptors['lipinski']:
                descriptors[method_name] = all_descriptors['lipinski'][method_name](molecule)
            elif method_name in all_descriptors['crippen']:
                descriptors[method_name] = all_descriptors['crippen'][method_name](molecule)
            elif method_name in all_descriptors['qed']:
                descriptors[method_name] = all_descriptors['qed'][method_name](molecule)
            elif method_name in all_descriptors['rdfreesasa']:
                descriptors[method_name] = all_descriptors['rdfreesasa'][method_name](molecule)
            elif method_name in all_descriptors['descriptor3d']:
                AllChem.EmbedMolecule(molecule, AllChem.ETKDG())
                descriptors[method_name] = all_descriptors['descriptor3d'][method_name](molecule)
            elif method_name in all_descriptors['rdmoldescriptors']:
                AllChem.EmbedMolecule(molecule, AllChem.ETKDG())
                descriptors[method_name] = all_descriptors['rdmoldescriptors'][method_name](molecule)

        if 'Image' in selected_options:
            img = Draw.MolToImage(molecule)
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
            descriptors['Image'] = img_str
        if 'FreeSASA' in selected_options:
            AllChem.EmbedMolecule(molecule, AllChem.ETKDG())
            radii = rdFreeSASA.classifyAtoms(molecule)
            sasa_opts = rdFreeSASA.SASAOpts(rdFreeSASA.SASAAlgorithm.ShrakeRupley, rdFreeSASA.SASAClassifier.OONS)
            descriptors['FreeSASA'] = rdFreeSASA.CalcSASA(molecule, radii, confIdx=-1, opts=sasa_opts)
    else:
        descriptors['Error'] = f"Invalid SMILES: {smiles}"

    return descriptors, img_str

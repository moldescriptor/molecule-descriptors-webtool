from flask import Flask, render_template, request, send_from_directory, render_template_string, Response
from rdkit import Chem, rdBase
from rdkit.Chem import Draw, Descriptors, Crippen, QED, rdFreeSASA, AllChem, Lipinski, Descriptors3D, rdMolDescriptors
from collections import OrderedDict
from io import StringIO
import csv
import re
import inspect
from io import BytesIO
import base64
import json
import os

app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 MB max-limit

def load_descriptor_synonyms(file_name="descriptor_synonyms.csv"):
    synonym_map = {}
    descriptor_synonyms = {}
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    csv_path = os.path.join(base_dir, 'data', file_name)
    
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Synonym CSV file not found at {csv_path}")
    
    with open(csv_path, 'r', newline='', encoding='utf-8') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            descriptor = row['Descriptor'].strip()
            synonyms = [row[f'Synonym {i}'].strip() for i in range(1, 6) if row.get(f'Synonym {i}')]
            synonyms = [syn for syn in synonyms if syn]
            
            # **Always add the descriptor, even if synonyms are empty**
            descriptor_synonyms[descriptor] = synonyms
            
            for synonym in synonyms:
                synonym_lower = synonym.lower()
                synonym_map[synonym_lower] = descriptor
            descriptor_lower = descriptor.lower()
            synonym_map[descriptor_lower] = descriptor
    return synonym_map, descriptor_synonyms

# **Initialize synonym_map and descriptor_synonyms globally**
synonym_map, descriptor_synonyms = load_descriptor_synonyms()

@app.route('/')
def index():
    all_descriptors = get_all_descriptors()
    descriptorsAndSynonyms = build_descriptors_and_synonyms(descriptor_synonyms)
    
    try:
        json_descriptorsAndSynonyms = json.dumps(descriptorsAndSynonyms)
    except TypeError as e:
        print(f"Serialization error: {e}")
        json_descriptorsAndSynonyms = '{}'
    
    return render_template('index.html', 
                           all_descriptors=all_descriptors, 
                           descriptorsAndSynonyms=json_descriptorsAndSynonyms)

@app.context_processor
def inject_defaults():
    return {'result_available': False}

@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    global_descriptors = get_all_descriptors()
    selected_options = request.form.getlist('displayOptions')

    final_selected_options = []
    for option in selected_options:
        option_lower = option.lower()
        if option_lower in synonym_map:
            method_name = synonym_map[option_lower]
            final_selected_options.append(method_name)
        else:
            final_selected_options.append(option)

    exclude_invalid = request.form.get('excludeInvalid') == 'true'

    file = request.files.get('csvFile')
    if file and file.filename != '':
        if not file.filename.endswith('.csv'):
            all_descriptors = get_all_descriptors()
            descriptorsAndSynonyms = build_descriptors_and_synonyms(descriptor_synonyms)
            return render_template('index.html', 
                                   error="Invalid file type! Please upload a .csv file.",
                                   all_descriptors=global_descriptors,
                                   descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms))
        smiles_list = []
        csv_file = csv.reader(file.stream.read().decode("utf-8").splitlines())
        for row in csv_file:
            if row:
                smiles_list.append(row[0])
    else:
        smiles_input = request.form.get('inputField', '')
        smiles_list = [s.strip() for s in re.split(r'[,.]', smiles_input) if s.strip() != '']

    if exclude_invalid:
        smiles_list = [smiles for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    computed_descriptors = []
    for smiles in smiles_list:
        descriptors, img_str = compute_descriptors(smiles, final_selected_options)
        computed_descriptors.append(descriptors)

    descriptorsAndSynonyms = build_descriptors_and_synonyms(descriptor_synonyms)

    return render_template('index.html', 
                           descriptors_list=computed_descriptors, 
                           all_descriptors=global_descriptors, 
                           descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms),
                           result_available=True)

def build_descriptors_and_synonyms(descriptor_synonyms):
    return descriptor_synonyms

def generate_csv(data):
    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(data.keys())
    writer.writerow(data.values())
    return output.getvalue()

def get_all_descriptors():
    all_descriptors = {
        'chem': {name: func for name, func in inspect.getmembers(Descriptors, inspect.isfunction) if filter_method(name)},
        'lipinski': {name: func for name, func in inspect.getmembers(Lipinski, inspect.isfunction) if filter_method(name)},
        'crippen': {name: func for name, func in inspect.getmembers(Crippen, inspect.isfunction) if filter_method(name)},
        'qed': {name: func for name, func in inspect.getmembers(QED, inspect.isfunction) if filter_method(name)},
        'rdfreesasa': {name: func for name, func in inspect.getmembers(rdFreeSASA, inspect.isfunction) if filter_method(name)},
        'descriptor3d': {name: func for name, func in inspect.getmembers(Descriptors3D, inspect.isfunction) if filter_method(name)},
        'rdmoldescriptors': {name: func for name, func in inspect.getmembers(rdMolDescriptors, callable) if filter_method(name) and name != 'CalcMolDescriptors'}
    }

    # Manually added methods
    all_descriptors['rdfreesasa']["FreeSASA"] = "FreeSASA"

    return all_descriptors

def filter_method(name):
    not_working_methods = [
        "rundoctest", "auto", 'namedtuple', 'setdescriptordersion', '_init', '_readpatts', 
        'ads', 'AtomPairsParameters', 'CalcChiNn', 'CalcChiNv', 'CustomProp_VSA_', 
        'GetAtomFeatures', 'GetAtomPairAtomCode', 'GetAtomPairCode', 
        'GetHashedMorganFingerprint', 'GetMorganFingerprint', 'GetUSR', 
        'MakePropertyRangeQuery', 'NumRotatableBondsOptions', 'Properties', 
        'PropertyFunctor', 'PropertyRangeQuery'
    ]
    for not_working in not_working_methods:
        if not_working.lower() in name.lower():
            return False
    return True

all_descriptors = get_all_descriptors()

def compute_descriptors(smiles, selected_options):
    molecule = Chem.MolFromSmiles(smiles)
    descriptors = {}
    img_str = None

    if molecule is not None:
        descriptors['SMILES'] = smiles

        if 'CalcMolDescriptors' in selected_options:
            full_descriptors = Descriptors.CalcMolDescriptors(molecule)
            for key, value in full_descriptors.items():
                descriptors[key] = value
            selected_options = [opt for opt in selected_options if opt != 'CalcMolDescriptors']
        
        for option in selected_options:
            method_name = option
            if method_name in all_descriptors['chem']:
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

def get_atom_counts(molecule):
    atom_counts = OrderedDict()
    mol_with_hydrogens = Chem.AddHs(molecule)
    for atom in mol_with_hydrogens.GetAtoms():
        symbol = atom.GetSymbol()
        key = "Number of " + symbol + " atoms"
        atom_counts[key] = atom_counts.get(key, 0) + 1
    
    atom_counts["Number of atoms total:"] = sum(atom_counts.values())
    
    return atom_counts

@app.route('/feedback')
def feedback():
    return render_template('feedback.html')

@app.route('/draw')
def draw():
    return render_template('draw.html')

#@app.route('/docs')
#def docs():
#    return render_template('docs.html')

@app.route('/about')
def about():
    return render_template('about.html')

@app.route('/editor')
def editor():
    dir_path = os.path.dirname(os.path.realpath(__file__))
    file_path = os.path.join(dir_path, 'frontend', 'editor.html')

    with open(file_path, "r") as f:
        content = f.read()
    return render_template_string(content)

@app.route('/editor/gui/<path:filename>')
def editor_resources(filename):
    return send_from_directory('frontend/gui', filename)

@app.route('/download_csv', methods=['POST'])
def download_csv():
    selected_options = request.form.getlist('displayOptions')
    exclude_invalid = request.form.get('excludeInvalid') == 'true'
    
    smiles_input = request.form.get('inputField', '')
    smiles_list = [s.strip() for s in smiles_input.split(',')]

    if exclude_invalid:
        smiles_list = [smiles for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]
    
    output = StringIO()
    writer = csv.writer(output)

    rdkit_version = rdBase.rdkitVersion

    all_descriptors = [compute_descriptors(smiles, selected_options)[0] for smiles in smiles_list]

    all_keys = set()
    for desc in all_descriptors:
        all_keys.update(desc.keys())

    all_keys.discard('Image')
    all_keys.discard('SMILES')

    all_keys = sorted(list(all_keys))

    atom_columns = sorted([key for key in all_keys if key.startswith("Number of")])
    
    other_keys = [key for key in all_keys if key not in atom_columns and key != 'Error']
    ordered_keys = ['SMILES'] + atom_columns + other_keys + (['Error'] if 'Error' in all_keys else [])

    writer.writerow(ordered_keys)

    for desc in all_descriptors:
        writer.writerow([desc.get(key, '') for key in ordered_keys])

    csv_data = output.getvalue()

    return Response(
        csv_data,
        mimetype="text/csv",
        headers={"Content-disposition": f"attachment; filename=data_RDKit_{rdkit_version}.csv"}
    )

if __name__ == '__main__':
    app.run(debug=True)

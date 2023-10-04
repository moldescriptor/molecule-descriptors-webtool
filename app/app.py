from flask import Flask, render_template, request, send_from_directory, render_template_string, Response
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, QED, rdFreeSASA, AllChem
from collections import OrderedDict
from io import StringIO
import csv


from io import BytesIO
import base64

#from app.routes import main  # Replace 'your_folder_name' with the actual folder name where routes.py is located

#app = Flask(__name__)
app = Flask(__name__, template_folder='templates', static_folder='static')
#app = Flask(__name__, static_folder='app/frontend', static_url_path='/app/frontend')


#app.register_blueprint(main)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    smiles_input = request.form.get('inputField', '')
    selected_options = request.form.getlist('displayOptions')  # This line was missing
    smiles_list = [s.strip() for s in smiles_input.split(',')]
    
    all_descriptors = []
    for smiles in smiles_list:
        descriptors, img_str = compute_descriptors(smiles, selected_options)
        all_descriptors.append(descriptors)

    return render_template('index.html', descriptors_list=all_descriptors)



def generate_csv(data):
    output = StringIO()
    writer = csv.writer(output)

    # Write headers
    writer.writerow(data.keys())

    # Write data
    writer.writerow(data.values())

    return output.getvalue()


def compute_descriptors(smiles, selected_options):
    molecule = Chem.MolFromSmiles(smiles)
    descriptors = {}
    img_str = None  # Initialize img_str here

    if molecule is not None:
        descriptors['SMILES'] = smiles
        if 'Image' in selected_options:
            img = Draw.MolToImage(molecule)
            buffered = BytesIO()
            img.save(buffered, format="PNG")
            img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
            descriptors['Image'] = img_str
        if 'MolecularWeight' in selected_options:
            descriptors['MolecularWeight'] = Descriptors.ExactMolWt(molecule)
        if 'PSA' in selected_options:
            descriptors['PSA'] = Descriptors.TPSA(molecule)
        if 'clogP' in selected_options:
            descriptors['clogP'] = Crippen.MolLogP(molecule)
        if 'QED' in selected_options:
            descriptors['QED'] = QED.qed(molecule)
        if 'numberOfAtoms' in selected_options:
            atom_counts = get_atom_counts(molecule)
            descriptors.update(atom_counts)
                #descriptors['numberOfAtoms'] = molecule.GetNumAtoms()
        if 'FreeSASA' in selected_options:
            # 1. Generate 3D coordinates for the molecule
            AllChem.EmbedMolecule(molecule, AllChem.ETKDG())

            # 2. Classify atoms and get radii
            radii = rdFreeSASA.classifyAtoms(molecule)

            # 3. Define SASA options: Using ShrakeRupley algorithm and OONS classifier as an example
            sasa_opts = rdFreeSASA.SASAOpts(rdFreeSASA.SASAAlgorithm.ShrakeRupley, rdFreeSASA.SASAClassifier.OONS)

            # 4. Compute the SASA and store in descriptors dictionary
            descriptors['FreeSASA'] = rdFreeSASA.CalcSASA(molecule, radii, confIdx=-1, opts=sasa_opts)
    else:
        descriptors['Error'] = "Invalid SMILES"

    return descriptors, img_str


def get_atom_counts(molecule):
    atom_counts = OrderedDict()
    # Convert implicit hydrogens to explicit ones
    mol_with_hydrogens = Chem.AddHs(molecule)
    for atom in mol_with_hydrogens.GetAtoms():
        symbol = atom.GetSymbol()
        key = "Number of "+symbol+" atoms"
        atom_counts[key] = atom_counts.get(key, 0) + 1
    
    atom_counts["Number of atoms total:"] = sum(atom_counts.values())
    
    return atom_counts



@app.route('/draw')
def draw():
    return render_template('draw.html')

@app.route('/docs')
def docs():
    return render_template('docs.html')

@app.route('/about')
def about():
    return render_template('about.html')


import os

@app.route('/editor')
def editor():
    # Get the directory of the current script
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
    smiles_input = request.form.get('inputField', '')
    smiles_list = [s.strip() for s in smiles_input.split(',')]
    
    output = StringIO()
    writer = csv.writer(output)
    first = True

    for smiles in smiles_list:
        descriptors, _ = compute_descriptors(smiles, selected_options)
        descriptors.pop('Image', None)

        # Write headers only for the first SMILES
        if first:
            writer.writerow(descriptors.keys())
            first = False

        # Write data for each SMILES
        writer.writerow(descriptors.values())

    csv_data = output.getvalue()

    return Response(
        csv_data,
        mimetype="text/csv",
        headers={"Content-disposition": "attachment; filename=data.csv"}
    )



if __name__ == '__main__':
    app.run(debug=True)

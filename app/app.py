from flask import Flask, render_template, request
from rdkit import Chem
from rdkit.Chem import Descriptors

app = Flask(__name__)

@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')

@app.route('/identify_molecule', methods=['POST'])
def identify_molecule():
    smiles = request.form.get('inputField', '')
    selected_options = request.form.getlist('displayOptions')
    
    molecule = Chem.MolFromSmiles(smiles)
    
    descriptors = {}
    if molecule is not None:
        if 'MolecularWeight' in selected_options:
            descriptors['MolecularWeight'] = Descriptors.ExactMolWt(molecule)
        if 'PSA' in selected_options:
            descriptors['PSA'] = Descriptors.TPSA(molecule)
        # Add other descriptors here based on selected_options
    else:
        descriptors['Error'] = "Invalid SMILES"
    
    return render_template('index.html', descriptors=descriptors)

if __name__ == '__main__':
    app.run(debug=True)
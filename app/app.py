from flask import Flask, render_template, request, send_from_directory, render_template_string, Response, send_file
from rdkit import Chem, rdBase
import json
import os

from app.descriptors import get_all_descriptors, compute_descriptors
from app.tanimoto import calculate_tanimoto_similarity
from app.utils import load_descriptor_synonyms, build_descriptors_and_synonyms, process_selected_options, get_smiles_list, generate_csv_data

app = Flask(__name__, template_folder='templates', static_folder='static')
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024  # 16 MB max-limit

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

    final_selected_options = process_selected_options(selected_options, synonym_map)
    exclude_invalid = request.form.get('excludeInvalid') == 'true'

    smiles_list, error = get_smiles_list(request, exclude_invalid)
    if error:
        descriptorsAndSynonyms = build_descriptors_and_synonyms(descriptor_synonyms)
        return render_template('index.html', 
                               error=error,
                               all_descriptors=global_descriptors,
                               descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms))

    computed_descriptors = [compute_descriptors(smiles, final_selected_options)[0] for smiles in smiles_list]
    descriptorsAndSynonyms = build_descriptors_and_synonyms(descriptor_synonyms)

    return render_template('index.html', 
                           descriptors_list=computed_descriptors, 
                           all_descriptors=global_descriptors, 
                           descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms),
                           result_available=True)

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
    
    all_descriptors = [compute_descriptors(smiles, selected_options)[0] for smiles in smiles_list]
    csv_data = generate_csv_data(all_descriptors, exclude_keys=['Image', 'SMILES'])

    return Response(
        csv_data,
        mimetype="text/csv",
        headers={"Content-disposition": f"attachment; filename=data_RDKit_{rdBase.rdkitVersion}.csv"}
    )

@app.route('/tanimoto', methods=['POST'])
def tanimoto_similarity():
    template_smiles = request.form['template_smiles'].split()
    comparison_smiles = request.form['comparison_smiles'].split()

    # Calculate Tanimoto similarity and generate CSV files
    csv_files = calculate_tanimoto_similarity(template_smiles, comparison_smiles)

    if not csv_files:
        return "No valid SMILES strings or no Tanimoto scores calculated.", 400

    # Use absolute path to the CSV file
    return send_file(os.path.abspath(csv_files[0]), as_attachment=True, download_name="tanimoto_scores.csv")


if __name__ == '__main__':
    app.run(debug=True)

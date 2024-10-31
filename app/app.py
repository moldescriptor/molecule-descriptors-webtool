from flask import Flask, render_template, request, send_from_directory, render_template_string, Response
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
    calculate_tanimoto = 'calculateTanimoto' in request.form

    smiles_input = request.form.get('inputField', '')

    descriptorsAndSynonyms = build_descriptors_and_synonyms(descriptor_synonyms)

    if calculate_tanimoto:
        if '||' not in smiles_input:
            error = "Please separate template and comparison SMILES strings using '||'."
            return render_template('index.html', 
                                   error=error,
                                   all_descriptors=global_descriptors,
                                   descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms))
        template_str, comparison_str = smiles_input.split('||', 1)
        template_smiles = [s.strip() for s in template_str.strip().split(',')]
        comparison_smiles = [s.strip() for s in comparison_str.strip().split(',')]

        tanimoto_scores = calculate_tanimoto_similarity(template_smiles, comparison_smiles)
        
        print("Form Data:", request.form)
        print("calculate_tanimoto:", calculate_tanimoto)
        print("exclude_invalid:", exclude_invalid)

        if not tanimoto_scores:
            error = "No valid SMILES strings or no Tanimoto scores calculated."
            return render_template('index.html', 
                                   error=error,
                                   all_descriptors=global_descriptors,
                                   descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms))

        return render_template('index.html', 
                               tanimoto_scores=tanimoto_scores, 
                               all_descriptors=global_descriptors, 
                               descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms),
                               result_available=True)
    else:
        smiles_list, error = get_smiles_list(request, exclude_invalid)
        if error:
            return render_template('index.html', 
                                   error=error,
                                   all_descriptors=global_descriptors,
                                   descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms))

        computed_descriptors = [compute_descriptors(smiles, final_selected_options)[0] for smiles in smiles_list]

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

@app.route('/download_tanimoto_csv', methods=['POST'])
def download_tanimoto_csv():
    tanimoto_scores_json = request.form.get('tanimoto_scores')
    tanimoto_scores = json.loads(tanimoto_scores_json)

    csv_data = "Template SMILES,Comparison SMILES,Tanimoto Score\n"
    for item in tanimoto_scores:
        csv_data += f"{item['template_smile']},{item['comparison_smile']},{item['tanimoto_score']}\n"

    return Response(
        csv_data,
        mimetype="text/csv",
        headers={"Content-disposition": "attachment; filename=tanimoto_scores.csv"}
    )

if __name__ == '__main__':
    app.run(debug=True)

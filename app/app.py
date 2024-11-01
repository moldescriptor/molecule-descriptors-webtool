from flask import Flask, render_template, request, send_from_directory, render_template_string, Response
from rdkit import Chem, rdBase
from collections import defaultdict
import json
import os

from app.descriptors import get_all_descriptors, compute_descriptors
from app.tanimoto import calculate_tanimoto_similarity
from app.utils import load_descriptor_synonyms, build_descriptors_and_synonyms, process_selected_options, get_smiles_list, generate_csv_data
from app.charger import charge_molecules

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
                           descriptorsAndSynonyms=json_descriptorsAndSynonyms,
                           result_available=False)

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

    charge_smiles = 'chargeSmiles' in request.form
    charge_smiles_selected = charge_smiles

    ph_value = request.form.get('phValue', type=float)
    ph_range = request.form.get('phRange', type=float)

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

        smiles_info = {}
        for smiles in smiles_list:
            smiles_info[smiles] = {
                'Original_SMILES': smiles,
                'Charged': False
            }
        
        if charge_smiles:
            results = []
            for original_smiles in smiles_list:
                try:
                    charged_smiles_list = charge_molecules(
                        input_data=original_smiles,
                        ph=ph_value if ph_value is not None else 7.0,
                        ph_range=ph_range if ph_range is not None else 1.0,
                        is_string=True,
                        debug=False
                    )
                    if not charged_smiles_list:
                        charged_smiles_list = [original_smiles]
                except Exception as e:
                    error = f"Error charging molecules: {e}"
                    return render_template('index.html', 
                                           error=error,
                                           all_descriptors=global_descriptors,
                                           descriptorsAndSynonyms=json.dumps(descriptorsAndSynonyms))
                charged_smiles_list = list(set(charged_smiles_list))

                for charged_smiles in charged_smiles_list:
                    if charged_smiles not in smiles_info:
                        smiles_info[charged_smiles] = {
                            'Original_SMILES': original_smiles,
                            'Charged': True
                        }

                results.append({
                    'Original_SMILES': original_smiles,
                    'Charged_SMILES_List': charged_smiles_list
                })
        else:
            results = []

        descriptors_list = []
        for smiles, info in smiles_info.items():
            descriptor = compute_descriptors(smiles, final_selected_options)[0]
            if descriptor:
                descriptor['Original_SMILES'] = info['Original_SMILES']
                descriptor['Charged'] = info['Charged']
                descriptors_list.append(descriptor)
        grouped_descriptors = defaultdict(list)
        for descriptor in descriptors_list:
            original_smiles = descriptor['Original_SMILES']
            grouped_descriptors[original_smiles].append(descriptor)

        return render_template('index.html', 
                            grouped_descriptors=grouped_descriptors,
                            descriptors_list=descriptors_list,
                            charge_smiles_selected=charge_smiles_selected,
                            results=results,
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
    descriptors_list_json = request.form.get('descriptors_list')
    if not descriptors_list_json:
        return "No descriptors data provided.", 400
    try:
        descriptors_list = json.loads(descriptors_list_json)
    except json.JSONDecodeError:
        return "Invalid descriptors data.", 400

    exclude_keys = ['Image']
    csv_data = generate_csv_data(descriptors_list, exclude_keys=exclude_keys)
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

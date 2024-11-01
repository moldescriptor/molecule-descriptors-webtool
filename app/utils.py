import os
import csv
from rdkit import Chem
import re
from io import StringIO

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

def build_descriptors_and_synonyms(descriptor_synonyms):
    return descriptor_synonyms

def process_selected_options(selected_options, synonym_map):
    """
    Process selected options to account for synonyms and return final selected options.
    """
    final_selected_options = []
    for option in selected_options:
        option_lower = option.lower()
        method_name = synonym_map.get(option_lower, option)
        final_selected_options.append(method_name)
    return final_selected_options


def get_smiles_list(request, exclude_invalid):
    """
    Get SMILES list from either CSV file upload or input field.
    """
    file = request.files.get('csvFile')
    if file and file.filename != '':
        if not file.filename.endswith('.csv'):
            return None, "Invalid file type! Please upload a .csv file."
        
        csv_file = csv.reader(file.stream.read().decode("utf-8").splitlines())
        smiles_list = [row[0] for row in csv_file if row]
    else:
        smiles_input = request.form.get('inputField', '')
        smiles_list = [s.strip() for s in re.split(r'[,.]', smiles_input) if s.strip() != '']

    if exclude_invalid:
        smiles_list = [smiles for smiles in smiles_list if Chem.MolFromSmiles(smiles) is not None]

    return smiles_list, None

def generate_csv_data(descriptors_list, exclude_keys):
    """
    Generate CSV data from a list of descriptors.
    """
    all_keys = set()
    for desc in descriptors_list:
        all_keys.update(desc.keys())

    for key in exclude_keys:
        all_keys.discard(key)

    ordered_keys = ['Original_SMILES', 'SMILES', 'Charged'] + sorted(all_keys - {'Original_SMILES', 'SMILES', 'Charged'})

    output = StringIO()
    writer = csv.writer(output)
    writer.writerow(ordered_keys)

    for desc in descriptors_list:
        writer.writerow([desc.get(key, '') for key in ordered_keys])

    return output.getvalue()
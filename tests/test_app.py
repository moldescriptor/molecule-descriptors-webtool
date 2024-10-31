import os
from io import BytesIO
from app.descriptors import get_all_descriptors, compute_descriptors

def test_homepage(client):
    response = client.get('/')
    
    # Check if the status code is 200 (OK)
    assert response.status_code == 200
    
    # Ensure the template is rendered correctly with necessary content
    assert b"MolDescriptor" in response.data
    assert b"RDKit packages" in response.data


def test_download_csv(client):
    # Prepare data to send to the endpoint
    data = {
        'inputField': "CCO,CCN",  # SMILES strings
        'displayOptions': ['MolWt', 'MolLogP'],  # Descriptors to include
        'excludeInvalid': 'false'  # Don't exclude invalid SMILES
    }
    
    # Send POST request
    response = client.post('/download_csv', data=data)

    # Check if the status code is 200 (OK)
    assert response.status_code == 200

    # Ensure the response contains a CSV file as an attachment
    assert 'attachment; filename=data_RDKit' in response.headers['Content-Disposition']

    # Optional: Validate CSV content
    csv_content = response.data.decode('utf-8')
    assert "MolWt" in csv_content
    assert "MolLogP" in csv_content

def test_invalid_file_upload(client):
    data = {
        'displayOptions': ['MolWt'],
        'excludeInvalid': 'false',
        'csvFile': (BytesIO(b"not,a,valid,CSV"), 'invalid.txt')  # Place the file in the data dict directly
    }
    
    response = client.post('/identify_molecule', data=data, content_type='multipart/form-data')

    # Check if the status code is 200 (the form reloads with an error message)
    assert response.status_code == 200

    # Ensure the response contains an error message about invalid file type
    assert b"Invalid file type! Please upload a .csv file." in response.data

def test_get_all_descriptors():
    descriptors = get_all_descriptors()
    
    assert 'chem' in descriptors
    assert 'lipinski' in descriptors
    assert 'crippen' in descriptors
    assert 'qed' in descriptors
    assert 'rdfreesasa' in descriptors
    assert 'descriptor3d' in descriptors
    assert 'rdmoldescriptors' in descriptors

    # Check if some functions are present in the descriptors list
    assert 'MolWt' in descriptors['chem']
    assert 'TPSA' in descriptors['chem']
    assert 'NumHDonors' in descriptors['lipinski']

def test_compute_descriptors_valid_smiles():
    smiles = "CCO"  # Ethanol
    selected_options = ['MolWt', 'TPSA']
    
    descriptors, img_str = compute_descriptors(smiles, selected_options)
    
    assert descriptors['SMILES'] == smiles
    assert 'MolWt' in descriptors
    assert 'TPSA' in descriptors
    assert img_str is None  # No image selected

def test_compute_descriptors_invalid_smiles():
    smiles = "invalid_smiles"
    selected_options = ['MolWt', 'TPSA']
    
    descriptors, img_str = compute_descriptors(smiles, selected_options)
    
    assert 'Error' in descriptors
    assert descriptors['Error'] == f"Invalid SMILES: {smiles}"
    assert img_str is None
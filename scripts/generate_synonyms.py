import os
import csv
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from app.app import get_all_descriptors

def save_descriptor_names_to_csv(file_name="descriptor_synonyms_empty.csv"):
    all_descriptors = get_all_descriptors()
    descriptors_list = []

    for category, descriptors in all_descriptors.items():
        for descriptor_name in descriptors.keys():
            descriptors_list.append(descriptor_name)

    base_dir = os.path.abspath(os.path.dirname(__file__))
    csv_path = os.path.join(base_dir, '..', 'data', file_name)
    
    with open(csv_path, 'w', newline='', encoding='utf-8') as csvfile:
        fieldnames = ['Descriptor', 'Synonym 1', 'Synonym 2', 'Synonym 3', 'Synonym 4', 'Synonym 5']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for descriptor in descriptors_list:
            writer.writerow({
                'Descriptor': descriptor,
                'Synonym 1': '',
                'Synonym 2': '',
                'Synonym 3': '',
                'Synonym 4': '',
                'Synonym 5': ''
            })

if __name__ == "__main__":
    save_descriptor_names_to_csv()
    print("Empty synonym CSV file created successfully.")

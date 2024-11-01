# MolDescriptor

[MolDescriptor](https://moldescriptor.com/) is a simple, user-friendly web application where users can input molecules (in a SMILES string) and use RDkit, an open source Python package, to calculate various descriptors. The descriptors are then displayed in a table on the website. The result can also be downloaded as a CSV file.

The goal of this project is to utilize RDKit to calculate descriptors and make them easily accessible to users. The descriptors are organized into tabs based on the RDKit modules they belong to. Users can also search for specific descriptors using the search tab. The descriptors are calculated using RDKit and displayed in a table on the website. The end goal is making life easier for people with a chemical background who want to calculate descriptors for their molecules, without having to worry about learning Python. 


### Dependencies
We are using Poetry for dependency management. You can find the documentation for poetry [here](https://python-poetry.org/docs/).

Install the dependencies using Poetry:

```bash
poetry install
```

#### Linting with Ruff
This project uses Ruff for linting Python code. To ensure good coding practices, run the following command:

```bash
poetry run ruff check .
```

To auto-fix and format code:
```bash
poetry run ruff check --fix .
```

Update lock file:
```bash
poetry lock
```

## Features

- SMILES Input: Directly type in or paste the SMILES notation for your molecule.
- Descriptor Tabs: Organized by RDKit modules such as Chem, Lipinski, Crippen, QED, rdFreeSASA, and AllChem.
- Search Descriptors: If you know the name, just search for the descriptor you need.
- Result Visualization: View the computed descriptor values and the visual representation of the molecule.
- CSV Support: Upload a CSV containing multiple molecules to get descriptors for all at once or download results in CSV format.

# How to use MolDescriptor
You can find an in depth tutorial on how to use MolDescriptor [here.](https://moldescriptor.github.io/moldescriptor-docs/)

## Getting Started
1. Visit the MolDescriptor [website](https://moldescriptor.com/).
2. In the input field, enter the SMILES notation of the molecule you wish to analyze.
3. Choose the descriptors you are interested in by navigating the tabs.
4. (Optional) Use the search tab to find and select specific descriptors.
5. Click on "Get Result" to compute and display the descriptors.
6. (Optional) Upload a CSV file containing SMILES notations for batch processing or download your results as a CSV.

## Running MolDescriptor locally

To run MolDescriptor on your local machine, follow these steps:

1. **Clone the repository**

2. **Install dependencies using Poetry**:
    ```bash
    poetry install
    ```

3. **Run the application**:
- You can start MolDescriptor by executing the `run.py` script.

This will start the Flask application, and you can access it by navigating to `http://127.0.0.1:5000` in your web browser.

## The contents of the project

- ```/app```
This directory contains the Flask application and is the main entry point. In particular, ```app.py``` is the main file where the Flask app is created.

- ```/static```  
This directory is used to store static files, such as CSS files, JavaScript files, and images. These files don't change and are served as-is to the client.
For example, if you have a CSS file to style your web application, you would place it in the static directory.


- ```/templates```
This directory is used to store HTML templates. Flask uses a templating engine called Jinja2 that allows you to embed Python code within HTML files. This is how dynamic content is rendered.
For example, if you want to display the molecular descriptors in an HTML page, you'd have an HTML template in the templates directory that Flask fills with the appropriate data before sending it to the client.

## Documentation

Please visit the [documentation repo](https://github.com/moldescriptor/moldescriptor-docs) for documentation source material

- The ```config.py``` file is where you store configuration settings for your Flask app.

- ```requirements.txt``` is the file that holds the dependencies (packages) that are needed. This includes RDkit and Flask.

- ```run.py``` is the file that allows the whole app to run

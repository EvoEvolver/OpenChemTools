from typing import Dict, Any

from tools.toolset import toolset


@toolset.add()
def smiles_to_name(smiles: str) -> Dict[str, Any]:
    """
    Convert SMILES string to chemical names (IUPAC and common).
    :param smiles: SMILES string representation of the molecule
    :return: dictionary containing IUPAC name, common name, and molecular formula
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdMolDescriptors

        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES string"}

        # Get molecular formula
        formula = rdMolDescriptors.CalcMolFormula(mol)

        # Note: RDKit doesn't have built-in IUPAC naming
        # For production use, consider integrating with external services like:
        # - PubChem REST API
        # - ChemSpider API
        # - Or specialized naming libraries

        return {
            "smiles": smiles,
            "molecular_formula": formula,
            "molecular_weight": rdMolDescriptors.CalcExactMolWt(mol),
            "message": "SMILES parsed successfully. For IUPAC/common names, integration with external naming services is recommended."
        }
    except ImportError:
        return {"error": "RDKit not installed. Please install with: pip install rdkit"}
    except Exception as e:
        return {"error": f"SMILES conversion failed: {str(e)}"}


import requests

@toolset.add()
def name_to_smiles(name: str) -> Dict[str, Any]:
    """
    Convert chemical name to SMILES string using PubChem REST API.
    :param name: chemical name (IUPAC or common name)
    :return: dictionary containing SMILES string and molecular properties
    """
    try:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON"
        response = requests.get(url, timeout=10)
        response.raise_for_status()
        data = response.json()
        properties = data.get("PropertyTable", {}).get("Properties", [])
        if properties and "IsomericSMILES" in properties[0]:
            smiles = properties[0]["IsomericSMILES"]
            return {
                "input_name": name,
                "smiles": smiles,
                "message": "SMILES retrieved successfully from PubChem."
            }
        else:
            return {
                "input_name": name,
                "error": "SMILES not found for the given name in PubChem."
            }
    except requests.RequestException as e:
        return {"error": f"PubChem request failed: {str(e)}"}
    except Exception as e:
        return {"error": f"Name conversion failed: {str(e)}"}
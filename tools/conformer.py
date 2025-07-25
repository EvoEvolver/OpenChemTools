from typing import Dict, Any

from tools.toolset import toolset


@toolset.add("generate_3d_conformer")
def generate_3d_conformer(smiles: str, num_conformers: int = 1, optimize: bool = True) -> Dict[str, Any]:
    """
    Create 3D molecular structures from 2D SMILES representations.
    :param smiles: SMILES string of the molecule
    :param num_conformers: number of conformers to generate, default is 1
    :param optimize: whether to optimize the geometry using force field, default is True
    :return: dictionary containing 3D coordinates and conformer information
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, rdMolDescriptors
        import numpy as np

        # Parse SMILES and add hydrogens
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {"error": "Invalid SMILES string"}

        mol = Chem.AddHs(mol)

        # Generate 3D conformers
        conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, randomSeed=42)

        if len(conformer_ids) == 0:
            return {"error": "Failed to generate 3D conformers"}

        conformers_data = []

        for conf_id in conformer_ids:
            # Optimize geometry if requested
            if optimize:
                AllChem.MMFFOptimizeMolecule(mol, confId=conf_id)

            # Extract coordinates
            conf = mol.GetConformer(conf_id)
            atoms = []
            coordinates = []

            for atom in mol.GetAtoms():
                pos = conf.GetAtomPosition(atom.GetIdx())
                symbol = atom.GetSymbol()
                atoms.append(symbol)
                coordinates.append([pos.x, pos.y, pos.z])

            conformers_data.append({
                "conformer_id": int(conf_id),
                "atoms": atoms,
                "coordinates": coordinates,
                "energy": None  # Would need force field calculation for energy
            })

        # Try to convert to ASE Atoms object if ASE is available
        ase_objects = []
        try:
            from ase import Atoms
            for conf_data in conformers_data:
                ase_mol = Atoms(symbols=conf_data["atoms"], positions=conf_data["coordinates"])
                ase_objects.append({
                    "conformer_id": conf_data["conformer_id"],
                    "ase_atoms": "ASE Atoms object created successfully"
                })
        except ImportError:
            ase_objects = [{"message": "ASE not available. Install with: pip install ase"}]

        return {
            "smiles": smiles,
            "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
            "num_conformers_generated": len(conformers_data),
            "conformers": conformers_data,
            "ase_compatibility": ase_objects,
            "optimized": optimize
        }

    except ImportError as e:
        return {"error": f"Required package not installed: {str(e)}. Install with: pip install rdkit"}
    except Exception as e:
        return {"error": f"3D conformer generation failed: {str(e)}"}
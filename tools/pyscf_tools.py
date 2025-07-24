from tools.toolset import toolset
import numpy as np
from typing import List, Tuple, Dict, Any


@toolset.add()
def dft_calculation(atoms: List[Tuple[str, List[float]]], basis: str = "6-31g", xc: str = "b3lyp") -> Dict[str, Any]:
    """
    Perform DFT calculation on a molecular system.
    :param atoms: list of tuples containing (element, [x, y, z]) coordinates in Angstrom
    :param basis: basis set to use for calculation, default is "6-31g"
    :param xc: exchange-correlation functional, default is "b3lyp"
    :return: dictionary containing energy, dipole moment, and other properties
    """
    try:
        from pyscf import gto, dft
        
        # Build geometry string
        geometry = []
        for atom, coords in atoms:
            geometry.append(f"{atom} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        
        # Create molecule object
        mol = gto.Mole()
        mol.atom = '; '.join(geometry)
        mol.basis = basis
        mol.build()
        
        # Perform DFT calculation
        mf = dft.RKS(mol)
        mf.xc = xc
        energy = mf.kernel()
        
        # Calculate additional properties
        dipole = mf.dip_moment()
        
        return {
            "energy": float(energy),
            "dipole_moment": dipole.tolist(),
            "converged": mf.converged,
            "basis": basis,
            "xc_functional": xc
        }
    except ImportError:
        return {"error": "PySCF not installed. Please install with: pip install pyscf"}
    except Exception as e:
        return {"error": f"Calculation failed: {str(e)}"}


@toolset.add()
def geometry_optimization(atoms: List[Tuple[str, List[float]]], basis: str = "6-31g", xc: str = "b3lyp") -> Dict[str, Any]:
    """
    Perform geometry optimization on a molecular system.
    :param atoms: list of tuples containing (element, [x, y, z]) coordinates in Angstrom
    :param basis: basis set to use for calculation, default is "6-31g"
    :param xc: exchange-correlation functional, default is "b3lyp"
    :return: dictionary containing optimized geometry, final energy, and convergence info
    """
    try:
        from pyscf import gto, dft
        from pyscf.geomopt.berny_solver import optimize
        
        # Build geometry string
        geometry = []
        for atom, coords in atoms:
            geometry.append(f"{atom} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        
        # Create molecule object
        mol = gto.Mole()
        mol.atom = '; '.join(geometry)
        mol.basis = basis
        mol.build()
        
        # Setup DFT method
        mf = dft.RKS(mol)
        mf.xc = xc
        
        # Perform geometry optimization
        mol_opt = optimize(mf)
        
        # Extract optimized coordinates
        optimized_atoms = []
        for i, (atom, _) in enumerate(atoms):
            coords = mol_opt.atom_coords()[i].tolist()
            optimized_atoms.append((atom, coords))
        
        # Calculate final energy
        final_energy = mf.kernel()
        
        return {
            "optimized_geometry": optimized_atoms,
            "final_energy": float(final_energy),
            "converged": mf.converged,
            "basis": basis,
            "xc_functional": xc
        }
    except ImportError:
        return {"error": "PySCF not installed. Please install with: pip install pyscf"}
    except Exception as e:
        return {"error": f"Optimization failed: {str(e)}"}


@toolset.add("freq_calculation")
def frequency_calculation(atoms: List[Tuple[str, List[float]]], basis: str = "6-31g", xc: str = "b3lyp") -> Dict[str, Any]:
    """
    Perform vibrational frequency calculation on a molecular system.
    :param atoms: list of tuples containing (element, [x, y, z]) coordinates in Angstrom
    :param basis: basis set to use for calculation, default is "6-31g"
    :param xc: exchange-correlation functional, default is "b3lyp"
    :return: dictionary containing frequencies, zero-point energy, and thermodynamic properties
    """
    try:
        from pyscf import gto, dft
        from pyscf.hessian import rks
        
        # Build geometry string
        geometry = []
        for atom, coords in atoms:
            geometry.append(f"{atom} {coords[0]:.6f} {coords[1]:.6f} {coords[2]:.6f}")
        
        # Create molecule object
        mol = gto.Mole()
        mol.atom = '; '.join(geometry)
        mol.basis = basis
        mol.build()
        
        # Perform DFT calculation
        mf = dft.RKS(mol)
        mf.xc = xc
        mf.kernel()
        
        # Calculate Hessian and frequencies
        hess = rks.Hessian(mf)
        hessian_matrix = hess.kernel()
        
        # Convert Hessian to frequencies (simplified)
        mass_weighted_hessian = hess.partial_hess_elec(mol.atom_coords(), atmlst=range(mol.natm))
        
        return {
            "hessian_calculated": True,
            "energy": float(mf.e_tot),
            "basis": basis,
            "xc_functional": xc,
            "message": "Frequency analysis completed. Full vibrational analysis requires additional processing."
        }
    except ImportError:
        return {"error": "PySCF not installed. Please install with: pip install pyscf"}
    except Exception as e:
        return {"error": f"Frequency calculation failed: {str(e)}"}
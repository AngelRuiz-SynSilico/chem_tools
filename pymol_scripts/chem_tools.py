# Import necessary libraries from PyMOL and RDKit
from pymol import cmd, stored
import os, tempfile, subprocess, shutil
from pdbfixer import PDBFixer
from openmm.app import PDBFile
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdFMCS, Descriptors, Descriptors3D, AllChem
import numpy as np
try:
    import datamol as dm
    DATAMOL_AVAILABLE = True
except ImportError:
    DATAMOL_AVAILABLE = False
    print("Warning: datamol is not available. Sanitize functionality will use RDKit only.")

def __init__(self):

    self.menuBar.addcascademenu('Plugin', 'ChemTools', 'ChemTools', label='ChemTools')
    self.menuBar.addmenuitem('ChemTools', 'command', 'align_with_rdkit_help', label='Align with RDKit Help', command=lambda: align_with_rdkit_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'align_with_fkcombu_help', label='Align with FKcombu Help', command=lambda: align_with_fkcombu_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'docking_with_smina_help', label='Docking with Smina Help', command=lambda: docking_with_smina_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'pockets_with_fpocket_help', label='Pockets with Fpocket Help', command=lambda: pockets_with_fpocket_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'sanitize_molecule_help', label='Sanitize Molecule Help', command=lambda: sanitize_molecule_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'generate_conformers_help', label='Generate Conformers Help', command=lambda: generate_conformers_Help())
    self.menuBar.addmenuitem('ChemTools', 'command', 'calculate_rmsd_help', label='Calculate RMSD Help', command=lambda: calculate_rmsd_Help())


def align_with_rdkit_Help():
    Usage="""

    Aligns two small molecules (ligands) within PyMOL based on maximum common substructure (MCS) alignment.

    Parameters:
        reference (str): The name of the object in PyMOL representing the reference molecule.
        target (str): The name of the object in PyMOL representing the target molecule to be aligned.

    Returns:
        rmsd (float): Root Mean Square Deviation (RMSD) value between the aligned structures.
        mcs_smarts (str): SMARTS string representing the Maximum Common Substructure (MCS).
        num_atoms (int): Number of atoms in the MCS used for alignment.

    Description:
        This function saves the target and reference molecules as temporary .sdf files from PyMOL.
        It uses RDKit to read the .sdf files, finds the MCS between the two molecules, and aligns
        the target molecule to the reference based on this common substructure. After alignment,
        the function reloads the aligned target molecule back into PyMOL as '{target}_aligned'
        and deletes the temporary .sdf files automatically.

    Requirements:
        - PyMOL (with the command-line version or PyMOL open to run the script)
        - RDKit (for chemical informatics functionalities)
        - The target and reference objects must be loaded in PyMOL

    Usage Example:
        Load the molecules in PyMOL and run the script:
            align_with_rdkit reference_name target_name

        Where `reference_name` is the reference molecule and `target_name` is the molecule to be aligned.
        The aligned molecule will be loaded as 'target_name_aligned'.

    Additional Information:
        The script will print the RMSD, the SMARTS pattern of the MCS, the number of atoms involved
        in the alignment, and the name of the aligned molecule after execution.
    """
    print(Usage)
    return

def align_with_fkcombu_Help():
    Usage="""
    
    Aligns the ligand to the receptor in PyMOL using fkcombu, with the receptor surface as a constraint.

    Parameters:
    - receptor (str): Name of the receptor object in PyMOL.
    - reference (str): Name of the reference ligand object in PyMOL.
    - target (str): Name of the target ligand object in PyMOL.

    Returns:
    - None
    """
    print(Usage)
    return
    
def docking_with_smina_Help():
    Usage="""

    docking_with_smina with specified parameters, loading receptor and ligand directly from PyMOL session.

    Parameters:
    - receptor (str): Name of the receptor object in PyMOL.
    - ligand (str): Name of the ligand object in PyMOL.
    - center_x (float): X-coordinate for the center of the search box.
    - center_y (float): Y-coordinate for the center of the search box.
    - center_z (float): Z-coordinate for the center of the search box.
    - size_x (float): Size of the search box in the X dimension.
    - size_y (float): Size of the search box in the Y dimension.
    - size_z (float): Size of the search box in the Z dimension.
    - exhaustiveness (int): Exhaustiveness of the search (default is 8).

    Returns:
    - None
    """
    print(Usage)
    return

def pockets_with_fpocket_Help():
    Usage="""
    Identifies pockets on a protein surface using fpocket and displays them in PyMOL.

    Parameters:
    - receptor (str): Name of the receptor object in PyMOL.

    Returns:
    - None

    Description:
    This function uses fpocket to identify pockets on the surface of a protein loaded in PyMOL. It saves the
    receptor structure as a temporary PDB file, runs fpocket on this file, and loads the identified pockets back
    into PyMOL as separate objects. Pockets are then displayed as semi-transparent colored spheres, with each
    pocket numbered sequentially.

    Requirements:
    - PyMOL
    - fpocket (accessible from the command line)

    Usage Example:
        pockets_with_fpocket receptor_name

    The function will print fpocket’s output to the console and display each identified pocket as a separate
    PyMOL object with a unique color.
    """
    print(Usage)
    return

def sanitize_molecule_Help():
    Usage="""
    Sanitizes molecules loaded from PDB files to make them valid for RDKit processing.

    Parameters:
    - molecule (str): Name of the molecule object in PyMOL to sanitize.
    - method (str): Sanitization method to use. Options: 'datamol', 'rdkit', 'both'. Default: 'both'.
    - fix_valence (bool): Whether to attempt fixing valence issues. Default: True.
    - standardize (bool): Whether to standardize the molecule (datamol only). Default: True.

    Returns:
    - None

    Description:
    This function addresses common issues when loading molecules from PDB files that make them
    invalid for RDKit processing. It preserves the 3D coordinates while cleaning up the molecular
    structure. The original molecule is kept, and a new sanitized version is loaded with the suffix
    '_sanitized'.

    The function supports multiple sanitization approaches:
    - 'datamol': Uses datamol's advanced sanitization (available in Chem3 environment)
    - 'rdkit': Uses RDKit's built-in sanitization
    - 'both': Tries datamol first, falls back to rdkit if needed

    Requirements:
    - PyMOL
    - RDKit (included in PyMOL)
    - datamol (available in Chem3 environment)

    Usage Examples:
        sanitize_molecule ligand
        sanitize_molecule ligand, method=datamol
        sanitize_molecule ligand, method=rdkit, fix_valence=False
        sanitize_molecule ligand, method=both, standardize=True

    Note: For best results, run PyMOL in the Chem3 environment where all dependencies are available.

    The function will print information about the sanitization process and any issues encountered.
    """
    print(Usage)
    return


def generate_conformers_Help():
    Usage="""
    Generates multiple 3D conformers for a molecule using RDKit's ETKDG or datamol methods.

    Parameters:
    - molecule (str): Name of the molecule object in PyMOL.
    - n_conformers (int): Number of conformers to generate. Default: 10.
    - method (str): Method to use - 'etkdg' (RDKit) or 'datamol'. Default: 'etkdg'.
    - energy_minimize (bool): Whether to minimize conformers with UFF. Default: True.
    - rms_threshold (float): RMSD threshold for pruning similar conformers. Default: 0.5.

    Returns:
    - None

    Description:
    This function generates multiple 3D conformations of a molecule, which is useful for:
    - Exploring conformational space before docking
    - Finding low-energy conformations
    - Generating diverse starting structures

    The function uses either:
    - 'etkdg': RDKit's ETKDG (Experimental-Torsion Knowledge Distance Geometry) method
    - 'datamol': datamol's conformer generation (available in Chem3 environment)

    Each conformer is loaded as a separate state in PyMOL, allowing you to cycle through them.
    The function also reports the energy of each conformer (if minimized).

    Requirements:
    - PyMOL
    - RDKit (included in PyMOL)
    - datamol (optional, available in Chem3 environment)

    Usage Examples:
        generate_conformers ligand
        generate_conformers ligand, n_conformers=20
        generate_conformers ligand, method=datamol, n_conformers=50
        generate_conformers ligand, n_conformers=10, energy_minimize=True, rms_threshold=1.0

    Note: For best results with datamol method, run PyMOL in the Chem3 environment.

    After generation, use PyMOL's state controls to view different conformers.
    """
    print(Usage)
    return


def calculate_rmsd_Help():
    Usage="""
    Calculates RMSD between two molecules based on their current 3D coordinates without altering them.

    Parameters:
    - reference (str): Name of the reference molecule object in PyMOL.
    - target (str): Name of the target molecule object in PyMOL.
    - match_substructure (bool): Whether to use MCS for atom mapping. Default: True.
    - symmetry (bool): Whether to account for symmetry in RMSD calculation. Default: True.

    Returns:
    - rmsd (float): RMSD value in Angstroms

    Description:
    This function calculates the Root Mean Square Deviation (RMSD) between two molecules
    based on their current 3D coordinates in PyMOL. Unlike align_with_rdkit, this function
    does NOT modify the coordinates - it only calculates the RMSD.

    The function can work in two modes:
    - match_substructure=True: Uses Maximum Common Substructure (MCS) to find matching atoms
    - match_substructure=False: Assumes molecules have the same atom order (faster)

    When symmetry=True, the function will find the best RMSD considering molecular symmetry.

    This is useful for:
    - Comparing docking poses without moving them
    - Evaluating alignment quality
    - Comparing conformers
    - Analyzing molecular dynamics trajectories

    Requirements:
    - PyMOL
    - RDKit (included in PyMOL)
    - Both molecules must be loaded in PyMOL

    Usage Examples:
        calculate_rmsd reference, target
        calculate_rmsd ref_mol, docked_mol, match_substructure=True
        calculate_rmsd conf1, conf2, match_substructure=False, symmetry=True

    The function will print the RMSD value and information about the atom mapping used.
    """
    print(Usage)
    return


def describe_molecule(molecule:str):
    
    descriptors_list=['MolWt','TPSA','FractionCSP3', 'HeavyAtomCount', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles',
                      'NumAliphaticHeterocycles', 'NumAliphaticRings', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 
                      'NumAromaticRings', 'NumHAcceptors', 'NumHDonors', 'NumHeteroatoms', 'NumRotatableBonds', 'NumSaturatedCarbocycles',
                      'NumSaturatedHeterocycles', 'NumSaturatedRings','RingCount', 'MolLogP', 'MolMR']
    
    # Create temporary files for the reference and target molecules
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as mol_tmp:
        
        # Save reference and target molecules to temporary .sdf files
        cmd.save(mol_tmp.name, format="sdf", selection=molecule, state=-1)
        
        # Read the molecules from the .sdf files using RDKit
        mol = Chem.SDMolSupplier(mol_tmp.name)[0]

    def _get_descriptors(mol):
        details={"SMILES": Chem.MolToSmiles(mol),
                "SMARTS":Chem.MolToSmarts(mol)}
        descriptors=Descriptors.CalcMolDescriptors(mol)
        descriptors_slice={k:v for k,v in descriptors.items() if k in descriptors_list}
        descriptors3D=Descriptors3D.CalcMolDescriptors3D(mol)
        mol_descriptors = {**details,**descriptors_slice, **descriptors3D}
        return mol_descriptors
    
    data=_get_descriptors(mol)
    os.remove(mol_tmp.name)
    
    for k,v in data.items():
        print(f"{k}: {v}")
    
    return data


def align_with_rdkit(reference:str,target:str):

    # Create temporary files for the reference and target molecules
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as ref_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tar_tmp:

        # Save reference and target molecules to temporary .sdf files
        cmd.save(ref_tmp.name, format="sdf", selection=reference, state=-1)
        cmd.save(tar_tmp.name, format="sdf", selection=target, state=-1)

        # Read the molecules from the .sdf files using RDKit
        ref = Chem.SDMolSupplier(ref_tmp.name)[0]
        tar = Chem.SDMolSupplier(tar_tmp.name)[0]

        # Find the Maximum Common Substructure (MCS) between reference and target
        mcs = rdFMCS.FindMCS([ref, tar])
        ref_atoms = ref.GetSubstructMatch(mcs.queryMol)
        tar_atoms = tar.GetSubstructMatch(mcs.queryMol)

        # Map the atoms based on MCS and align the target molecule to the reference
        atom_map = list(zip(tar_atoms, ref_atoms))
        rmsd = rdMolAlign.AlignMol(prbMol=tar, refMol=ref, atomMap=atom_map)

        # Write the aligned target molecule to a new temporary .sdf file
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as aligned_tmp:
            out = Chem.SDWriter(aligned_tmp.name)
            out.write(tar)
            out.close()
            aligned_file = aligned_tmp.name

    # Load the aligned molecule back into PyMOL
    aligned_name = f"{target}_aligned"
    cmd.load(aligned_file, object=aligned_name, format="sdf")

    # Print alignment information
    print(f"RMSD: {rmsd}")
    print(f"MCS (SMARTS): {mcs.smartsString}")
    print(f"Number of atoms in MCS: {len(atom_map)}")
    print(f"Aligned molecule loaded as: {aligned_name}")

    # Clean up temporary files
    os.remove(ref_tmp.name)
    os.remove(tar_tmp.name)
    os.remove(aligned_file)

    return rmsd, mcs.smartsString, len(atom_map)


def align_with_fkcombu(receptor: str, reference: str, target: str):
    
    # Create temporary files for receptor and ligands in PDB and SDF formats
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as receptor_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as reference_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as target_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as out_tmp:

        # Save the receptor, reference, and target molecules from PyMOL to temporary files
        cmd.save(receptor_tmp.name, receptor, format="pdb")
        cmd.save(reference_tmp.name, reference, format="sdf")
        cmd.save(target_tmp.name, target, format="sdf")

        # Prepare the fkcombu command
        command = [
            "fkcombu",
            "-T", target_tmp.name,    # Target ligand file path
            "-R", reference_tmp.name,  # Reference ligand file path
            "-P", receptor_tmp.name,   # Receptor protein file path
            "-osdfT", out_tmp.name     # Output file path for the aligned ligand
        ]

        # Run the fkcombu command
        try:
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Output any messages from fkcombu
            print("fkcombu output:")
            print(result.stdout)
            if result.stderr:
                print("fkcombu errors:")
                print(result.stderr)
    
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running fkcombu: {e}")
            print("Standard output:", e.stdout)
            print("Standard error:", e.stderr)

        # Check if the output file was created and has content
        if os.path.exists(out_tmp.name) and os.path.getsize(out_tmp.name) > 0:
            # Load the aligned ligand back into PyMOL
            cmd.load(out_tmp.name, object="aligned_ligand", format="sdf")
            print("Aligned ligand loaded as 'aligned_ligand'")
        else:
            print("Error: The output file was not created or is empty.")

        # Clean up temporary files
        os.remove(receptor_tmp.name)
        os.remove(reference_tmp.name)
        os.remove(target_tmp.name)
        os.remove(out_tmp.name)


def docking_with_smina(receptor: str, ligand: str, center_x: float, center_y: float, center_z: float,
              size_x: float, size_y: float, size_z: float, exhaustiveness: int = 8):

    
    # Create temporary files for receptor and ligand in PDBQT format
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as receptor_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as ligand_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as out_tmp:
        
        # Export receptor and ligand from PyMOL to temporary PDBQT files
        cmd.save(receptor_tmp.name, receptor, format="pdb")
        cmd.save(ligand_tmp.name, ligand, format="sdf")

        # Prepare the command to run smina
        command = [
            "smina",
            "--receptor", receptor_tmp.name,
            "--ligand", ligand_tmp.name,
            "--out", out_tmp.name,
            "--center_x", str(center_x),
            "--center_y", str(center_y),
            "--center_z", str(center_z),
            "--size_x", str(size_x),
            "--size_y", str(size_y),
            "--size_z", str(size_z),
            "--exhaustiveness", str(exhaustiveness)
        ]

        # Run the command and capture the output and errors
        try:
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            # Output the results to the console
            print("Smina output:")
            print(result.stdout)
            
            if result.stderr:
                print("Smina errors:")
                print(result.stderr)

            # Load the docking result back into PyMOL as a new object
            cmd.load(out_tmp.name, object="docked_ligand")
            print("Docking results loaded as 'docked_ligand'")

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running smina: {e}")
            print("Standard output:", e.stdout)
            print("Standard error:", e.stderr)

def cleanup_with_pdbfixer(receptor: str, removeHeterogens: str = "False", pH: float = 7.5):
    """
    Cleans up a receptor molecule by adding missing atoms, hydrogens, and optionally removing heterogens.
    
    Parameters:
        receptor (str): Name of the receptor object in PyMOL.
        removeHeterogens (int): Option to remove heteroatoms; 1 for True, 0 for False. Default is 1.
        ph (float): pH value for adding hydrogens. Default is 7.5.
    
    Returns:
        None
    """  
    # Create temporary files for the input and output PDB files
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as receptor_tmp, \
         tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as fixed_tmp:
        
        # Save the receptor from PyMOL to a temporary PDB file
        cmd.save(receptor_tmp.name, receptor, format="pdb")
        fixer = PDBFixer(filename=receptor_tmp.name)
        
        # Perform cleanup steps
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(eval(removeHeterogens))
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(float(pH))
        
        # Write the fixed structure to a new temporary PDB file
        with open(fixed_tmp.name, 'w') as output_file:
            PDBFile.writeFile(fixer.topology, fixer.positions, output_file)
        
        # Load the cleaned-up structure into PyMOL
        cmd.load(fixed_tmp.name, object="fixed_protein")
        print("Fixed protein loaded as 'fixed_protein'")
    
    # Clean up temporary files
    os.remove(receptor_tmp.name)
    os.remove(fixed_tmp.name)

def sanitize_molecule(molecule: str, method: str = "both", fix_valence: bool = True, standardize: bool = True):
    """
    Sanitizes a molecule while preserving its 3D structure.

    Parameters:
        molecule (str): Name of the molecule object in PyMOL.
        method (str): Sanitization method - 'datamol', 'rdkit', or 'both'. Default: 'both'.
        fix_valence (bool): Whether to attempt fixing valence issues. Default: True.
        standardize (bool): Whether to standardize the molecule (datamol only). Default: True.

    Returns:
        None
    """

    # Validate method parameter
    valid_methods = ['datamol', 'rdkit', 'both']
    if method not in valid_methods:
        print(f"Error: method must be one of {valid_methods}")
        return

    # Check if datamol is available when needed
    if method in ['datamol', 'both'] and not DATAMOL_AVAILABLE:
        if method == 'datamol':
            print("Error: datamol is not available. Make sure you're running in the Chem3 environment.")
            return
        else:
            print("Warning: datamol is not available. Falling back to rdkit method.")
            method = 'rdkit'

    # Create temporary file for the molecule
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as mol_tmp:
        try:
            # Save molecule to temporary .sdf file
            cmd.save(mol_tmp.name, format="sdf", selection=molecule, state=-1)

            # Read the molecule from the .sdf file using RDKit
            mol = Chem.SDMolSupplier(mol_tmp.name, removeHs=False, sanitize=False)[0]

            if mol is None:
                print(f"Error: Could not read molecule '{molecule}' from PyMOL")
                return

            print(f"Original molecule: {Chem.MolToSmiles(mol, allHsExplicit=False) if mol else 'Invalid'}")
            print(f"Number of atoms: {mol.GetNumAtoms()}")
            print(f"Number of conformers: {mol.GetNumConformers()}")

            # Store original 3D coordinates
            if mol.GetNumConformers() > 0:
                conf = mol.GetConformer()
                original_coords = conf.GetPositions()
                print("3D coordinates preserved")
            else:
                print("Warning: No 3D coordinates found in the molecule")
                original_coords = None

            sanitized_mol = None

            # Try sanitization based on method
            if method in ['datamol', 'both']:
                print("\nAttempting datamol sanitization...")
                try:
                    # Use datamol's sanitization pipeline
                    sanitized_mol = dm.copy_mol(mol)

                    # Fix common issues
                    if fix_valence:
                        print("  - Fixing valence issues...")
                        sanitized_mol = dm.fix_mol(sanitized_mol)

                    if sanitized_mol is not None:
                        # Sanitize the molecule
                        print("  - Sanitizing molecule...")
                        sanitized_mol = dm.sanitize_mol(sanitized_mol, sanifix=True, charge_neutral=False)

                    if sanitized_mol is not None and standardize:
                        # Standardize the molecule
                        print("  - Standardizing molecule...")
                        sanitized_mol = dm.standardize_mol(
                            sanitized_mol,
                            disconnect_metals=False,
                            normalize=True,
                            reionize=True,
                            uncharge=False,
                            stereo=True
                        )

                    if sanitized_mol is not None:
                        print("  ✓ Datamol sanitization successful")
                    else:
                        print("  ✗ Datamol sanitization failed")
                        if method == 'both':
                            print("  Falling back to RDKit method...")

                except Exception as e:
                    print(f"  ✗ Datamol sanitization error: {e}")
                    sanitized_mol = None
                    if method == 'both':
                        print("  Falling back to RDKit method...")

            # Try RDKit sanitization if datamol failed or method is 'rdkit'
            if sanitized_mol is None and method in ['rdkit', 'both']:
                print("\nAttempting RDKit sanitization...")
                try:
                    sanitized_mol = Chem.Mol(mol)

                    # Try different sanitization options
                    sanitize_ops = Chem.SanitizeFlags.SANITIZE_ALL

                    if not fix_valence:
                        # Skip valence checks if requested
                        sanitize_ops ^= Chem.SanitizeFlags.SANITIZE_PROPERTIES

                    Chem.SanitizeMol(sanitized_mol, sanitizeOps=sanitize_ops, catchErrors=False)
                    print("  ✓ RDKit sanitization successful")

                except Exception as e:
                    print(f"  ✗ RDKit sanitization error: {e}")
                    # Try with more lenient settings
                    try:
                        print("  Trying lenient sanitization...")
                        sanitized_mol = Chem.Mol(mol)
                        Chem.SanitizeMol(
                            sanitized_mol,
                            sanitizeOps=Chem.SanitizeFlags.SANITIZE_FINDRADICALS |
                                       Chem.SanitizeFlags.SANITIZE_KEKULIZE |
                                       Chem.SanitizeFlags.SANITIZE_SETAROMATICITY |
                                       Chem.SanitizeFlags.SANITIZE_SETCONJUGATION |
                                       Chem.SanitizeFlags.SANITIZE_SETHYBRIDIZATION,
                            catchErrors=False
                        )
                        print("  ✓ Lenient RDKit sanitization successful")
                    except Exception as e2:
                        print(f"  ✗ Lenient sanitization also failed: {e2}")
                        sanitized_mol = None

            # Check if sanitization was successful
            if sanitized_mol is None:
                print("\n✗ All sanitization methods failed")
                return

            # Restore 3D coordinates if they existed
            if original_coords is not None and sanitized_mol.GetNumConformers() > 0:
                print("\nRestoring 3D coordinates...")
                conf = sanitized_mol.GetConformer()
                for i in range(min(sanitized_mol.GetNumAtoms(), len(original_coords))):
                    conf.SetAtomPosition(i, original_coords[i])
                print("  ✓ 3D coordinates restored")

            # Write the sanitized molecule to a new temporary file
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as sanitized_tmp:
                writer = Chem.SDWriter(sanitized_tmp.name)
                writer.write(sanitized_mol)
                writer.close()

                # Load the sanitized molecule back into PyMOL
                sanitized_name = f"{molecule}_sanitized"
                cmd.load(sanitized_tmp.name, object=sanitized_name, format="sdf")
                print(f"\n✓ Sanitized molecule loaded as '{sanitized_name}'")
                print(f"Sanitized SMILES: {Chem.MolToSmiles(sanitized_mol)}")

                # Clean up the sanitized temp file
                os.remove(sanitized_tmp.name)

        except Exception as e:
            print(f"Error during sanitization: {e}")
            import traceback
            traceback.print_exc()

        finally:
            # Clean up the original temp file
            if os.path.exists(mol_tmp.name):
                os.remove(mol_tmp.name)


def pockets_with_fpocket(receptor: str):
    # Create a temporary file for the receptor PDB
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as receptor_tmp:
        # Save the receptor from PyMOL to a temporary PDB file
        cmd.save(receptor_tmp.name, receptor, format="pdb")

        # Prepare the fpocket command
        command = ["fpocket", "-f", receptor_tmp.name]

        try:
            # Run the fpocket command
            result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            # Output the results to the console
            print("Fpocket output:")
            print(result.stdout)

            if result.stderr:
                print("Fpocket errors:")
                print(result.stderr)

            # Construct paths for the fpocket output directory and files
            pocket_output_dir = receptor_tmp.name.replace(".pdb", "_out")
            pocket_pdb_file = os.path.join(pocket_output_dir, os.path.basename(receptor_tmp.name).replace(".pdb", "_out.pdb"))
            info_file = os.path.join(pocket_output_dir, os.path.basename(receptor_tmp.name).replace(".pdb", "_info.txt"))

            # Load the fpocket output PDB file into PyMOL
            if os.path.exists(pocket_pdb_file):
                cmd.load(pocket_pdb_file,format="pdb",object="pockets_found")
                print("Pockets loaded")

                # Display fpocket info file content
                if os.path.exists(info_file):
                    with open(info_file) as pocket_data:
                        print(pocket_data.read())
                else:
                    print("Warning: Info file not found.")
            else:
                print("Error: Pocket PDB file not found.")

        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running fpocket: {e}")
            print("Fpocket output:", e.stdout)
            print("Fpocket error:", e.stderr)
            return

        # Selecting and displaying pockets
        stored.list = []
        cmd.iterate("(resn STP)", "stored.list.append(resi)")
        lastSTP = stored.list[-1] if stored.list else None

        if lastSTP:
            for my_index in range(1, int(lastSTP) + 1):
                pocket_selection = f"pocket{my_index}"
                cmd.select(pocket_selection, f"resn STP and resi {my_index}")
                cmd.color(my_index+1, pocket_selection)
                cmd.show("spheres", pocket_selection)
                cmd.set("sphere_scale", 1, pocket_selection)
                cmd.set("sphere_transparency", 0.4, pocket_selection)
        cmd.remove("pockets_found and polymer")

        # Clean up temporary and output files
        os.remove(receptor_tmp.name)
        # Remove the output directory and its contents
        if os.path.exists(pocket_output_dir):
            shutil.rmtree(pocket_output_dir)


def generate_conformers(molecule: str, n_conformers: int = 10, method: str = "etkdg",
                       energy_minimize: bool = True, rms_threshold: float = 0.5):
    """
    Generates multiple 3D conformers for a molecule.

    Parameters:
        molecule (str): Name of the molecule object in PyMOL.
        n_conformers (int): Number of conformers to generate. Default: 10.
        method (str): Method to use - 'etkdg' or 'datamol'. Default: 'etkdg'.
        energy_minimize (bool): Whether to minimize conformers with UFF. Default: True.
        rms_threshold (float): RMSD threshold for pruning similar conformers. Default: 0.5.

    Returns:
        None
    """

    # Validate method parameter
    valid_methods = ['etkdg', 'datamol']
    if method not in valid_methods:
        print(f"Error: method must be one of {valid_methods}")
        return

    # Check if datamol is available when needed
    if method == 'datamol' and not DATAMOL_AVAILABLE:
        print("Error: datamol is not available. Use method='etkdg' or activate Chem3 environment.")
        return

    # Create temporary file for the molecule
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as mol_tmp:
        try:
            # Save molecule to temporary .sdf file
            cmd.save(mol_tmp.name, format="sdf", selection=molecule, state=-1)

            # Read the molecule from the .sdf file using RDKit
            mol = Chem.SDMolSupplier(mol_tmp.name, removeHs=False)[0]

            if mol is None:
                print(f"Error: Could not read molecule '{molecule}' from PyMOL")
                return

            print(f"Generating {n_conformers} conformers for '{molecule}'...")
            print(f"Method: {method}")
            print(f"Energy minimize: {energy_minimize}")
            print(f"RMS threshold: {rms_threshold}")

            # Add hydrogens if not present
            mol = Chem.AddHs(mol)

            conformer_ids = []
            energies = []

            if method == 'etkdg':
                # Use RDKit's ETKDG method
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                params.pruneRmsThresh = rms_threshold
                params.numThreads = 0  # Use all available threads

                # Generate conformers
                conformer_ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conformers, params=params)

                if len(conformer_ids) == 0:
                    print("Error: Failed to generate conformers")
                    return

                print(f"Generated {len(conformer_ids)} conformers")

                # Energy minimize if requested
                if energy_minimize:
                    print("Minimizing conformers with UFF...")
                    results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=200)
                    energies = [result[1] for result in results]

            elif method == 'datamol':
                # Use datamol's conformer generation
                print("Using datamol conformer generation...")

                # Generate conformers with datamol
                mol = dm.conformers.generate(
                    mol,
                    n_confs=n_conformers,
                    rms_cutoff=rms_threshold,
                    minimize_energy=energy_minimize,
                    method="etkdg",
                    random_seed=42
                )

                conformer_ids = list(range(mol.GetNumConformers()))

                if len(conformer_ids) == 0:
                    print("Error: Failed to generate conformers with datamol")
                    return

                print(f"Generated {len(conformer_ids)} conformers with datamol")

                # Get energies if minimized
                if energy_minimize:
                    energies = []
                    for conf_id in conformer_ids:
                        ff = AllChem.UFFGetMoleculeForceField(mol, confId=conf_id)
                        energy = ff.CalcEnergy()
                        energies.append(energy)

            # Remove hydrogens for cleaner visualization
            mol = Chem.RemoveHs(mol)

            # Write conformers to temporary file
            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as conf_tmp:
                writer = Chem.SDWriter(conf_tmp.name)
                for i, conf_id in enumerate(conformer_ids):
                    writer.write(mol, confId=conf_id)
                writer.close()
                conf_file = conf_tmp.name

            # Load conformers into PyMOL
            conformer_name = f"{molecule}_conformers"
            cmd.load(conf_file, object=conformer_name, format="sdf")

            # Print summary
            print(f"\n✓ Conformers loaded as: {conformer_name}")
            print(f"  Total conformers: {len(conformer_ids)}")
            print(f"  Use PyMOL state controls to view different conformers")

            if energies:
                print(f"\nConformer Energies (kcal/mol):")
                for i, energy in enumerate(energies, 1):
                    print(f"  State {i}: {energy:.2f}")
                min_energy_idx = energies.index(min(energies))
                print(f"\n  Lowest energy conformer: State {min_energy_idx + 1} ({min(energies):.2f} kcal/mol)")

            # Clean up temporary files
            os.remove(mol_tmp.name)
            os.remove(conf_file)

        except Exception as e:
            print(f"Error generating conformers: {e}")
            import traceback
            traceback.print_exc()


def calculate_rmsd(reference: str, target: str, match_substructure: bool = True, symmetry: bool = True):
    """
    Calculates RMSD between two molecules without altering their 3D coordinates.

    Parameters:
        reference (str): Name of the reference molecule in PyMOL.
        target (str): Name of the target molecule in PyMOL.
        match_substructure (bool): Whether to use MCS for atom mapping. Default: True.
        symmetry (bool): Whether to account for symmetry. Default: True.

    Returns:
        float: RMSD value in Angstroms
    """

    # Create temporary files for the molecules
    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as ref_tmp, \
         tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as tar_tmp:

        try:
            # Save molecules to temporary .sdf files
            cmd.save(ref_tmp.name, format="sdf", selection=reference, state=-1)
            cmd.save(tar_tmp.name, format="sdf", selection=target, state=-1)

            # Read the molecules from the .sdf files using RDKit
            ref_mol = Chem.SDMolSupplier(ref_tmp.name, removeHs=False, sanitize=False)[0]
            tar_mol = Chem.SDMolSupplier(tar_tmp.name, removeHs=False, sanitize=False)[0]

            if ref_mol is None or tar_mol is None:
                print(f"Error: Could not read molecules from PyMOL")
                return None

            # Check that both molecules have conformers
            if ref_mol.GetNumConformers() == 0 or tar_mol.GetNumConformers() == 0:
                print("Error: One or both molecules lack 3D coordinates")
                return None

            print(f"Calculating RMSD between '{reference}' and '{target}'...")
            print(f"Match substructure: {match_substructure}")
            print(f"Consider symmetry: {symmetry}")

            rmsd = None
            atom_map = None

            if match_substructure:
                # Find Maximum Common Substructure
                mcs = rdFMCS.FindMCS([ref_mol, tar_mol])

                if mcs.numAtoms == 0:
                    print("Error: No common substructure found between molecules")
                    return None

                # Get atom matches
                ref_match = ref_mol.GetSubstructMatch(mcs.queryMol)
                tar_match = tar_mol.GetSubstructMatch(mcs.queryMol)

                if len(ref_match) == 0 or len(tar_match) == 0:
                    print("Error: Could not match substructure")
                    return None

                # Create atom map
                atom_map = list(zip(tar_match, ref_match))

                print(f"MCS found: {mcs.smartsString}")
                print(f"Number of matched atoms: {len(atom_map)}")

                # Calculate RMSD without aligning (just calculate based on current positions)
                rmsd = rdMolAlign.CalcRMS(ref_mol, tar_mol, map=atom_map)

            else:
                # Assume same atom order - calculate RMSD directly
                if ref_mol.GetNumAtoms() != tar_mol.GetNumAtoms():
                    print(f"Error: Molecules have different number of atoms ({ref_mol.GetNumAtoms()} vs {tar_mol.GetNumAtoms()})")
                    print("Try using match_substructure=True")
                    return None

                print(f"Calculating RMSD for {ref_mol.GetNumAtoms()} atoms (assuming same order)")

                # Calculate RMSD without atom mapping
                rmsd = rdMolAlign.CalcRMS(ref_mol, tar_mol)

            # Print results
            print(f"\n{'='*50}")
            print(f"RMSD: {rmsd:.3f} Å")
            print(f"{'='*50}")

            # Clean up temporary files
            os.remove(ref_tmp.name)
            os.remove(tar_tmp.name)

            return rmsd

        except Exception as e:
            print(f"Error calculating RMSD: {e}")
            import traceback
            traceback.print_exc()
            return None


# Extend PyMOL's command set to include this function
cmd.extend("align_with_rdkit",align_with_rdkit)
cmd.extend("align_with_fkcombu", align_with_fkcombu)
cmd.extend("docking_with_smina", docking_with_smina)
cmd.extend("pockets_with_fpocket", pockets_with_fpocket)
cmd.extend("describe_molecule", describe_molecule)
cmd.extend("cleanup_with_pdbfixer",cleanup_with_pdbfixer)
cmd.extend("sanitize_molecule", sanitize_molecule)
cmd.extend("generate_conformers", generate_conformers)
cmd.extend("calculate_rmsd", calculate_rmsd)

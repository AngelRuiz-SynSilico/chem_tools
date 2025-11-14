# Changes to chem_tools.py

## Summary

Added a new `sanitize_molecule` function to clean up molecules loaded from PDB files that may be invalid for RDKit processing, while preserving their 3D structures.

## New Functionality

### Function: `sanitize_molecule`

**Purpose**: Sanitize molecules from PDB files to make them valid for RDKit while preserving 3D coordinates.

**Signature**:
```python
sanitize_molecule(molecule: str, method: str = "both", fix_valence: bool = True, standardize: bool = True)
```

**Parameters**:
- `molecule` (str): Name of the molecule object in PyMOL
- `method` (str): Sanitization method - 'datamol', 'rdkit', or 'both' (default: 'both')
- `fix_valence` (bool): Whether to attempt fixing valence issues (default: True)
- `standardize` (bool): Whether to standardize the molecule using datamol (default: True)

**Usage Examples**:
```python
# Basic usage
sanitize_molecule ligand

# Use specific method
sanitize_molecule ligand, method=datamol
sanitize_molecule ligand, method=rdkit

# Custom options
sanitize_molecule ligand, method=both, fix_valence=True, standardize=True
```

**Output**: Creates a new PyMOL object with suffix `_sanitized` (e.g., `ligand_sanitized`)

## Implementation Details

### Dependencies Added
- **datamol** (optional): For advanced sanitization features
  - Install with: `pip install datamol`
  - If not installed, falls back to RDKit-only sanitization

### Code Changes

1. **Import Section** (lines 1-13):
   - Added conditional import for datamol
   - Added `DATAMOL_AVAILABLE` flag to check if datamol is installed

2. **Menu Integration** (line 22):
   - Added menu item for sanitize_molecule help

3. **Help Function** (lines 129-169):
   - Added `sanitize_molecule_Help()` function with comprehensive documentation

4. **Main Function** (lines 388-586):
   - Added `sanitize_molecule()` function with:
     - Parameter validation
     - Datamol sanitization pipeline
     - RDKit sanitization pipeline
     - 3D coordinate preservation
     - Error handling and fallback mechanisms
     - Detailed console output

5. **Command Extension** (line 637):
   - Registered `sanitize_molecule` as a PyMOL command

## Features

### ✓ Multiple Sanitization Methods
- **Datamol**: Advanced sanitization with fix_mol, sanitize_mol, and standardize_mol
- **RDKit**: Built-in RDKit sanitization with configurable flags
- **Both**: Tries datamol first, falls back to RDKit if needed

### ✓ 3D Coordinate Preservation
- Extracts original 3D coordinates before sanitization
- Restores coordinates after sanitization
- Maintains spatial information

### ✓ Comprehensive Error Handling
- Validates input parameters
- Checks for datamol availability
- Provides fallback mechanisms
- Detailed error messages

### ✓ Detailed Output
- Shows original molecule information
- Reports sanitization progress
- Indicates success/failure for each step
- Displays final SMILES

### ✓ Non-Destructive
- Original molecule remains unchanged
- Creates new sanitized version with `_sanitized` suffix

## What Gets Fixed

The sanitization process addresses:

1. **Valence Issues**: Incorrect valence states on atoms
2. **Bond Orders**: Missing or incorrect bond order assignments
3. **Aromaticity**: Improper aromaticity detection
4. **Charges**: Incorrect formal charges
5. **Stereochemistry**: Stereochemical inconsistencies
6. **Functional Groups**: Non-standard functional group representations
7. **Metal Coordination**: Metal-ligand bond issues

## Datamol Sanitization Pipeline

When using datamol method:

1. **Copy Molecule**: Creates a copy to preserve original
2. **Fix Valence** (if enabled): Uses `dm.fix_mol()` to fix valence issues
3. **Sanitize**: Uses `dm.sanitize_mol()` with sanifix enabled
4. **Standardize** (if enabled): Uses `dm.standardize_mol()` with:
   - Normalization of functional groups
   - Reionization for correct charge states
   - Stereochemistry correction
5. **Restore Coordinates**: Applies original 3D coordinates

## RDKit Sanitization Pipeline

When using rdkit method:

1. **Copy Molecule**: Creates a copy to preserve original
2. **Sanitize**: Applies RDKit sanitization flags
3. **Fallback**: If strict sanitization fails, tries lenient mode with:
   - SANITIZE_FINDRADICALS
   - SANITIZE_KEKULIZE
   - SANITIZE_SETAROMATICITY
   - SANITIZE_SETCONJUGATION
   - SANITIZE_SETHYBRIDIZATION
4. **Restore Coordinates**: Applies original 3D coordinates

## Files Added

1. **SANITIZE_README.md**: Comprehensive documentation
2. **test_sanitize.py**: Test script and usage examples
3. **example_sanitize_usage.pml**: PyMOL script with examples
4. **CHANGES.md**: This file

## Backward Compatibility

- All existing functions remain unchanged
- No breaking changes to existing API
- New functionality is purely additive

## Testing Recommendations

1. Test with molecules from PDB files that previously failed RDKit validation
2. Verify 3D coordinates are preserved
3. Test with both datamol and rdkit methods
4. Test fallback mechanism when datamol is not installed
5. Verify integration with existing functions (describe_molecule, align_with_rdkit)

## Future Enhancements

Potential improvements for future versions:

1. Batch sanitization of multiple molecules
2. Custom sanitization rules
3. Option to save sanitization report
4. Integration with molecular validation tools
5. Support for additional file formats

## References

- [Datamol Documentation](https://docs.datamol.io/)
- [Datamol GitHub](https://github.com/datamol-io/datamol)
- [RDKit Documentation](https://www.rdkit.org/docs/)
- [PyMOL Wiki](https://pymolwiki.org/)


# Interactive Molecular Structure Alignment Application

A PyQt5-based interactive tool for aligning and visualizing molecular analogs with customizable highlighting, transformations, and high-resolution export capabilities.

## Table of Contents

* [Installation](#installation)
* [Quick Start](#quick-start)
* [Features](#features)
* [User Interface](#user-interface)
* [Workflows](#workflows)
* [File Formats](#file-formats)
* [Troubleshooting](#troubleshooting)
* [Tips and Best Practices](#tips-and-best-practices)

## Installation

### Requirements

* Python 3.7 or higher
* RDKit (chemistry toolkit)
* PyQt5 (GUI framework)
* Pillow (image processing)

### Step-by-Step Installation

#### Using Conda (Recommended)
```bash
# Create a new conda environment
conda create -n mol_aligner python=3.9

# Activate the environment
conda activate mol_aligner

# Install RDKit from conda-forge
conda install -c conda-forge rdkit=2023.03.1

# Install PyQt5 and Pillow
conda install pyqt pillow
# Or using pip
pip install rdkit>=2023.03.1 PyQt5 pillow
```

**Note:** Installing RDKit via pip can be challenging on some systems. Conda installation is strongly recommended.

### Platform-Specific Notes

#### macOS
* No additional steps required
* Font rendering uses system fonts (Helvetica, Arial)

#### Linux
* Install font packages for better rendering:
```bash
sudo apt-get install fonts-dejavu-core  # Ubuntu/Debian
sudo yum install dejavu-sans-fonts      # CentOS/RHEL
```

#### Windows
* No additional steps required
* Arial font is used by default

## Quick Start

### Launching the Application
```bash
python align-2d-struct-highlight_interactive.py
```

### Basic Workflow

1. **Load molecules:** Click "Browse File" and select your SDF or SMILES file
2. **Set reference:** Enter the SMILES string of your reference/core structure
3. **Configure display:** Set number of columns (default: 4) and font size (default: 24)
4. **Click "Load & Align":** Molecules will be aligned and displayed
5. **Edit as needed:** Select molecules, toggle highlights, rotate/mirror structures
6. **Export:** Click "💾 Export High-Res PNG" to save your grid

## Features

### Core Features

* **Molecular Alignment:** Automatically aligns analogs based on Maximum Common Substructure (MCS)
* **Interactive Highlighting:** Click atoms to toggle highlights on/off
* **Structure Transformations:** Rotate (90°) and mirror structures while keeping atom labels readable
* **Flexible Display:** Adjustable grid columns (1-10) and font sizes (16-48pt)
* **Property Display:** Show any molecular property as label (Catalog ID, Price, MW, LogP, etc.)
* **High-Resolution Export:** Export publication-quality PNG at 300 DPI

### Supported File Formats

#### Input Files

* **SDF (Structure Data File):** .sdf extension
* **SMILES:** .smi or .smiles extension

#### SMILES File Format

Tab or space-delimited format:
```text
SMILES_STRING    Catalog_ID    Property1    Property2
c1ccccc1         Benzene       10.50        150.2
CCO              Ethanol       5.00         78.1
```

* First column: SMILES string (required)
* Second column: Catalog ID (optional, defaults to "Mol_N")
* Additional columns: Any properties you want to display

#### SDF Files

Standard SDF format with properties embedded:
* Catalog ID
* Price, USD
* Molecular Weight
* LogP
* Any custom properties

## User Interface

### Main Window Components

#### 1. Toolbar

| Button | Function | Shortcut |
|--------|----------|----------|
| **✏️ Edit Highlights Mode** | Toggle between selection and atom editing modes | Click to toggle |
| **📝 Atom List...** | Open dialog with checkboxes for all atoms | - |
| **Reset Highlights** | Restore original highlights (non-matching atoms) | - |
| **Clear All Highlights** | Remove all highlights | - |
| **Rotate ↻ 90°** | Rotate selected molecule 90° clockwise | - |
| **Rotate ↺ 90°** | Rotate selected molecule 90° counter-clockwise | - |
| **Mirror Horizontal ↔** | Flip selected molecule left-to-right | - |
| **Mirror Vertical ↕** | Flip selected molecule top-to-bottom | - |
| **💾 Export High-Res PNG** | Save grid as high-resolution image | - |
| **Reset Selected** | Reset all transformations for selected molecule | - |

#### 2. Input Section

* **File:** Path to SDF or SMILES file
* **Browse File:** Open file dialog
* **Reference:** SMILES string of reference structure
* **Columns:** Number of columns in grid (1-10)
* **Load & Align:** Process and display molecules

#### 3. Label Settings

* **Font Size:** Dropdown with sizes 16-48pt
* **Second Line:** Dropdown to select which property to display
  * "Catalog ID only" - single line label
  * Any property from your file

#### 4. Grid Display

* Scrollable area showing aligned molecules
* Click molecule to select (blue border)
* In Edit Mode: Click atoms to toggle highlights

#### 5. Status Bar

* Shows current operation status
* Displays MCS statistics after loading
* Shows atom selection feedback

### Visual Feedback

* **Selection:** Blue border (3px) around selected molecule
* **Default:** Gray border (1px) around unselected molecules
* **Edit Mode:** Info banner turns blue, cursor becomes crosshair
* **Selection Mode:** Gray info banner, cursor is pointing hand
* **Highlights:** Light green circles around highlighted atoms

## Workflows

### Workflow 1: Basic Structure Alignment

**Goal:** Align a series of molecular analogs and highlight variable regions

1. Prepare your data:
   * SDF file with multiple analogs
   * Reference SMILES (the core scaffold)
2. Launch the application
3. Click "Browse File" and select your SDF
4. Enter reference SMILES in "Reference" field
5. Set "Columns" to desired value (e.g., 4 for 4×N grid)
6. Click "Load & Align"
7. Review alignment:
   * Status bar shows MCS statistics
   * Green highlights show atoms different from reference
8. Export: Click "💾 Export High-Res PNG"

### Workflow 2: Custom Highlight Editing

**Goal:** Manually adjust which atoms are highlighted

**Method 1: Direct Atom Clicking**
1. Load and align molecules (see Workflow 1)
2. Click "✏️ Edit Highlights Mode" button in toolbar
3. Click directly on atoms to toggle highlights on/off
4. Status bar confirms each atom toggle
5. Click "✏️ Edit Highlights Mode" again to exit

**Method 2: Atom List Dialog**
1. Load and align molecules
2. Click a molecule to select it
3. Click "📝 Atom List..." in toolbar
4. Check/uncheck atoms in the list
5. Use "Select All", "Select None", or "Invert" buttons
6. Click "OK" to apply

**Method 3: Bulk Operations**
1. Select a molecule
2. Click "Reset Highlights" to restore original highlights
3. Click "Clear All Highlights" to remove all highlights

### Workflow 3: Structure Transformations

**Goal:** Adjust molecule orientation for better visual comparison

1. Load and align molecules
2. Click a molecule to select it (blue border appears)
3. Apply transformations:
   * **Rotate:** Click "Rotate ↻ 90°" or "Rotate ↺ 90°"
   * **Mirror:** Click "Mirror Horizontal ↔" or "Mirror Vertical ↕"
4. Transformations can be combined:
   * Example: Rotate 90° + Mirror Horizontal
5. **Important:** Atom labels remain upright and readable
6. To undo: Click "Reset Selected"

### Workflow 4: Customizing Labels

**Goal:** Display specific molecular properties

1. Load molecules with properties
2. After loading, check "Second Line" dropdown
   * It auto-populates with available properties
3. Select desired property:
   * "Catalog ID only" - single line
   * "Price, USD" - shows price
   * "Molecular Weight" - shows MW
   * Any other property from your file
4. Adjust "Font Size" for readability (try 32 or 36)
5. Changes apply immediately to all molecules
6. Export preserves your settings

### Workflow 5: Publication-Quality Export

**Goal:** Create high-resolution image for papers/presentations

1. Complete your analysis (alignment, highlights, labels)
2. Configure display:
   * Font Size: 32-40pt for publication
   * Second Line: Choose relevant property
   * Columns: Adjust for best layout
3. Click "💾 Export High-Res PNG"
4. Choose save location
5. Result:
   * 700×700 pixels per molecule (2× display size)
   * 300 DPI resolution
   * Fonts scaled 2× for clarity
   * All transformations and highlights preserved

## File Formats

### Creating SMILES Files

**Basic format:**
```text
c1ccccc1         Benzene      10.50
c1ccc(O)cc1      Phenol       15.20
c1ccc(N)cc1      Aniline      18.00
```

**With multiple properties:**
```text
SMILES           ID           Price    MW      LogP
c1ccccc1         BEN-001      10.50    78.11   2.13
c1ccc(O)cc1      PHE-002      15.20    94.11   1.46
c1ccc(N)cc1      ANI-003      18.00    93.13   0.90
```

**Naming convention:**
* Use .smi or .smiles extension
* Tab or space-delimited
* First column must be valid SMILES
* Comments: Lines starting with # are ignored

### SDF File Properties

The application reads all properties from SDF files. Common properties:

* **Catalog ID** - Molecule identifier
* **Price, USD** - Pricing information
* **Molecular Weight** - Calculated or experimental MW
* **LogP** - Partition coefficient
* **SMILES** - SMILES representation
* **Custom properties** - Any field in your SDF

**Tip:** Use RDKit to add custom properties to SDF:
```python
from rdkit import Chem

mol = Chem.MolFromSmiles('c1ccccc1')
mol.SetProp('Catalog ID', 'BEN-001')
mol.SetProp('Price, USD', '10.50')
mol.SetProp('Custom Property', 'Value')

writer = Chem.SDWriter('output.sdf')
writer.write(mol)
writer.close()
```

## Troubleshooting

### Performance Issues

#### Slow loading with many molecules

* MCS calculation is O(n) for n molecules
* For >100 molecules, expect 10-30 second load time
* Timeout is 10 seconds per molecule for MCS

## License
(c) Claude and Andrii

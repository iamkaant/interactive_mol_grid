"""
Interactive Molecular Structure Alignment and Highlighting Application
Requires: rdkit, PyQt5, pillow
Install: pip install rdkit PyQt5 pillow
"""

from rdkit import Chem, Geometry
from rdkit.Chem import AllChem, Draw, rdFMCS
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import rdDepictor
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QScrollArea, QLabel, QPushButton,
                             QFileDialog, QMessageBox, QLineEdit, QSpinBox,
                             QToolBar, QAction, QFrame, QGridLayout, QGroupBox,
                             QDialog, QListWidget, QListWidgetItem, QCheckBox,
                             QDialogButtonBox, QComboBox)
from PyQt5.QtCore import Qt, QPoint, QRect, QPointF
from PyQt5.QtGui import QPixmap, QImage, QPainter, QPen, QColor, QCursor
from io import BytesIO
from PIL import Image
import numpy as np
import math
import sys


class MoleculeData:
    """Stores molecule and its rendering state"""
    def __init__(self, mol, match_atoms, all_atoms, catalog_id=""):
        self.mol = mol
        self.original_mol = Chem.Mol(mol)  # Deep copy
        self.match_atoms = match_atoms  # Atoms matching reference
        self.all_atoms = all_atoms  # All atom indices
        self.highlighted_atoms = set(all_atoms) - set(match_atoms)  # Non-matching atoms
        self.catalog_id = catalog_id
        self.rotation = 0  # 0, 90, 180, 270
        self.mirror_h = False
        self.mirror_v = False

    def reset_highlights(self):
        """Reset to original highlights"""
        self.highlighted_atoms = set(self.all_atoms) - set(self.match_atoms)

    def toggle_highlight(self, atom_idx):
        """Toggle highlight for an atom"""
        if atom_idx in self.highlighted_atoms:
            self.highlighted_atoms.remove(atom_idx)
        else:
            self.highlighted_atoms.add(atom_idx)

    def get_transformed_mol(self):
        """Get molecule with transformed 2D coordinates"""
        # Create a copy of the molecule
        mol_copy = Chem.Mol(self.original_mol)

        if not mol_copy.GetNumConformers():
            return mol_copy

        conf = mol_copy.GetConformer()

        # Apply transformations to coordinates
        for atom_idx in range(mol_copy.GetNumAtoms()):
            pos = conf.GetAtomPosition(atom_idx)
            x, y = pos.x, pos.y

            # Apply rotation
            if self.rotation == 90:
                x, y = y, -x
            elif self.rotation == 180:
                x, y = -x, -y
            elif self.rotation == 270:
                x, y = -y, x

            # Apply mirroring
            if self.mirror_h:
                x = -x
            if self.mirror_v:
                y = -y

            conf.SetAtomPosition(atom_idx, Geometry.Point3D(x, y, 0))

        return mol_copy


class AtomSelectorDialog(QDialog):
    """Dialog for selecting atoms to highlight"""
    def __init__(self, mol_data, parent=None):
        super().__init__(parent)
        self.mol_data = mol_data
        self.setWindowTitle("Select Atoms to Highlight")
        self.setModal(True)
        self.resize(400, 500)

        layout = QVBoxLayout(self)

        # Instructions
        label = QLabel("Check atoms to highlight:")
        layout.addWidget(label)

        # Atom list with checkboxes
        self.list_widget = QListWidget()
        mol = mol_data.mol

        for atom_idx in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(atom_idx)
            symbol = atom.GetSymbol()
            item_text = f"Atom {atom_idx}: {symbol}"

            item = QListWidgetItem(item_text)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)

            # Set check state based on current highlights
            if atom_idx in mol_data.highlighted_atoms:
                item.setCheckState(Qt.Checked)
            else:
                item.setCheckState(Qt.Unchecked)

            self.list_widget.addItem(item)

        layout.addWidget(self.list_widget)

        # Buttons
        btn_layout = QHBoxLayout()

        btn_all = QPushButton("Select All")
        btn_all.clicked.connect(self.select_all)
        btn_layout.addWidget(btn_all)

        btn_none = QPushButton("Select None")
        btn_none.clicked.connect(self.select_none)
        btn_layout.addWidget(btn_none)

        btn_invert = QPushButton("Invert")
        btn_invert.clicked.connect(self.invert_selection)
        btn_layout.addWidget(btn_invert)

        layout.addLayout(btn_layout)

        # OK/Cancel buttons
        button_box = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        button_box.accepted.connect(self.accept)
        button_box.rejected.connect(self.reject)
        layout.addWidget(button_box)

    def select_all(self):
        """Select all atoms"""
        for i in range(self.list_widget.count()):
            self.list_widget.item(i).setCheckState(Qt.Checked)

    def select_none(self):
        """Deselect all atoms"""
        for i in range(self.list_widget.count()):
            self.list_widget.item(i).setCheckState(Qt.Unchecked)

    def invert_selection(self):
        """Invert selection"""
        for i in range(self.list_widget.count()):
            item = self.list_widget.item(i)
            if item.checkState() == Qt.Checked:
                item.setCheckState(Qt.Unchecked)
            else:
                item.setCheckState(Qt.Checked)

    def get_selected_atoms(self):
        """Get set of selected atom indices"""
        selected = set()
        for i in range(self.list_widget.count()):
            if self.list_widget.item(i).checkState() == Qt.Checked:
                selected.add(i)
        return selected


class MoleculeWidget(QLabel):
    """Interactive widget for displaying and manipulating molecule"""
    def __init__(self, mol_data, img_size=(400, 400), main_window=None, parent=None):
        super().__init__(parent)
        self.mol_data = mol_data
        self.img_size = img_size
        self.main_window = main_window  # Reference to main window
        self.selected = False
        self.atom_positions = {}  # Maps atom_idx to screen QPointF

        self.setFrameStyle(QFrame.Box | QFrame.Plain)
        self.setLineWidth(2)
        self.update_display()
        self.update_cursor()

    def update_display(self):
        """Render the molecule with current settings"""
        # Get transformed molecule (with proper coordinate transformations)
        mol = self.mol_data.get_transformed_mol()

        # Render molecule with RDKit (atom labels stay upright)
        d2d = rdMolDraw2D.MolDraw2DCairo(self.img_size[0], self.img_size[1])

        # Set drawing options for thicker bonds
        draw_options = d2d.drawOptions()
        draw_options.bondLineWidth = 3.0  # Increased from default 2.0

        # Prepare highlighting
        highlight_atoms = list(self.mol_data.highlighted_atoms)
        highlight_colors = {idx: (0.596, 0.984, 0.596) for idx in highlight_atoms}  # Light green
        highlight_radii = {idx: 0.8 for idx in highlight_atoms}

        # Draw molecule
        d2d.DrawMolecule(mol,
                        highlightAtoms=highlight_atoms,
                        highlightAtomColors=highlight_colors,
                        highlightAtomRadii=highlight_radii)
        d2d.FinishDrawing()

        # Store atom positions for click detection
        self.atom_positions = {}
        if mol.GetNumConformers() > 0:
            # Get molecule bounds to calculate transformation
            conf = mol.GetConformer()
            min_x = min_y = float('inf')
            max_x = max_y = float('-inf')

            for atom_idx in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(atom_idx)
                min_x = min(min_x, pos.x)
                max_x = max(max_x, pos.x)
                min_y = min(min_y, pos.y)
                max_y = max(max_y, pos.y)

            # Calculate scale and offset (RDKit centers and scales to fit)
            mol_width = max_x - min_x
            mol_height = max_y - min_y

            # Add padding (RDKit uses ~10% padding)
            padding = 0.1
            scale_x = self.img_size[0] * (1 - 2*padding) / mol_width if mol_width > 0 else 1
            scale_y = self.img_size[1] * (1 - 2*padding) / mol_height if mol_height > 0 else 1
            scale = min(scale_x, scale_y)

            center_x = (min_x + max_x) / 2
            center_y = (min_y + max_y) / 2

            # Calculate screen positions
            for atom_idx in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(atom_idx)
                # Transform to screen coordinates
                # Note: Invert Y because RDKit Y-axis goes up, screen Y-axis goes down
                screen_x = (pos.x - center_x) * scale + self.img_size[0] / 2
                screen_y = -(pos.y - center_y) * scale + self.img_size[1] / 2
                self.atom_positions[atom_idx] = QPointF(screen_x, screen_y)

        # Get image
        mol_img_bytes = d2d.GetDrawingText()
        bio = BytesIO(mol_img_bytes)
        pil_img = Image.open(bio).convert('RGB')

        # Add text labels
        pil_img = self._add_text_labels(pil_img)

        # Convert to QPixmap
        img_bytes = BytesIO()
        pil_img.save(img_bytes, format='PNG')
        img_bytes.seek(0)

        qimage = QImage.fromData(img_bytes.read())
        pixmap = QPixmap.fromImage(qimage)

        self.setPixmap(pixmap)
        self.setFixedSize(pixmap.size())

        # Update border color based on selection
        if self.selected:
            self.setStyleSheet("border: 3px solid #2196F3;")
        else:
            self.setStyleSheet("border: 1px solid #CCCCCC;")

    def _add_text_labels(self, img):
        """Add catalog ID and price labels"""
        from PIL import ImageDraw, ImageFont

        # Get font size from main window
        if self.main_window:
            font_size = int(self.main_window.font_size_combo.currentText())
        else:
            font_size = 24

        # Try to load a nice font - try multiple common fonts
        font = None
        font_names = ["arial.ttf", "Arial.ttf", "DejaVuSans.ttf", "Helvetica.ttf"]
        for font_name in font_names:
            try:
                font = ImageFont.truetype(font_name, font_size)
                break
            except:
                continue

        # Fallback - create a basic font at larger size
        if font is None:
            try:
                font = ImageFont.load_default()
                # For default font, we need to manually scale
                # Create a scaled version by adjusting image
            except:
                font = ImageFont.load_default()

        line1 = self.mol_data.catalog_id

        # Get second line based on property selection
        line2 = None
        if self.main_window:
            selected_prop = self.main_window.property_combo.currentData()
            if selected_prop != "catalog_only":
                # Get the property value
                if self.mol_data.mol.HasProp(selected_prop):
                    prop_value = self.mol_data.mol.GetProp(selected_prop)
                    if selected_prop == "Price, USD":
                        line2 = f'Price: ${prop_value}'
                    else:
                        line2 = f'{selected_prop}: {prop_value}'
                else:
                    line2 = None

        # Calculate dimensions with the properly sized font
        draw = ImageDraw.Draw(img)
        bbox1 = draw.textbbox((0, 0), line1, font=font)
        w1, h1 = bbox1[2] - bbox1[0], bbox1[3] - bbox1[1]

        if line2:
            bbox2 = draw.textbbox((0, 0), line2, font=font)
            w2, h2 = bbox2[2] - bbox2[0], bbox2[3] - bbox2[1]
            text_height = h1 + h2 + 20
        else:
            text_height = h1 + 10

        # Create new image with space for text
        new_img = Image.new('RGB', (img.width, img.height + text_height), 'white')
        new_img.paste(img, (0, 0))

        # Draw text
        draw = ImageDraw.Draw(new_img)
        x1 = (img.width - w1) // 2
        y1 = img.height + 5
        draw.text((x1, y1), line1, font=font, fill='black')

        if line2:
            x2 = (img.width - w2) // 2
            y2 = y1 + h1 + 5
            draw.text((x2, y2), line2, font=font, fill='black')

        return new_img

    def mousePressEvent(self, event):
        """Handle click to select molecule or toggle atom highlight"""
        if event.button() == Qt.LeftButton:
            # Check if we're in edit highlights mode
            if self.main_window and self.main_window.edit_highlights_mode:
                # Find clicked atom and toggle highlight
                click_pos = QPointF(event.pos())
                atom_idx = self._find_atom_at_position(click_pos)

                if atom_idx is not None:
                    self.mol_data.toggle_highlight(atom_idx)
                    self.update_display()
                    if self.main_window:
                        highlight_status = "highlighted" if atom_idx in self.mol_data.highlighted_atoms else "unhighlighted"
                        self.main_window.statusBar().showMessage(
                            f"Atom {atom_idx} ({self.mol_data.mol.GetAtomWithIdx(atom_idx).GetSymbol()}) {highlight_status} in {self.mol_data.catalog_id}"
                        )
                else:
                    if self.main_window:
                        self.main_window.statusBar().showMessage("No atom at click position")
            else:
                # Select this molecule
                if self.main_window:
                    self.main_window.select_molecule(self)

    def _find_atom_at_position(self, click_pos):
        """Find atom index at click position (within 20 pixel radius)"""
        click_radius = 20

        for atom_idx, atom_screen_pos in self.atom_positions.items():
            dx = click_pos.x() - atom_screen_pos.x()
            dy = click_pos.y() - atom_screen_pos.y()
            distance = math.sqrt(dx*dx + dy*dy)

            if distance <= click_radius:
                return atom_idx

        return None

    def set_selected(self, selected):
        """Set selection state"""
        self.selected = selected
        self.update_display()

    def update_cursor(self):
        """Update cursor based on edit mode"""
        if self.main_window and self.main_window.edit_highlights_mode:
            self.setCursor(QCursor(Qt.CrossCursor))
        else:
            self.setCursor(QCursor(Qt.PointingHandCursor))


class MoleculeAlignmentApp(QMainWindow):
    """Main application window"""
    def __init__(self):
        super().__init__()
        self.molecules = []  # List of MoleculeWidget
        self.selected_molecule = None
        self.reference_mol = None
        self.edit_highlights_mode = False  # Toggle for highlight editing mode

        self.init_ui()

    def init_ui(self):
        """Initialize user interface"""
        self.setWindowTitle("Interactive Molecular Structure Alignment")
        self.setGeometry(100, 100, 1400, 900)

        # Create toolbar
        self.create_toolbar()

        # Main widget
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)

        # Input section
        input_group = QGroupBox("Input Files")
        input_layout = QHBoxLayout()

        self.file_path = QLineEdit()
        self.file_path.setPlaceholderText("SDF or SMILES file path...")
        btn_file = QPushButton("Browse File")
        btn_file.clicked.connect(self.browse_file)

        self.ref_smiles = QLineEdit()
        self.ref_smiles.setPlaceholderText("Reference SMILES...")

        btn_load = QPushButton("Load & Align")
        btn_load.clicked.connect(self.load_molecules)
        btn_load.setStyleSheet("background-color: #4CAF50; color: white; font-weight: bold;")

        # Columns selector
        self.columns_spin = QSpinBox()
        self.columns_spin.setMinimum(1)
        self.columns_spin.setMaximum(10)
        self.columns_spin.setValue(4)
        self.columns_spin.setToolTip("Number of columns in grid")

        input_layout.addWidget(QLabel("File:"))
        input_layout.addWidget(self.file_path)
        input_layout.addWidget(btn_file)
        input_layout.addWidget(QLabel("Reference:"))
        input_layout.addWidget(self.ref_smiles)
        input_layout.addWidget(QLabel("Columns:"))
        input_layout.addWidget(self.columns_spin)
        input_layout.addWidget(btn_load)
        input_group.setLayout(input_layout)
        main_layout.addWidget(input_group)

        # Label settings section
        label_group = QGroupBox("Label Settings")
        label_layout = QHBoxLayout()

        # Font size selector
        self.font_size_combo = QComboBox()
        self.font_size_combo.addItems(["16", "20", "24", "28", "32", "36", "40", "48"])
        self.font_size_combo.setCurrentText("24")
        self.font_size_combo.currentIndexChanged.connect(self.update_all_displays)

        # Property selector for second line
        self.property_combo = QComboBox()
        self.property_combo.addItem("Catalog ID only", "catalog_only")
        self.property_combo.addItem("Price, USD", "Price, USD")
        self.property_combo.setCurrentIndex(1)  # Default to showing price
        self.property_combo.currentIndexChanged.connect(self.update_all_displays)

        label_layout.addWidget(QLabel("Font Size:"))
        label_layout.addWidget(self.font_size_combo)
        label_layout.addWidget(QLabel("Second Line:"))
        label_layout.addWidget(self.property_combo)
        label_layout.addStretch()
        label_group.setLayout(label_layout)
        main_layout.addWidget(label_group)

        # Info label
        self.info_label = QLabel("💡 Tip: Click molecule to select | Enable 'Edit Highlights Mode' then click atoms")
        self.info_label.setStyleSheet("color: #666; font-style: italic; padding: 5px;")
        main_layout.addWidget(self.info_label)

        # Molecule grid with scroll area
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumHeight(600)

        self.grid_widget = QWidget()
        self.grid_layout = QGridLayout(self.grid_widget)
        self.grid_layout.setSpacing(10)
        scroll.setWidget(self.grid_widget)

        main_layout.addWidget(scroll)

        # Status bar
        self.statusBar().showMessage("Load an SDF or SMILES file and reference SMILES to begin")

    def create_toolbar(self):
        """Create toolbar with actions"""
        toolbar = QToolBar()
        self.addToolBar(toolbar)

        # Mode toggle
        self.act_edit_mode = QAction("✏️ Edit Highlights Mode", self)
        self.act_edit_mode.setCheckable(True)
        self.act_edit_mode.setToolTip("Toggle: Click atoms to add/remove highlights")
        self.act_edit_mode.triggered.connect(self.toggle_edit_mode)
        toolbar.addAction(self.act_edit_mode)

        toolbar.addSeparator()

        # Highlight controls
        act_atom_selector = QAction("📝 Atom List...", self)
        act_atom_selector.triggered.connect(self.open_atom_selector)
        toolbar.addAction(act_atom_selector)

        act_reset_hl = QAction("Reset Highlights", self)
        act_reset_hl.triggered.connect(self.reset_highlights)
        toolbar.addAction(act_reset_hl)

        act_clear_hl = QAction("Clear All Highlights", self)
        act_clear_hl.triggered.connect(self.clear_highlights)
        toolbar.addAction(act_clear_hl)

        toolbar.addSeparator()

        # Rotation controls
        act_rot_cw = QAction("Rotate ↻ 90°", self)
        act_rot_cw.triggered.connect(lambda: self.rotate_selected(90))
        toolbar.addAction(act_rot_cw)

        act_rot_ccw = QAction("Rotate ↺ 90°", self)
        act_rot_ccw.triggered.connect(lambda: self.rotate_selected(-90))
        toolbar.addAction(act_rot_ccw)

        toolbar.addSeparator()

        # Mirror controls
        act_mirror_h = QAction("Mirror Horizontal ↔", self)
        act_mirror_h.triggered.connect(lambda: self.mirror_selected('h'))
        toolbar.addAction(act_mirror_h)

        act_mirror_v = QAction("Mirror Vertical ↕", self)
        act_mirror_v.triggered.connect(lambda: self.mirror_selected('v'))
        toolbar.addAction(act_mirror_v)

        toolbar.addSeparator()

        # Export
        act_export = QAction("💾 Export High-Res PNG", self)
        act_export.triggered.connect(self.export_image)
        toolbar.addAction(act_export)

        # Reset
        act_reset = QAction("Reset Selected", self)
        act_reset.triggered.connect(self.reset_selected)
        toolbar.addAction(act_reset)

    def browse_file(self):
        """Browse for SDF or SMILES file"""
        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Select Molecule File",
            "",
            "Molecule Files (*.sdf *.smi *.smiles);;SDF Files (*.sdf);;SMILES Files (*.smi *.smiles);;All Files (*)"
        )
        if filename:
            self.file_path.setText(filename)

    def load_molecules(self):
        """Load and align molecules from SDF or SMILES file"""
        file_path = self.file_path.text()
        ref_smiles = self.ref_smiles.text()

        if not file_path or not ref_smiles:
            QMessageBox.warning(self, "Input Required", "Please provide both molecule file and reference SMILES")
            return

        try:
            # Clear existing molecules
            for i in reversed(range(self.grid_layout.count())):
                self.grid_layout.itemAt(i).widget().setParent(None)
            self.molecules.clear()

            # Load molecules based on file extension
            mols = []
            if file_path.lower().endswith('.sdf'):
                mols = self._load_sdf(file_path)
            elif file_path.lower().endswith(('.smi', '.smiles')):
                mols = self._load_smiles(file_path)
            else:
                QMessageBox.warning(self, "Error", "Unsupported file format. Use .sdf, .smi, or .smiles")
                return

            if not mols:
                QMessageBox.warning(self, "Error", "No valid molecules found in file")
                return

            # Parse reference
            self.reference_mol = Chem.MolFromSmiles(ref_smiles)
            if self.reference_mol is None:
                QMessageBox.warning(self, "Error", "Invalid reference SMILES")
                return

            # Align and display molecules
            self._align_and_display(mols)

            self.statusBar().showMessage(f"Loaded {len(mols)} molecules")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load molecules: {str(e)}")

    def _load_sdf(self, file_path):
        """Load molecules from SDF file"""
        suppl = Chem.SDMolSupplier(file_path)
        return [mol for mol in suppl if mol is not None]

    def _load_smiles(self, file_path):
        """Load molecules from SMILES file with header support"""
        mols = []
        headers = None

        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue

                # Split by tab or space
                if '\t' in line:
                    parts = line.split('\t')
                else:
                    parts = line.split()

                if not parts:
                    continue

                # First non-comment line is the header
                if headers is None:
                    headers = parts
                    # Validate: first column should look like a header (not a SMILES)
                    # If it looks like SMILES, treat it as data instead
                    first_col_upper = parts[0].upper()
                    if first_col_upper in ['SMILES', 'SMI', 'SMILE'] or not self._is_likely_smiles(parts[0]):
                        # This is a header row
                        continue
                    else:
                        # No header, use default column names
                        headers = ['SMILES', 'Catalog ID'] + [f'Property_{i}' for i in range(2, len(parts))]
                        # Process this line as data (fall through)

                # Parse SMILES (first column)
                smiles = parts[0]
                mol = Chem.MolFromSmiles(smiles)

                if mol is not None:
                    # Set Catalog ID (second column, or default)
                    if len(parts) > 1 and parts[1].strip():
                        mol.SetProp("Catalog ID", parts[1].strip())
                    else:
                        mol.SetProp("Catalog ID", f"Mol_{line_num}")

                    # Set additional properties using header names
                    for col_idx in range(2, len(parts)):
                        if col_idx < len(headers):
                            prop_name = headers[col_idx].strip()
                            prop_value = parts[col_idx].strip()
                            if prop_value:  # Only set non-empty values
                                mol.SetProp(prop_name, prop_value)

                    mols.append(mol)

        return mols

    def _is_likely_smiles(self, text):
        """Check if text looks like a SMILES string"""
        # Simple heuristic: SMILES typically contain characters like c, C, (, ), =, #
        smiles_chars = set('cCnNoOpPsSFfIiBbr[]()=#@+-/\\123456789')
        text_chars = set(text)
        # If more than 50% of characters are SMILES-like, probably a SMILES
        common = len(text_chars & smiles_chars)
        return common / len(text_chars) > 0.5 if text_chars else False

    def _align_and_display(self, mols):
        """Align molecules to reference and display in grid"""
        # Find MCS between reference and each molecule
        match_list = []
        mcs_sizes = []
        rdDepictor.SetPreferCoordGen(True)

        for mol in mols:
            mcs = rdFMCS.FindMCS([self.reference_mol, mol],
                                maximizeBonds=True,
                                ringMatchesRingOnly=True,
                                completeRingsOnly=True,
                                timeout=10)  # 10 second timeout

            if mcs and mcs.numAtoms > 0:
                # Get the MCS as a molecule
                mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
                if mcs_mol:
                    # Find which atoms in the molecule match the MCS
                    match_atoms = mol.GetSubstructMatch(mcs_mol)
                    match_list.append(match_atoms)
                    mcs_sizes.append(mcs.numAtoms)
                else:
                    match_list.append(())
                    mcs_sizes.append(0)
            else:
                # No MCS found, use empty tuple
                match_list.append(())
                mcs_sizes.append(0)

        # Check if we found any matches
        if not any(match_list):
            QMessageBox.warning(self, "Warning",
                              "No common substructure found between reference and molecules.\n"
                              "Molecules will be displayed without alignment.")
        else:
            avg_mcs = sum(mcs_sizes) / len(mcs_sizes) if mcs_sizes else 0
            self.statusBar().showMessage(
                f"Found MCS matches: avg {avg_mcs:.1f} atoms, "
                f"min {min(mcs_sizes) if mcs_sizes else 0}, "
                f"max {max(mcs_sizes) if mcs_sizes else 0}"
            )

        # Align structures based on first molecule's match
        if match_list[0]:
            AllChem.Compute2DCoords(mols[0], useRingTemplates=True)
            coords = [mols[0].GetConformer().GetAtomPosition(x) for x in match_list[0]]
            coords2D = [Geometry.Point2D(pt.x, pt.y) for pt in coords]

            for mol_idx, mol in enumerate(mols[1:], start=1):
                if match_list[mol_idx] and len(match_list[mol_idx]) == len(coords2D):
                    coord_dict = {}
                    for i, coord in enumerate(coords2D):
                        try:
                            coord_dict[match_list[mol_idx][i]] = coord
                        except IndexError:
                            continue
                    AllChem.Compute2DCoords(mol, coordMap=coord_dict)
                else:
                    # No match or size mismatch, compute coords normally
                    AllChem.Compute2DCoords(mol, useRingTemplates=True)
        else:
            # No matches at all, just compute coords for all
            for mol in mols:
                AllChem.Compute2DCoords(mol)

        # Create molecule widgets
        columns = self.columns_spin.value()

        # Collect all unique properties from molecules
        all_properties = set()
        for mol in mols:
            all_properties.update(mol.GetPropsAsDict().keys())

        # Update property dropdown with available properties
        current_selection = self.property_combo.currentData()
        self.property_combo.blockSignals(True)  # Prevent triggering update during population
        self.property_combo.clear()
        self.property_combo.addItem("Catalog ID only", "catalog_only")

        # Add all available properties
        for prop in sorted(all_properties):
            if prop not in ["Catalog ID", "_Name"]:  # Skip catalog ID (shown as line 1) and internal props
                self.property_combo.addItem(prop, prop)

        # Restore previous selection if it exists
        idx = self.property_combo.findData(current_selection)
        if idx >= 0:
            self.property_combo.setCurrentIndex(idx)
        else:
            # Default to Price, USD if available
            idx = self.property_combo.findData("Price, USD")
            if idx >= 0:
                self.property_combo.setCurrentIndex(idx)

        self.property_combo.blockSignals(False)

        for idx, mol in enumerate(mols):
            catalog_id = mol.GetProp("Catalog ID") if mol.HasProp("Catalog ID") else f"Mol {idx+1}"

            all_atoms = list(range(mol.GetNumAtoms()))
            mol_data = MoleculeData(mol, match_list[idx], all_atoms, catalog_id)

            mol_widget = MoleculeWidget(mol_data, img_size=(350, 350), main_window=self)
            self.molecules.append(mol_widget)

            row = idx // columns
            col = idx % columns
            self.grid_layout.addWidget(mol_widget, row, col)

    def select_molecule(self, molecule_widget):
        """Select a molecule for editing"""
        if self.selected_molecule:
            self.selected_molecule.set_selected(False)

        self.selected_molecule = molecule_widget
        molecule_widget.set_selected(True)

        self.statusBar().showMessage(f"Selected: {molecule_widget.mol_data.catalog_id}")

    def update_all_displays(self):
        """Update all molecule displays (when font size or property changes)"""
        if not self.molecules:
            return

        # Update all molecules
        for mol_widget in self.molecules:
            mol_widget.update_display()

        # Update status bar to confirm
        font_size = self.font_size_combo.currentText()
        self.statusBar().showMessage(f"Updated labels: Font size {font_size}pt")

    def toggle_edit_mode(self, checked):
        """Toggle between selection mode and highlight editing mode"""
        self.edit_highlights_mode = checked

        # Update all molecule cursors
        for mol_widget in self.molecules:
            mol_widget.update_cursor()

        # Update info label
        if checked:
            self.info_label.setText("✏️ EDIT MODE: Click atoms to add/remove highlights | Click button again to exit")
            self.info_label.setStyleSheet("color: #2196F3; font-weight: bold; padding: 5px; background-color: #E3F2FD;")
            self.statusBar().showMessage("Edit Highlights Mode: Click on atoms to toggle highlights")
        else:
            self.info_label.setText("💡 Tip: Click molecule to select | Enable 'Edit Highlights Mode' then click atoms")
            self.info_label.setStyleSheet("color: #666; font-style: italic; padding: 5px;")
            self.statusBar().showMessage("Selection Mode: Click molecules to select them")

    def open_atom_selector(self):
        """Open dialog to select atoms for highlighting"""
        if not self.selected_molecule:
            QMessageBox.information(self, "Info", "Please select a molecule first")
            return

        dialog = AtomSelectorDialog(self.selected_molecule.mol_data, self)
        if dialog.exec_() == QDialog.Accepted:
            # Update highlights based on selection
            selected_atoms = dialog.get_selected_atoms()
            self.selected_molecule.mol_data.highlighted_atoms = selected_atoms
            self.selected_molecule.update_display()

    def reset_highlights(self):
        """Reset highlights to original for selected molecule"""
        if not self.selected_molecule:
            QMessageBox.information(self, "Info", "Please select a molecule first")
            return

        self.selected_molecule.mol_data.reset_highlights()
        self.selected_molecule.update_display()

    def clear_highlights(self):
        """Clear all highlights for selected molecule"""
        if not self.selected_molecule:
            QMessageBox.information(self, "Info", "Please select a molecule first")
            return

        self.selected_molecule.mol_data.highlighted_atoms.clear()
        self.selected_molecule.update_display()

    def rotate_selected(self, degrees):
        """Rotate selected molecule"""
        if not self.selected_molecule:
            QMessageBox.information(self, "Info", "Please select a molecule first")
            return

        mol_data = self.selected_molecule.mol_data
        mol_data.rotation = (mol_data.rotation + degrees) % 360
        self.selected_molecule.update_display()

    def mirror_selected(self, direction):
        """Mirror selected molecule"""
        if not self.selected_molecule:
            QMessageBox.information(self, "Info", "Please select a molecule first")
            return

        mol_data = self.selected_molecule.mol_data
        if direction == 'h':
            mol_data.mirror_h = not mol_data.mirror_h
        else:
            mol_data.mirror_v = not mol_data.mirror_v

        self.selected_molecule.update_display()

    def reset_selected(self):
        """Reset all transformations for selected molecule"""
        if not self.selected_molecule:
            QMessageBox.information(self, "Info", "Please select a molecule first")
            return

        mol_data = self.selected_molecule.mol_data
        mol_data.rotation = 0
        mol_data.mirror_h = False
        mol_data.mirror_v = False
        mol_data.reset_highlights()
        self.selected_molecule.update_display()

    def export_image(self):
        """Export high-resolution grid image"""
        if not self.molecules:
            QMessageBox.warning(self, "Warning", "No molecules to export")
            return

        filename, _ = QFileDialog.getSaveFileName(self, "Save Image", "molecular_grid.png", "PNG Files (*.png)")
        if not filename:
            return

        try:
            # Render at high resolution (2x)
            img_size = (700, 700)  # High-res size
            columns = self.columns_spin.value()
            rows = (len(self.molecules) + columns - 1) // columns
            padding = 20

            # Render each molecule at high resolution
            high_res_images = []
            for mol_widget in self.molecules:
                mol_data = mol_widget.mol_data

                # Get transformed molecule
                mol = mol_data.get_transformed_mol()

                # Render molecule
                d2d = rdMolDraw2D.MolDraw2DCairo(img_size[0], img_size[1])

                # Set drawing options for thicker bonds in export
                draw_options = d2d.drawOptions()
                draw_options.bondLineWidth = 4.0  # Thicker for high-res export

                highlight_atoms = list(mol_data.highlighted_atoms)
                highlight_colors = {idx: (0.596, 0.984, 0.596) for idx in highlight_atoms}
                highlight_radii = {idx: 0.8 for idx in highlight_atoms}

                d2d.DrawMolecule(mol,
                                highlightAtoms=highlight_atoms,
                                highlightAtomColors=highlight_colors,
                                highlightAtomRadii=highlight_radii)
                d2d.FinishDrawing()

                mol_img_bytes = d2d.GetDrawingText()
                bio = BytesIO(mol_img_bytes)
                pil_img = Image.open(bio).convert('RGB')

                # Add text
                from PIL import ImageDraw, ImageFont

                font_size = int(self.font_size_combo.currentText())
                font_size_export = int(font_size * 2.0)  # Scale up 2x for high-res (700 vs 350)

                # Try to load font with multiple fallbacks
                font = None
                font_names = ["arial.ttf", "Arial.ttf", "DejaVuSans.ttf", "Helvetica.ttf"]
                for font_name in font_names:
                    try:
                        font = ImageFont.truetype(font_name, font_size_export)
                        break
                    except:
                        continue

                if font is None:
                    font = ImageFont.load_default()

                line1 = mol_data.catalog_id

                # Get second line based on property selection
                line2 = None
                selected_prop = self.property_combo.currentData()
                if selected_prop != "catalog_only":
                    if mol_data.mol.HasProp(selected_prop):
                        prop_value = mol_data.mol.GetProp(selected_prop)
                        if selected_prop == "Price, USD":
                            line2 = f'Price: ${prop_value}'
                        else:
                            line2 = f'{selected_prop}: {prop_value}'

                draw = ImageDraw.Draw(pil_img)
                bbox1 = draw.textbbox((0, 0), line1, font=font)
                w1, h1 = bbox1[2] - bbox1[0], bbox1[3] - bbox1[1]

                if line2:
                    bbox2 = draw.textbbox((0, 0), line2, font=font)
                    w2, h2 = bbox2[2] - bbox2[0], bbox2[3] - bbox2[1]
                    text_height = h1 + h2 + 40
                else:
                    text_height = h1 + 20

                new_img = Image.new('RGB', (pil_img.width, pil_img.height + text_height), 'white')
                new_img.paste(pil_img, (0, 0))

                draw = ImageDraw.Draw(new_img)
                x1 = (pil_img.width - w1) // 2
                y1 = pil_img.height + 10
                draw.text((x1, y1), line1, font=font, fill='black')

                if line2:
                    x2 = (pil_img.width - w2) // 2
                    y2 = y1 + h1 + 10
                    draw.text((x2, y2), line2, font=font, fill='black')

                high_res_images.append(new_img)

            # Create grid
            img_width = img_size[0]
            img_height = high_res_images[0].height

            grid_width = columns * (img_width + padding) + padding
            grid_height = rows * (img_height + padding) + padding
            grid_img = Image.new('RGB', (grid_width, grid_height), 'white')

            for idx, img in enumerate(high_res_images):
                x = (idx % columns) * (img_width + padding) + padding
                y = (idx // columns) * (img_height + padding) + padding
                grid_img.paste(img, (x, y))

            grid_img.save(filename, dpi=(600, 600))

            QMessageBox.information(self, "Success", f"High-resolution image saved to:\n{filename}\n\nResolution: 600 DPI")
            self.statusBar().showMessage(f"Exported to {filename} (600 DPI)")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to export image: {str(e)}")


def main():
    app = QApplication(sys.argv)
    window = MoleculeAlignmentApp()
    window.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()

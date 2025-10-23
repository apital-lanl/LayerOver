import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import os
import copy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

class DIW_PSPP:
    def __init__(self, root):
        self.root = root
        self.root.title("DIW PSPP Application")
        self.root.geometry("1200x800")
        
        # Initialize dictionaries and variables
        self.user_defined_structure = "S3HS"
        self.parsed_user_structure = {"ALL": {
            "layer_type": "standard",
            "nozzle_size": 0.25,
            "material": "LL20",
            "angular_offset": 0,
            "lateral_offset": 0,
            "stochastic_modifiers": {},
            "type_modifiers": {}
        }}
        self.current_target_structure = copy.deepcopy(self.parsed_user_structure)
        self.session_considered_structures = {}
        self.structure_data_elements = {}
        self.selected_structure_data = {}
        
        # Material viscosity mapping
        self.material_viscosity = {
            "LL20": 10000,
            "LL50": 25000,
            "LL60": 30000,
            "LL70": 35000
        }
        
        # Create main panels
        self.create_panels()
        
        # Create menu
        self.create_menu()
        
        # Parse initial structure
        self.parse_user_defined_structure()

    def create_panels(self):
        # Create PanedWindow for three main panels
        main_paned = ttk.PanedWindow(self.root, orient=tk.HORIZONTAL)
        main_paned.pack(fill=tk.BOTH, expand=True)
        
        # Left panel - Controls (with scrollbar)
        left_container = ttk.Frame(main_paned)
        main_paned.add(left_container, weight=1)
        
        # Create canvas and scrollbar for left panel
        left_canvas = tk.Canvas(left_container)
        left_scrollbar = ttk.Scrollbar(left_container, orient="vertical", command=left_canvas.yview)
        left_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        left_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Create frame inside canvas for content
        left_frame = ttk.Frame(left_canvas)
        left_canvas.create_window((0, 0), window=left_frame, anchor="nw")
        
        # Configure canvas scrolling
        left_frame.bind("<Configure>", lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))
        left_canvas.configure(yscrollcommand=left_scrollbar.set)
        
        # Middle panel - Structure Data (with scrollbar)
        middle_container = ttk.LabelFrame(main_paned, text="Structure Data")
        main_paned.add(middle_container, weight=1)
        
        # Create canvas and scrollbar for middle panel
        middle_canvas = tk.Canvas(middle_container)
        middle_scrollbar = ttk.Scrollbar(middle_container, orient="vertical", command=middle_canvas.yview)
        middle_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        middle_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Create frame inside canvas for content
        middle_frame = ttk.Frame(middle_canvas)
        middle_canvas.create_window((0, 0), window=middle_frame, anchor="nw")
        
        # Configure canvas scrolling
        middle_frame.bind("<Configure>", lambda e: middle_canvas.configure(scrollregion=middle_canvas.bbox("all")))
        middle_canvas.configure(yscrollcommand=middle_scrollbar.set)
        
        # Right panel - Data View (with scrollbar)
        right_container = ttk.LabelFrame(main_paned, text="Data View")
        main_paned.add(right_container, weight=2)
        
        # Create canvas and scrollbar for right panel
        right_canvas = tk.Canvas(right_container)
        right_scrollbar = ttk.Scrollbar(right_container, orient="vertical", command=right_canvas.yview)
        right_scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        right_canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        
        # Create frame inside canvas for content
        right_frame = ttk.Frame(right_canvas)
        right_canvas.create_window((0, 0), window=right_frame, anchor="nw")
        
        # Configure canvas scrolling
        right_frame.bind("<Configure>", lambda e: right_canvas.configure(scrollregion=right_canvas.bbox("all")))
        right_canvas.configure(yscrollcommand=right_scrollbar.set)
        
        # Setup left panel sections
        self.setup_left_panel(left_frame)
        
        # Setup middle panel
        self.setup_middle_panel(middle_frame)
        
        # Setup right panel
        self.setup_right_panel(right_frame)
        
        # Store canvas references for later use
        self.left_canvas = left_canvas
        self.middle_canvas = middle_canvas
        self.right_canvas = right_canvas

    def setup_left_panel(self, parent):
        # Create sections using LabelFrames
        process_frame = ttk.LabelFrame(parent, text="Process")
        process_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        structure_frame = ttk.LabelFrame(parent, text="Structure")
        structure_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        properties_frame = ttk.LabelFrame(parent, text="Properties")
        properties_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        data_controls_frame = ttk.LabelFrame(parent, text="Data Controls")
        data_controls_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Setup Process section
        self.setup_process_section(process_frame)
        
        # Setup Structure section
        self.setup_structure_section(structure_frame)
        
        # Setup Properties section
        self.setup_properties_section(properties_frame)
        
        # Setup Data Controls section
        self.setup_data_controls_section(data_controls_frame)

    def setup_process_section(self, parent):
        # Material selection
        material_frame = ttk.Frame(parent)
        material_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(material_frame, text="Material:").pack(side=tk.LEFT)
        
        self.material_var = tk.StringVar(value="LL20")
        material_combo = ttk.Combobox(material_frame, textvariable=self.material_var, 
                                      values=["LL20", "LL50", "LL60", "LL70"], state="readonly")
        material_combo.pack(side=tk.LEFT, padx=5)
        material_combo.bind("<<ComboboxSelected>>", self.update_viscosity)
        
        # Machine selection
        machine_frame = ttk.Frame(parent)
        machine_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(machine_frame, text="Machine:").pack(side=tk.LEFT)
        
        self.machine_var = tk.StringVar(value="Aerotech 5-axis")
        machine_combo = ttk.Combobox(machine_frame, textvariable=self.machine_var, 
                                     values=["Aerotech 5-axis", "Aerotech 3-axis", "Aerotech \"5-1\"-axis"], state="readonly")
        machine_combo.pack(side=tk.LEFT, padx=5)
        
        # Viscosity subsection
        viscosity_frame = ttk.LabelFrame(parent, text="Viscosity")
        viscosity_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Viscosity value
        visc_value_frame = ttk.Frame(viscosity_frame)
        visc_value_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(visc_value_frame, text="Value:").pack(side=tk.LEFT)
        
        self.viscosity_var = tk.IntVar(value=self.material_viscosity["LL20"])
        viscosity_spinbox = ttk.Spinbox(visc_value_frame, from_=1000, to=100000, increment=100, textvariable=self.viscosity_var)
        viscosity_spinbox.pack(side=tk.LEFT, padx=5)
        
        # Viscosity modifier
        visc_mod_frame = ttk.Frame(viscosity_frame)
        visc_mod_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(visc_mod_frame, text="Modifier:").pack(side=tk.LEFT)
        
        self.viscosity_mod_var = tk.IntVar(value=0)
        viscosity_mod_spinbox = ttk.Spinbox(visc_mod_frame, from_=-1000, to=1000, increment=10, textvariable=self.viscosity_mod_var)
        viscosity_mod_spinbox.pack(side=tk.LEFT, padx=5)
        
        # Add from data button and display
        visc_data_frame = ttk.Frame(viscosity_frame)
        visc_data_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(visc_data_frame, text="Add from data", command=self.add_viscosity_from_data).pack(side=tk.LEFT, padx=5)
        
        self.viscosity_file_var = tk.StringVar(value="No file selected")
        ttk.Label(visc_data_frame, textvariable=self.viscosity_file_var).pack(side=tk.LEFT, padx=5)

    def setup_structure_section(self, parent):
        # Structure Name
        name_frame = ttk.Frame(parent)
        name_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(name_frame, text="Structure Name:").pack(side=tk.LEFT)
        
        self.structure_name_var = tk.StringVar(value=self.user_defined_structure)
        structure_entry = ttk.Entry(name_frame, textvariable=self.structure_name_var)
        structure_entry.pack(side=tk.LEFT, padx=5, fill=tk.X, expand=True)
        structure_entry.bind("<Return>", lambda e: self.parse_user_defined_structure())
        
        ttk.Button(name_frame, text="Apply", command=self.parse_user_defined_structure).pack(side=tk.LEFT, padx=5)
        
        # Create notebook for layers
        self.layer_notebook = ttk.Notebook(parent)
        self.layer_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Initial layer tab will be created when parse_user_defined_structure is called

    def setup_properties_section(self, parent):
        # Thickness and Density
        basic_props_frame = ttk.Frame(parent)
        basic_props_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(basic_props_frame, text="Thickness (mm):").pack(side=tk.LEFT)
        self.thickness_var = tk.DoubleVar(value=1.0)
        thickness_spinbox = ttk.Spinbox(basic_props_frame, from_=0.1, to=10.0, increment=0.1, textvariable=self.thickness_var)
        thickness_spinbox.pack(side=tk.LEFT, padx=5)
        thickness_spinbox.bind("<Return>", lambda e: self.update_mechanical_response_labels())
        
        ttk.Label(basic_props_frame, text="Density (g/cc):").pack(side=tk.LEFT, padx=10)
        self.density_var = tk.DoubleVar(value=1.2)
        density_spinbox = ttk.Spinbox(basic_props_frame, from_=0.1, to=10.0, increment=0.1, textvariable=self.density_var)
        density_spinbox.pack(side=tk.LEFT, padx=5)
        
        # Mechanical Response section
        mech_frame = ttk.LabelFrame(parent, text="Mechanical Response (stress/mm)")
        mech_frame.pack(fill=tk.X, padx=5, pady=5)
        
        # Create 10 spinboxes for mechanical response
        self.mech_response_vars = []
        self.mech_response_mm_labels = []
        
        mech_spinbox_frame = ttk.Frame(mech_frame)
        mech_spinbox_frame.pack(fill=tk.X, padx=5, pady=5)
        
        percentages = [1, 5, 10, 25, 50, 75, 100, 125, 140, 150]
        
        for i, percent in enumerate(percentages):
            frame = ttk.Frame(mech_spinbox_frame)
            frame.pack(side=tk.LEFT, padx=2)
            
            var = tk.DoubleVar(value=0.0)
            self.mech_response_vars.append(var)
            
            spinbox = ttk.Spinbox(frame, from_=0, to=1000, increment=0.1, width=5, textvariable=var)
            spinbox.pack(side=tk.TOP)
            
            ttk.Label(frame, text=f"{percent}%").pack(side=tk.TOP)
            
            mm_label = ttk.Label(frame, text=f"{percent/100 * self.thickness_var.get():.2f}")
            mm_label.pack(side=tk.TOP)
            self.mech_response_mm_labels.append(mm_label)
        
        # Update mechanical response labels when thickness changes
        self.thickness_var.trace_add("write", lambda *args: self.update_mechanical_response_labels())

    def setup_data_controls_section(self, parent):
        # Experimental Data Filepath
        filepath_frame = ttk.Frame(parent)
        filepath_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Label(filepath_frame, text="Exp. Data Filepath:").pack(side=tk.LEFT)
        
        self.data_filepath_var = tk.StringVar(value="No file selected")
        ttk.Label(filepath_frame, textvariable=self.data_filepath_var).pack(side=tk.LEFT, padx=5)
        
        ttk.Button(filepath_frame, text="Change Data Filepath", command=self.change_data_filepath).pack(side=tk.LEFT, padx=5)
        
        # Checkboxes for display options
        display_frame = ttk.Frame(parent)
        display_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.show_closest_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(display_frame, text="Show closest matches", variable=self.show_closest_var).pack(anchor=tk.W)
        
        self.show_model_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(display_frame, text="Show model results", variable=self.show_model_var).pack(anchor=tk.W)
        
        # Model parameters subsection
        model_frame = ttk.LabelFrame(parent, text="Model parameters")
        model_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.opacity_model_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(model_frame, text="Opacity Model", variable=self.opacity_model_var).pack(anchor=tk.W)
        
        self.classification_model_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(model_frame, text="Classification Model", variable=self.classification_model_var).pack(anchor=tk.W)
        
        self.inverse_model_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(model_frame, text="Inverse Model", variable=self.inverse_model_var).pack(anchor=tk.W)
        
        self.custom_model_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(model_frame, text="Custom Model", variable=self.custom_model_var).pack(anchor=tk.W)
        
        # Data Options subsection
        data_options_frame = ttk.LabelFrame(parent, text="Data Options")
        data_options_frame.pack(fill=tk.X, padx=5, pady=5)
        
        self.plot_option_var = tk.StringVar(value="together")
        ttk.Radiobutton(data_options_frame, text="Plot Selected Together", variable=self.plot_option_var, value="together").pack(anchor=tk.W)
        ttk.Radiobutton(data_options_frame, text="Plot Selected Separately", variable=self.plot_option_var, value="separate").pack(anchor=tk.W)
        
        # Buttons for data actions
        buttons_frame = ttk.Frame(data_options_frame)
        buttons_frame.pack(fill=tk.X, padx=5, pady=5)
        
        ttk.Button(buttons_frame, text="View Selected", command=self.view_selected_data).pack(fill=tk.X, pady=2)
        ttk.Button(buttons_frame, text="Report Selected", command=self.report_selected_data).pack(fill=tk.X, pady=2)
        ttk.Button(buttons_frame, text="Open Selected in Analysis View", command=self.open_analysis_view).pack(fill=tk.X, pady=2)
        ttk.Button(buttons_frame, text="Refine Structure for Properties", command=self.refine_structure).pack(fill=tk.X, pady=2)
        ttk.Button(buttons_frame, text="Clear Data View", command=self.clear_data_view).pack(fill=tk.X, pady=2)

    def setup_middle_panel(self, parent):
        # Create Treeview for structure data
        self.structure_tree = ttk.Treeview(parent)
        self.structure_tree.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Configure treeview
        self.structure_tree["columns"] = ("type", "value")
        self.structure_tree.column("#0", width=150, minwidth=150)
        self.structure_tree.column("type", width=100, minwidth=100)
        self.structure_tree.column("value", width=100, minwidth=100)
        
        self.structure_tree.heading("#0", text="Item")
        self.structure_tree.heading("type", text="Type")
        self.structure_tree.heading("value", text="Value")
        
        # Create main sections
        self.exp_data_node = self.structure_tree.insert("", "end", text="Experimental Data", open=True)
        self.model_data_node = self.structure_tree.insert("", "end", text="Model Data", open=True)
        self.predictive_data_node = self.structure_tree.insert("", "end", text="Predictive Data", open=True)
        
        # Bind selection event
        self.structure_tree.bind("<<TreeviewSelect>>", self.on_tree_select)
        
        # Populate with sample data
        self.populate_experimental_tree()

    def setup_right_panel(self, parent):
        # Create notebook for data view tabs
        self.data_view_notebook = ttk.Notebook(parent)
        self.data_view_notebook.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)
        
        # Initial empty tab
        empty_frame = ttk.Frame(self.data_view_notebook)
        self.data_view_notebook.add(empty_frame, text="No Data")
        
        ttk.Label(empty_frame, text="Select data from the Structure Data panel to view").pack(expand=True)

    def create_menu(self):
        menubar = tk.Menu(self.root)
        
        # File menu
        file_menu = tk.Menu(menubar, tearoff=0)
        file_menu.add_command(label="New Project", command=self.new_project)
        file_menu.add_command(label="Open Project", command=self.open_project)
        file_menu.add_command(label="Save Project", command=self.save_project)
        file_menu.add_separator()
        file_menu.add_command(label="Export Data", command=self.export_data)
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=self.root.quit)
        menubar.add_cascade(label="File", menu=file_menu)
        
        # Refine menu
        refine_menu = tk.Menu(menubar, tearoff=0)
        refine_menu.add_command(label="Refine Structure", command=self.refine_structure)
        refine_menu.add_command(label="Optimize Parameters", command=self.optimize_parameters)
        menubar.add_cascade(label="Refine", menu=refine_menu)
        
        # Analyze menu
        analyze_menu = tk.Menu(menubar, tearoff=0)
        analyze_menu.add_command(label="Compare Structures", command=self.compare_structures)
        analyze_menu.add_command(label="Statistical Analysis", command=self.statistical_analysis)
        menubar.add_cascade(label="Analyze", menu=analyze_menu)
        
        # Report menu
        report_menu = tk.Menu(menubar, tearoff=0)
        report_menu.add_command(label="Generate Report", command=self.generate_report)
        report_menu.add_command(label="Export Figures", command=self.export_figures)
        menubar.add_cascade(label="Report", menu=report_menu)
        
        self.root.config(menu=menubar)

    def update_viscosity(self, event=None):
        material = self.material_var.get()
        self.viscosity_var.set(self.material_viscosity[material])

    def add_viscosity_from_data(self):
        filename = filedialog.askopenfilename(
            title="Select Viscosity Data File",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
        )
        if filename:
            self.viscosity_file_var.set(os.path.basename(filename))
            # Here you would add code to read the file and update viscosity values

    def parse_user_defined_structure(self):
        structure_name = self.structure_name_var.get()
        self.user_defined_structure = structure_name
        
        # This is a simplified parsing logic - in a real application, 
        # you would have more sophisticated parsing based on the structure name
        
        # For demonstration, let's say the first character indicates the number of layers
        try:
            num_layers = int(structure_name[0]) if structure_name[0].isdigit() else 1
        except:
            num_layers = 1
        
        # Reset the parsed structure
        self.parsed_user_structure = {}
        
        # Create layers
        for i in range(1, num_layers + 1):
            self.parsed_user_structure[f"Layer_{i}"] = {
                "layer_type": "standard",
                "nozzle_size": 0.25,
                "material": self.material_var.get(),
                "pitch": 0.5,
                "angular_offset": (i-1) * 45,
                "lateral_offset": 0,
                "stochastic_modifiers": {},
                "type_modifiers": {}
            }
        
        # Update current target structure
        self.current_target_structure = copy.deepcopy(self.parsed_user_structure)
        
        # Check if this is a new structure to add to session_considered_structures
        structure_hash = str(self.current_target_structure)  # Simple hash for demo
        if structure_hash not in [str(s) for s in self.session_considered_structures.values()]:
            self.session_considered_structures[structure_name] = copy.deepcopy(self.current_target_structure)
        
        # Update the layer tabs
        self.update_layer_tabs()
        
        return self.parsed_user_structure

    def update_layer_tabs(self):
        # Clear existing tabs
        for tab in self.layer_notebook.tabs():
            self.layer_notebook.forget(tab)
        
        # Create new tabs for each layer
        for layer_name, layer_data in self.current_target_structure.items():
            layer_frame = ttk.Frame(self.layer_notebook)
            self.layer_notebook.add(layer_frame, text=layer_name)
            
            # Create input fields for layer properties
            row = 0
            
            # Nozzle Size
            ttk.Label(layer_frame, text="Nozzle Size:").grid(row=row, column=0, sticky=tk.W, padx=5, pady=5)
            nozzle_var = tk.DoubleVar(value=layer_data["nozzle_size"])
            nozzle_spinbox = ttk.Spinbox(layer_frame, from_=0.1, to=1.0, increment=0.05, textvariable=nozzle_var)
            nozzle_spinbox.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
            row += 1
            
            # Material
            ttk.Label(layer_frame, text="Material:").grid(row=row, column=0, sticky=tk.W, padx=5, pady=5)
            material_var = tk.StringVar(value=layer_data["material"])
            material_combo = ttk.Combobox(layer_frame, textvariable=material_var, 
                                          values=["LL20", "LL50", "LL60", "LL70"], state="readonly")
            material_combo.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
            row += 1
            
            # Pitch
            ttk.Label(layer_frame, text="Pitch:").grid(row=row, column=0, sticky=tk.W, padx=5, pady=5)
            pitch_var = tk.DoubleVar(value=layer_data.get("pitch", 0.5))
            pitch_spinbox = ttk.Spinbox(layer_frame, from_=0.1, to=2.0, increment=0.1, textvariable=pitch_var)
            pitch_spinbox.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
            row += 1
            
            # Angular Offset
            ttk.Label(layer_frame, text="Angular Offset:").grid(row=row, column=0, sticky=tk.W, padx=5, pady=5)
            angle_var = tk.DoubleVar(value=layer_data["angular_offset"])
            angle_spinbox = ttk.Spinbox(layer_frame, from_=0, to=180, increment=5, textvariable=angle_var)
            angle_spinbox.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
            row += 1
            
            # Lateral Offset
            ttk.Label(layer_frame, text="Lateral Offset:").grid(row=row, column=0, sticky=tk.W, padx=5, pady=5)
            offset_var = tk.DoubleVar(value=layer_data["lateral_offset"])
            offset_spinbox = ttk.Spinbox(layer_frame, from_=0, to=1.0, increment=0.05, textvariable=offset_var)
            offset_spinbox.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
            row += 1
            
            # Layer Type
            ttk.Label(layer_frame, text="Layer Type:").grid(row=row, column=0, sticky=tk.W, padx=5, pady=5)
            layer_type_var = tk.StringVar(value=layer_data["layer_type"])
            layer_type_entry = ttk.Entry(layer_frame, textvariable=layer_type_var)
            layer_type_entry.grid(row=row, column=1, sticky=tk.W, padx=5, pady=5)
            row += 1
            
            # Save button
            ttk.Button(layer_frame, text="Apply Changes", 
                      command=lambda ln=layer_name, nv=nozzle_var, mv=material_var, pv=pitch_var, 
                                    av=angle_var, ov=offset_var, ltv=layer_type_var: 
                      self.update_layer_data(ln, nv, mv, pv, av, ov, ltv)).grid(row=row, column=0, columnspan=2, pady=10)

    def update_layer_data(self, layer_name, nozzle_var, material_var, pitch_var, angle_var, offset_var, layer_type_var):
        # Update the current target structure with values from the UI
        self.current_target_structure[layer_name]["nozzle_size"] = nozzle_var.get()
        self.current_target_structure[layer_name]["material"] = material_var.get()
        self.current_target_structure[layer_name]["pitch"] = pitch_var.get()
        self.current_target_structure[layer_name]["angular_offset"] = angle_var.get()
        self.current_target_structure[layer_name]["lateral_offset"] = offset_var.get()
        self.current_target_structure[layer_name]["layer_type"] = layer_type_var.get()
        
        # Check if this is a new structure to add to session_considered_structures
        structure_hash = str(self.current_target_structure)  # Simple hash for demo
        if structure_hash not in [str(s) for s in self.session_considered_structures.values()]:
            self.session_considered_structures[self.user_defined_structure] = copy.deepcopy(self.current_target_structure)
        
        messagebox.showinfo("Layer Updated", f"Layer {layer_name} has been updated")

    def update_mechanical_response_labels(self):
        thickness = self.thickness_var.get()
        percentages = [1, 5, 10, 25, 50, 75, 100, 125, 140, 150]
        
        for i, percent in enumerate(percentages):
            self.mech_response_mm_labels[i].config(text=f"{percent/100 * thickness:.2f}")

    def change_data_filepath(self):
        filename = filedialog.askopenfilename(
            title="Select Experimental Data File",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
        )
        if filename:
            self.data_filepath_var.set(os.path.basename(filename))
            # Here you would add code to read the file and update the data

    def populate_experimental_tree(self):
        # This function would look at the structure information and populate the tree
        # For demonstration, we'll just add some example items
        
        # Clear existing items
        for item in self.structure_tree.get_children(self.exp_data_node):
            self.structure_tree.delete(item)
        
        # Add example items
        for i in range(3):
            item_id = self.structure_tree.insert(self.exp_data_node, "end", text=f"Sample {i+1}", 
                                               values=("Experimental", f"Value {i+1}"))
            
            # Add some sub-items
            self.structure_tree.insert(item_id, "end", text=f"Property {i+1}.1", 
                                     values=("Stress-Strain", f"{i*10 + 5}"))
            self.structure_tree.insert(item_id, "end", text=f"Property {i+1}.2", 
                                     values=("Modulus", f"{i*100 + 50}"))
        
        # Update the structure_data_elements dictionary
        self.structure_data_elements["experimental"] = {
            "Sample 1": {"type": "Experimental", "value": "Value 1", "properties": {"Stress-Strain": 5, "Modulus": 50}},
            "Sample 2": {"type": "Experimental", "value": "Value 2", "properties": {"Stress-Strain": 15, "Modulus": 150}},
            "Sample 3": {"type": "Experimental", "value": "Value 3", "properties": {"Stress-Strain": 25, "Modulus": 250}}
        }

    def on_tree_select(self, event):
        # Get selected items
        selected_items = self.structure_tree.selection()
        
        # Clear previous selection
        self.selected_structure_data = {}
        
        # Add selected items to the selected_structure_data dictionary
        for item_id in selected_items:
            item_text = self.structure_tree.item(item_id, "text")
            item_values = self.structure_tree.item(item_id, "values")
            
            # Skip if it's a main category
            if item_id in [self.exp_data_node, self.model_data_node, self.predictive_data_node]:
                continue
            
            # Get parent if it's a sub-item
            parent_id = self.structure_tree.parent(item_id)
            if parent_id in [self.exp_data_node, self.model_data_node, self.predictive_data_node]:
                # It's a top-level item under a main category
                self.selected_structure_data[item_text] = {
                    "type": item_values[0] if item_values else "",
                    "value": item_values[1] if len(item_values) > 1 else "",
                    "properties": {}
                }
            else:
                # It's a sub-item
                parent_text = self.structure_tree.item(parent_id, "text")
                if parent_text not in self.selected_structure_data:
                    parent_values = self.structure_tree.item(parent_id, "values")
                    self.selected_structure_data[parent_text] = {
                        "type": parent_values[0] if parent_values else "",
                        "value": parent_values[1] if len(parent_values) > 1 else "",
                        "properties": {}
                    }
                
                # Add the property
                self.selected_structure_data[parent_text]["properties"][item_text] = {
                    "type": item_values[0] if item_values else "",
                    "value": item_values[1] if len(item_values) > 1 else ""
                }

    def view_selected_data(self):
        # Parse the selected data
        tabs_info = self.parse_selected_data()
        
        # Create the data views
        self.create_selected_data_views(tabs_info)

    def parse_selected_data(self):
        # Determine how many tabs to create based on plot option
        if not self.selected_structure_data:
            return [{"name": "No Data", "items": []}]
        
        if self.plot_option_var.get() == "together":
            # Create a single tab with all selected items
            return [{"name": "Combined Data", "items": list(self.selected_structure_data.keys())}]
        else:
            # Create a separate tab for each selected item
            return [{"name": item_name, "items": [item_name]} for item_name in self.selected_structure_data.keys()]

    def create_selected_data_views(self, tabs_info):
        # Clear existing tabs if needed (keep persistent tabs as required)
        if self.data_view_notebook.index("end") > 0:  # If there are tabs
            if self.data_view_notebook.tab(0, "text") == "No Data":
                self.data_view_notebook.forget(0)  # Remove the "No Data" tab
        
        # Create tabs based on the parsed information
        for tab_info in tabs_info:
            # Create a new frame for the tab
            tab_frame = ttk.Frame(self.data_view_notebook)
            self.data_view_notebook.add(tab_frame, text=tab_info["name"])
            
            if not tab_info["items"]:
                ttk.Label(tab_frame, text="No data selected").pack(expand=True)
                continue
            
            # Create a figure for plotting
            fig = plt.Figure(figsize=(6, 4), dpi=100)
            ax = fig.add_subplot(111)
            
            # Plot data for each item in this tab
            for item_name in tab_info["items"]:
                if item_name in self.selected_structure_data:
                    # Here you would get the actual data and plot it
                    # For demonstration, we'll just plot some random data
                    import numpy as np
                    x = np.linspace(0, 10, 100)
                    y = np.sin(x) + np.random.random(100) * 0.2
                    ax.plot(x, y, label=item_name)
            
            ax.set_xlabel('X Label')
            ax.set_ylabel('Y Label')
            ax.set_title(tab_info["name"])
            ax.legend()
            
            # Create a canvas to display the plot
            canvas = FigureCanvasTkAgg(fig, master=tab_frame)
            canvas.draw()
            canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    def report_selected_data(self):
        if not self.selected_structure_data:
            messagebox.showinfo("No Data", "No data selected to report")
            return
        
        # In a real application, this would generate a report
        # For demonstration, we'll just show a message
        items = ", ".join(self.selected_structure_data.keys())
        messagebox.showinfo("Report Generated", f"Report generated for: {items}")

    def open_analysis_view(self):
        if not self.selected_structure_data:
            messagebox.showinfo("No Data", "No data selected for analysis")
            return
        
        # In a real application, this would open a new window with analysis tools
        # For demonstration, we'll just show a message
        items = ", ".join(self.selected_structure_data.keys())
        messagebox.showinfo("Analysis View", f"Opening analysis view for: {items}")

    def refine_structure(self):
        if not self.selected_structure_data:
            messagebox.showinfo("No Data", "No data selected for structure refinement")
            return
        
        # In a real application, this would refine the structure based on properties
        # For demonstration, we'll just show a message
        messagebox.showinfo("Structure Refinement", "Structure refinement process started")

    def clear_data_view(self):
        # Remove all tabs from the data view notebook
        for i in range(self.data_view_notebook.index("end")):
            self.data_view_notebook.forget(0)
        
        # Add an empty tab
        empty_frame = ttk.Frame(self.data_view_notebook)
        self.data_view_notebook.add(empty_frame, text="No Data")
        ttk.Label(empty_frame, text="Select data from the Structure Data panel to view").pack(expand=True)

    # Menu functions
    def new_project(self):
        # Reset application state
        self.user_defined_structure = "S3HS"
        self.structure_name_var.set(self.user_defined_structure)
        self.parse_user_defined_structure()
        self.clear_data_view()
        messagebox.showinfo("New Project", "New project created")

    def open_project(self):
        filename = filedialog.askopenfilename(
            title="Open Project",
            filetypes=(("JSON files", "*.json"), ("All files", "*.*"))
        )
        if filename:
            # In a real application, this would load the project from a file
            messagebox.showinfo("Open Project", f"Project opened from {filename}")

    def save_project(self):
        filename = filedialog.asksaveasfilename(
            title="Save Project",
            defaultextension=".json",
            filetypes=(("JSON files", "*.json"), ("All files", "*.*"))
        )
        if filename:
            # In a real application, this would save the project to a file
            messagebox.showinfo("Save Project", f"Project saved to {filename}")

    def export_data(self):
        if not self.selected_structure_data:
            messagebox.showinfo("No Data", "No data selected to export")
            return
        
        filename = filedialog.asksaveasfilename(
            title="Export Data",
            defaultextension=".csv",
            filetypes=(("CSV files", "*.csv"), ("All files", "*.*"))
        )
        if filename:
            # In a real application, this would export the data to a file
            messagebox.showinfo("Export Data", f"Data exported to {filename}")

    def optimize_parameters(self):
        messagebox.showinfo("Optimize Parameters", "Parameter optimization started")

    def compare_structures(self):
        messagebox.showinfo("Compare Structures", "Structure comparison tool opened")

    def statistical_analysis(self):
        messagebox.showinfo("Statistical Analysis", "Statistical analysis tool opened")

    def generate_report(self):
        messagebox.showinfo("Generate Report", "Report generation started")

    def export_figures(self):
        messagebox.showinfo("Export Figures", "Figure export dialog opened")


if __name__ == "__main__":
    root = tk.Tk()
    app = DIW_PSPP(root)
    root.mainloop()

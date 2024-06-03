from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import transaction

from interaction.models import ChemokinePartnerPair, ChemokinePartnerInteraction
from structure.models import Structure
from partner.models import Partner

import os
import matplotlib
matplotlib.use('Agg')
import sys
import warnings
import io
import json
import pandas as pd
import numpy as np
from tqdm import tqdm
from pathlib import Path
from urllib.parse import urljoin
from Bio import PDB
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import MDAnalysis as mda
import prolif as plf
from rdkit import Chem
from IPython.display import HTML
import numpy as np
import pandas as pd
import re
from matplotlib.colors import ListedColormap
import matplotlib.pyplot as plt, mpld3
from rdkit import Chem

# custom vis.js network interface, because pyvis implementation
# hasn't been updated in a while and doesn't support some of the
# arguments needed to make this figure
class Network:
    _JS_TEMPLATE = """
        var nodes, edges;
        function drawGraph(_id, nodes, edges, options) {
            var container = document.getElementById(_id);
            nodes = new vis.DataSet(nodes);
            edges = new vis.DataSet(edges);
            var data = {nodes: nodes, edges: edges};
            var network = new vis.Network(container, data, options);
            %(post_initialization)s
            return network;
        }
        nodes = %(nodes)s;
        edges = %(edges)s;
        network = drawGraph('%(div_id)s', nodes, edges, %(options)s);
    """
    _HTML_TEMPLATE = """
        <script type="text/javascript" src="https://unpkg.com/vis-network@9.0.4/dist/vis-network.min.js"></script>
        <link href="https://unpkg.com/vis-network@9.0.4/dist/dist/vis-network.min.css" rel="stylesheet" type="text/css" />
        <style type="text/css">
            body { padding: 0; margin: 0; }
        </style>
        <div id="mynetwork"></div>
        <div id="networklegend"></div>
        <script type="text/javascript">
            %(js)s
        </script>
    """
    def __init__(self):
        self.nodes = []
        self.edges = []
        self.post_initialization = ""

    def add_node(self, _id, **kwargs):
        self.nodes.append({"id": _id, **kwargs})

    def add_edge(self, _from, to, **kwargs):
        self.edges.append({"from": _from, "to": to, **kwargs})
    
    def _get_js(self, width="100%", height="800px", div_id="mynetwork",
        fontsize=20):
        """Returns the JavaScript code to draw the network"""
        self.width = width
        self.height = height
        options = {
            "width": width,
            "height": height,
        }
        options.update(self.options)
        js = self._JS_TEMPLATE % dict(div_id=div_id,
                                      nodes=json.dumps(self.nodes),
                                      edges=json.dumps(self.edges),
                                      options=json.dumps(options),
                                      post_initialization=self.post_initialization,
                                     )
        return js
    
    def _get_html(self, **kwargs):
        """Returns the HTML code to draw the network"""
        return self._HTML_TEMPLATE % dict(js=self._get_js(**kwargs))

    def display(self, **kwargs):
        """Prepare and display the network"""
        html = self._get_html(**kwargs)
        iframe = ('<iframe width="{width}" height="{height}" frameborder="0" '
                  'srcdoc="{doc}"></iframe>')
        return HTML(iframe.format(width=self.width, height=self.height,
                                  doc=escape(html)))
    
    def show(self, filename, **kwargs):
        """Save and display the network"""
        html = self._get_html(**kwargs)
        with open(filename, "w") as f:
            f.write(html)
        iframe = ('<iframe width="{width}" height="{height}" frameborder="0" '
                  'src="{filename}"></iframe>')
        return HTML(iframe.format(width=self.width, height="800px" if "%" in self.height else self.height,
                                  filename=filename))


class Command(BaseBuild):
    help = 'Splits PDB into chemokine and partner files and calculates interactions with ProLIF'

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")

    def handle(self, *args, **options):
        #if options['purge']:
        #    PDBData.objects.all().delete()

        # Calculate interactions from structures
        self.prolif()


    def prolif(self):
        '''Use ProLIF to calculate interactions between chemokine and partner pdb objects.
            Interactions are stored in a table as csv file.'''
        error_log_path = os.path.join(settings.DATA_DIR, 'error_log.txt')
        output_directory = os.path.join(settings.DATA_DIR, 'interactions_output')
        pdb_directory = os.path.join(settings.DATA_DIR, 'prepared_pdbs')

        with open(error_log_path, 'w') as error_log:
            warnings.filterwarnings('ignore')

            for structure in Structure.objects.all():
                pdb_id = structure.pdb_id.lower()

                # Define the directory where chains for this structure are stored
                chain_dir = os.path.join(pdb_directory, f"{pdb_id}_prepared")  # Assuming chains are stored in 'prepared_pdbs/{pdb_id}/'

                # Check if chain directory exists
                if not os.path.isdir(chain_dir):
                    print(f"-----Skipping {pdb_id} as the chain directory does not exist.-----")
                    continue
                try:
                    structure_file = os.path.join(settings.DATA_DIR, 'prepared_pdbs', f"{pdb_id}_prepared.pdb")
                    rdkit_struct = Chem.MolFromPDBFile(structure_file, removeHs=False)
                    struct_mol = plf.Molecule(rdkit_struct)
                    if struct_mol is None:
                        error_message = f"Failed to parse partner chain for {pdb_id}. Skipping.\n"
                        print(error_message)
                        error_log.write(error_message)
                        continue
                    
                    # Calculate interactions
                    fp = plf.Fingerprint()
                    fp.run_from_iterable([struct_mol], struct_mol)

                    # Convert interactions to DataFrame and save as CSV
                    df = fp.to_dataframe()
                    output_csv_path = os.path.join(chain_dir, f"{pdb_id}_all_interactions.csv")
                    df.to_csv(output_csv_path, index=False)
                    print(f"Interactions saved to {output_csv_path}")
                except Exception as e:
                    error_message = f"Error processing {pdb_id}: {e}\n"
                    print(error_message)
                    error_log.write(error_message)


#                    ### Create interaction network plot ###
#                    def sort_ifp_columns(df, order=["ligand", "protein", "interaction"]):
#                        df = df.copy()
#                        cols = np.array([(
#                                plf.ResidueId.from_string(i[0]),
#                                plf.ResidueId.from_string(i[1]),
#                                i[2])
#                            for i in df.columns.values],
#                            dtype=[("ligand", object), ("protein", object), ("interaction", object)])
#                        df.columns = df.columns[np.argsort(cols, order=order)]
#                        df = df.reorder_levels(order, axis=1)
#                        return df
#
#                    df_network = df.copy()
#                    df_network = sort_ifp_columns(df_network)
#                    df_network = df.droplevel("interaction", axis=1)
#
#                    contacts = np.array([list(i) for i in df_network.columns.values])
#                    contacts = np.unique(np.sort(contacts, axis=1), axis=0)
#
#                    net = Network()
#                    colors = {  "CK": "#41ae76",
#                                "Partner": "#3690c0",}
#                    rid_cache = {}
#                    nodes = np.unique(contacts)
#
#                    # Creating nodes
#                    for r in nodes:
#                        rid = plf.ResidueId.from_string(r)
#                        rid_cache[r] = rid
#                        label = f"{rid.name}{rid.number}"
#                        if rid.chain == chemokine_chain_id:
#                            title = "CK"
#                            label = f"CK\n{label}"
#                            color = colors.get("CK", "grey")
#                            net.add_node(r, label=label, title=title, shape="box",
#                                        color=color, dtype="chemokine", margin=20)
#                        else:
#                            if rid.chain == partner_chain_id:
#                                title = "Partner"
#                                color = colors.get("Partner", "grey")
#                                net.add_node(r, label=label, title=title, shape="circle",
#                                            color=color, dtype="partner")
#                            else:
#                                title = "Partner"
#                                label = f"{title}\n{label}"
#                                net.add_node(r, label=label, title=title, shape="circle",
#                                            color="#fec44f", dtype="partner")
#                
#                
#                    def get_interaction_type(dataframe, residue1, residue2):
#                        """
#                        Function to find the interaction type for a specified pair of residues.
#
#                        :param dataframe: The dataframe containing the interaction data.
#                        :param residue1: The residue to match in the first row.
#                        :param residue2: The residue to match in the second row.
#                        :return: The interaction type if a match is found, otherwise a message indicating no match.
#                        """
#                        for col in dataframe.columns:
#                            if dataframe.iloc[0, col] == residue1 and dataframe.iloc[1, col] == residue2:
#                                return dataframe.iloc[2, col]
#                            if dataframe.iloc[0, col] == residue2 and dataframe.iloc[1, col] == residue1:
#                                return dataframe.iloc[2, col]
#                        return None
#
#                    interaction_colors = {
#                        "Hydrophobic": "#50C878",
#                        "PiStacking": "#800080",
#                        "HBDonor": "#FFA500",
#                        "HBAcceptor": "#FFFF00",
#                        "Cationic": "#FF0000",
#                        }
#
#                    # Creating edges
#                    data = pd.read_csv(os.path.join(pdb_output_directory, f"{pdb_id}_interactions.csv"), header=None)
#
#                    for r1, r2 in contacts:
#                        rid1 = rid_cache[r1]
#                        rid2 = rid_cache[r2]
#                        label = f"{r1} - {r2}"
#                        interaction_type = get_interaction_type(data, r1, r2)
#                        edge_color = interaction_colors.get(interaction_type, "#a9a9a9")
#                        same_chain = rid1.chain == rid2.chain
#                        width = 3
#                        length = 100 if same_chain else 300
#                        dashes = [10] if same_chain else False
#
#                        # Add the edge with the specified color
#                        net.add_edge(r1, r2, length=length, title=label, color=edge_color,
#                                    width=width, selectionWidth=6, dashes=dashes)
#
#                    for i in range(len(nodes)):
#                        node = net.nodes[i]
#                        nid = node["id"]
#                        n_neighbors = 0
#                        for edge in net.edges:
#                            if nid in [edge["from"], edge["to"]]:
#                                n_neighbors += 1
#                        node["margin"] = max(5, n_neighbors*3.5)
#
#                    net.options = {
#                        "nodes": {
#                            "font": {
#                                "size": 36,
#                            },
#                        },
#                        "edges": {
#                            "smooth": {
#                                "type": "continuous",
#                            },
#                        },
#                        "physics": {
#                            "hierarchicalRepulsion": {
#                                "avoidOverlap": 0.8,
#                                "springConstant": 0.005,
#                            },
#                            "solver": "hierarchicalRepulsion",
#                            "minVelocity": 3,
#                        },
#                        "interaction": {
#                            "hover": True,
#                            "multiselect": True,
#                        },
#                    }
#                    net.post_initialization = """
#                    network.on("stabilizationIterationsDone", function () {
#                        network.setOptions( { physics: false } );
#                    });
#                    """ 
#                    # Save the interaction network plot in HTML
#                    html_file_path = os.path.join(pdb_output_directory, f"{pdb_id}_network.html")
#                    net.show(html_file_path)
#                
#                
#                            
#                    ### Create interaction matrix plot ###
#                    def extract_number(residue_name):
#                        match = re.search(r'\d+', residue_name)
#                        return int(match.group()) if match else 0
#
#                    # Read the CSV file
#                    interactions_df = data
#                    # Extract unique protein and ligand residues
#                    protein_residues = interactions_df.iloc[0].unique()
#                    ligand_residues = interactions_df.iloc[1].unique()
#
#                    # Sorting the protein and ligand residues
#                    sorted_protein_residues = sorted(protein_residues, key=extract_number)
#                    sorted_ligand_residues = sorted(ligand_residues, key=extract_number)
#
#                    # Update the index mappings
#                    protein_indices = {residue: i for i, residue in enumerate(sorted_protein_residues)}
#                    ligand_indices = {residue: i for i, residue in enumerate(sorted_ligand_residues)}
#
#                    # Identifying unique interaction types
#                    interaction_types = interactions_df.iloc[2].unique()
#
#                    # Specifying HEX color codes for each interaction type
#                    interaction_colors_hex = {
#                        "Hydrophobic": "#50C878",
#                        "PiStacking": "#800080",
#                        "HBDonor": "#FFA500",
#                        "HBAcceptor": "#FFFF00",
#                        "Cationic": "#FF0000",
#                        "Default": "#FFFFFF"
#                    }
#
#                    # Map each interaction type to a color from the specified HEX codes
#                    interaction_colors = {interaction: interaction_colors_hex.get(interaction, interaction_colors_hex["Default"]) for interaction in interaction_types}
#
#                    # Initializing the interaction matrix
#                    interaction_matrix_colored = np.full((len(sorted_ligand_residues), len(sorted_protein_residues)), -1, dtype=int)
#
#                    # Filling the matrix
#                    for col in interactions_df.columns:
#                        protein_residue = interactions_df.iloc[0, col]
#                        ligand_residue = interactions_df.iloc[1, col]
#                        interaction_type = interactions_df.iloc[2, col]
#
#                        if protein_residue in protein_indices and ligand_residue in ligand_indices:
#                            interaction_matrix_colored[ligand_indices[ligand_residue], protein_indices[protein_residue]] = list(interaction_colors.keys()).index(interaction_type)
#
#                    # Creating a colormap
#                    colors_list = [interaction_colors_hex["Default"]] + list(interaction_colors.values())
#                    cmap = ListedColormap(colors_list)
#
#                    # Plotting
#                    plt.figure(figsize=(12, 8))
#                    plt.imshow(interaction_matrix_colored, cmap=cmap, aspect='auto')
#                    cbar = plt.colorbar(ticks=range(-1, len(interaction_colors)), label='Interaction Types')
#                    cbar.ax.set_yticklabels(['No Interaction'] + list(interaction_colors.keys()))
#                    plt.xticks(ticks=np.arange(len(sorted_protein_residues)), labels=sorted_protein_residues, rotation=90)
#                    plt.yticks(ticks=np.arange(len(sorted_ligand_residues)), labels=sorted_ligand_residues)
#                    plt.title(f'Protein-Ligand Interaction Matrix of {pdb_id}')
#                    plt.xlabel('Protein Residues')
#                    plt.ylabel('Ligand Residues')
#                    plt.grid(False)
#
#                    # Save the plot as a PNG file
#                    interaction_matrix_png = os.path.join(settings.DATA_DIR, 'interactions_output', pdb_id, f"{pdb_id}_interaction_matrix.png")
#                    plt.savefig(interaction_matrix_png, bbox_inches='tight')
#
#                    # Close the plot
#                    plt.close()             
#
#                    print(f"-----Interactions succesfully calculated for {pdb_id}-----")
#
#                except Exception as e:
#                    error_message = f"-----An error occurred with {pdb_id}: {e}-----\n"
#                    print(error_message)
#                    error_log.write(error_message)
#                    continue
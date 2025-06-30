import os
import csv
import requests
from django.conf import settings
from django.core.management.base import BaseCommand
from protein.models import GPCRProtein

# Path to the CSV file in the data directory
gpcr_protein_source_file = os.path.join(settings.DATA_DIR, 'protein_data', 'gpcr_protein_data.csv')

# Default values for other GPCR attributes
default_gpcr_data = {
    "gpcr_class": "Class A",
    "family": "Chemokine",
    "endogenous_ligand": "Unknown",
    "g_protein_coupling": "G_i/o",
    "signaling_pathways": {},
}


def fetch_residues_from_gpcrdb(gene_name):
    """Fetches residues data from GPCRdb for a given gene name."""
    url = f"https://gpcrdb.org/services/residues/{gene_name}/"
    headers = {"accept": "application/json"}
    
    try:
        response = requests.get(url, headers=headers)
        response.raise_for_status()
        
        residues = response.json()
        
        # Convert residues list to a dictionary format
        residue_dict = {
            residue["sequence_number"]: {
                "amino_acid": residue["amino_acid"],
                "protein_segment": residue["protein_segment"],
                "generic_number": residue["display_generic_number"]
            }
            for residue in residues
        }
        return residue_dict
    
    except requests.exceptions.RequestException as e:
        print(f"Error fetching residues for gene name {gene_name}: {e}")
        return None


class Command(BaseCommand):
    help = "Populates GPCR proteins in the database with data from GPCRdb using a CSV file."

    def handle(self, *args, **options):
        GPCRProtein.objects.all().delete()
        
        # Verify the file exists before reading
        if not os.path.exists(gpcr_protein_source_file):
            self.stdout.write(self.style.ERROR(f"File not found: {gpcr_protein_source_file}"))
            return
        
        try:
            with open(gpcr_protein_source_file, mode='r', newline='', encoding='utf-8') as file:
                reader = csv.DictReader(file, delimiter=';')
                print("CSV Headers:", reader.fieldnames)  # Print headers to confirm names
                
                for row in reader:
                    common_name = row['name']
                    gene_name = row['gene_name']
                    uniprot_id = row['uniprot_id']
                    
                    # Fetch the residues dictionary from GPCRdb using gene name
                    residues = fetch_residues_from_gpcrdb(gene_name)

                    if residues is None:
                        self.stdout.write(self.style.WARNING(f"Skipping {common_name} due to missing residue data"))
                        continue

                    # Merge default data with specific protein name and UniProt ID
                    gpcr_protein_data = {
                        **default_gpcr_data,
                        "name": common_name,
                        "residues": residues,
                        "unp_accession": uniprot_id,
                    }
                    
                    # Create or update GPCRProtein entries
                    obj, created = GPCRProtein.objects.update_or_create(
                        name=common_name,
                        defaults=gpcr_protein_data,
                    )
                    
                    if created:
                        self.stdout.write(self.style.SUCCESS(f"Created {common_name} ({gene_name})"))
                    else:
                        self.stdout.write(self.style.SUCCESS(f"Updated {common_name} ({gene_name})"))
        
        except csv.Error as e:
            self.stdout.write(self.style.ERROR(f"CSV reading error: {e}"))

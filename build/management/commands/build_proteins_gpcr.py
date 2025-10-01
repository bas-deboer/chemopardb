import os
import csv
import requests
from django.conf import settings
from django.core.management.base import BaseCommand
from protein.models import GPCRProtein
from partner.models import Partner

# Path to the CSV file in the data directory
gpcr_protein_source_file = os.path.join(settings.DATA_DIR, 'protein_data', 'gpcr_protein_data.csv')

# Default values for other GPCR attributes
default_gpcr_data = {
    "gpcr_class": "Class A",
    "family": "Chemokine",
    "endogenous_ligand": "",
    "g_protein_coupling": "",
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
    help = "Load chemokine GPCRs from CSV, fetch GPCRdb residue generic numbers, store GPCR details, and create Partner entries."

    def handle(self, *args, **options):
        # Refresh GPCRProtein table to reflect current CSV + GPCRdb data
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

                    # Also create/update a Partner entry so GPCRs are available as partners in the app
                    try:
                        partner_defaults = {
                            "type": "GPCR",
                            "description": f"UniProt: {uniprot_id}; Gene: {gene_name}; Family: {default_gpcr_data.get('family')}",
                            "gpcrdb_url": f"https://gpcrdb.org/protein/{gene_name}/",
                            "uniprot_id": uniprot_id,
                        }
                        partner, p_created = Partner.objects.update_or_create(
                            name=common_name,
                            defaults=partner_defaults,
                        )
                        if p_created:
                            self.stdout.write(self.style.SUCCESS(f"Partner created for {common_name}"))
                        else:
                            self.stdout.write(self.style.SUCCESS(f"Partner updated for {common_name}"))
                    except Exception as e:
                        self.stdout.write(self.style.WARNING(f"Failed to upsert Partner for {common_name}: {e}"))
        
        except csv.Error as e:
            self.stdout.write(self.style.ERROR(f"CSV reading error: {e}"))

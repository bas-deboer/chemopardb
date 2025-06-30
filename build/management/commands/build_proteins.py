import os
import csv
import json
import logging
import pandas as pd
from urllib.request import urlopen

from tqdm import tqdm
from django.core.management.base import BaseCommand
from django.conf import settings
from django.db import IntegrityError

from Bio import Align

from protein.models import (
    Protein, Species, ProteinFamily, ProteinSequenceType,
    ProteinSource, ProteinAlias, Gene, SignalSequence
)
from residue.models import Residue
from common.models import WebResource, WebLink


class Command(BaseCommand):
    help = 'Imports proteins, filters by ALL_para_df_mod.csv, uses alignment columns as ccn_number, compares UniProt and alignment sequences'

    def __init__(self):
        super().__init__()
        self.protein_source_file = os.path.join(settings.DATA_DIR, 'protein_data', 'pfam_protein_data.csv')
        self.remote_uniprot_dir = 'https://rest.uniprot.org/uniprotkb/'
        self.logger = logging.getLogger(__name__)
        self.accession_filter = set()
        self.alignment_ccn_map = {}
        self.alignment_ccn_columns = []
        self.allpara_sequences = {}

    def add_arguments(self, parser):
        parser.add_argument("--debug", action="store_true", help="Enable debug mode")
        parser.add_argument("--test", action="store_true", help="Use test data")
        parser.add_argument("--purge", action="store_true", help="Purge existing data before import")

    def handle(self, *args, **options):
        if options['debug']:
            self.logger.setLevel(logging.DEBUG)

        if options['test']:
            self.protein_source_file = os.path.join(
                settings.DATA_DIR, 'protein_data', 'proteins_and_families_test.txt'
            )

        if options['purge']:
            self.purge_existing_data()

        self.load_alignment_ccn_map()

        try:
            self.create_parent_protein_family()
            self.create_proteins_and_families()
        except Exception as e:
            self.logger.error("Error during execution: %s", e)

    def purge_existing_data(self):
        self.logger.info("Purging existing data.")
        Protein.objects.all().delete()
        Species.objects.all().delete()
        ProteinFamily.objects.all().delete()

    def create_parent_protein_family(self):
        ProteinFamily.objects.get_or_create(slug='000', defaults={'name': 'Parent family'})

    def load_alignment_ccn_map(self):
        align_file = os.path.join(settings.DATA_DIR, 'protein_data', 'ALL_para_df_mod.csv')
        if not os.path.exists(align_file):
            self.logger.error("Alignment file %s not found!", align_file)
            raise Exception("Alignment file not found")

        df = pd.read_csv(align_file, delimiter=';')
        df.columns = df.columns.str.strip()
        #print("ALL_para_df_mod.csv columns:", list(df.columns))
        #print("First row sample:", df.iloc[0].to_dict())

        # The meta columns are: 'protein', 'uniprot', 'seq', 'class'
        meta_cols = ['protein', 'uniprot', 'seq', 'class']
        self.alignment_ccn_columns = [c for c in df.columns if c not in meta_cols]

        self.accession_filter = set(df['uniprot'].astype(str).str.strip())
        self.logger.info(f"Loaded {len(self.accession_filter)} accessions for filtering from ALL_para_df_mod.csv")

        for _, row in df.iterrows():
            acc = str(row['uniprot']).strip()
            # Build the sequence by concatenating all non-gap, non-null amino acids from the alignment columns
            sequence_letters = [
                str(row[col]).strip() for col in self.alignment_ccn_columns
                if pd.notnull(row[col]) and str(row[col]).strip() not in ('', '-')
            ]
            seq_val = "".join(sequence_letters)
            self.allpara_sequences[acc] = seq_val

            # Now, get the actual residue at each alignment position (None for gaps)
            pos_to_aa = [
                (str(row[col]).strip() if pd.notnull(row[col]) and str(row[col]).strip() not in ('', '-') else None)
                for col in self.alignment_ccn_columns
            ]
            self.alignment_ccn_map[acc] = pos_to_aa


    def create_proteins_and_families(self):
        self.logger.info("Creating proteins and families from file: %s", self.protein_source_file)

        sequence_type, _ = ProteinSequenceType.objects.get_or_create(
            slug='wt', defaults={'name': 'Wild-type'}
        )

        with open(self.protein_source_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            headers = next(reader)

            for row in tqdm(reader, desc="Processing proteins", unit="protein"):
                try:
                    slug = row[0].lower().replace(' ', '_')
                    family, _ = ProteinFamily.objects.get_or_create(slug=slug, defaults={'name': row[0]})
                    species_csv, _ = Species.objects.get_or_create(latin_name=row[2])

                    protein_name = f"{row[1]}_{species_csv.latin_name}"
                    accession = row[3].strip()

                    # Only process proteins in the accession filter
                    if self.accession_filter and accession not in self.accession_filter:
                        continue

                    uniprot_data = self.parse_uniprot_file(accession)
                    if not uniprot_data:
                        continue

                    # --- Sequence comparison ---
                    uniprot_seq = uniprot_data['sequence']
                    allpara_seq = self.allpara_sequences.get(accession, "")

                    if not allpara_seq:
                        self.logger.warning(f"No ALL_para_df_mod sequence found for accession {accession}")
                    else:
                        if uniprot_seq == allpara_seq:
                            self.logger.info(f"Sequences match for {accession} (length: {len(uniprot_seq)})")
                        else:
                            self.logger.warning(f"Sequence mismatch for {accession}: "
                                                f"UniProt length {len(uniprot_seq)}, ALL_para_df_mod length {len(allpara_seq)}")
                            if len(uniprot_seq) < 120 and len(allpara_seq) < 120:
                                self.logger.warning(f"UniProt:      {uniprot_seq}")
                                self.logger.warning(f"ALL_para_df_mod: {allpara_seq}")

                                # Align and print
                                aligner = Align.PairwiseAligner()
                                aligner.mode = 'global'
                                alignments = aligner.align(uniprot_seq, allpara_seq)
                                top_alignment = alignments[0]
                                # Bio.Align.Alignment objects have a nice __str__ for printing
                                print("Alignment between UniProt and ALL_para_df_mod sequences:")
                                print(top_alignment)


                    # Get CCN numbers from the alignment, not by MSA!
                    generic_numbers = self.get_ccn_map_from_alignment(accession, uniprot_seq)


                    self.create_protein_record(
                        protein_name, family, sequence_type, accession,
                        uniprot_data, generic_numbers
                    )

                except Exception as e:
                    self.logger.error("Error processing row %s: %s", row, e)

    def get_ccn_map_from_alignment(self, accession, sequence):
        """
        Returns a dict mapping sequence index (1-based) to CCN label (alignment column header)
        """
        ccn_map = {}
        pos_list = self.alignment_ccn_map.get(accession, [])
        seq_pos = 1  # sequence index, 1-based

        # Go through all alignment columns (which match the order in the CSV after meta columns)
        for align_col_idx, aa in enumerate(pos_list):
            if aa is None:  # This means a gap ('-') for this alignment column in this protein
                continue
            # Only map if the corresponding aligned character is not a gap
            ccn_map[seq_pos] = self.alignment_ccn_columns[align_col_idx]
            seq_pos += 1
            # Stop if we've reached the actual sequence length
            if seq_pos > len(sequence):
                break
        return ccn_map


    def create_protein_record(self, name, family, sequence_type, accession, uniprot_data, generic_numbers):
        try:
            source, _ = ProteinSource.objects.get_or_create(
                name=uniprot_data['source'], defaults={'name': uniprot_data['source']}
            )
            species, _ = Species.objects.get_or_create(
                latin_name=uniprot_data['species_latin_name'],
                defaults={'common_name': uniprot_data['species_common_name']}
            )

            protein = Protein(
                name=name,
                gene_name=name,
                subfamily=family,
                type=name.split("_")[0],
                species=species,
                full_name=uniprot_data.get('full_name', name),
                uniprot_id=accession,
                sequence=uniprot_data['sequence'],
                sequence_type=sequence_type,
                source=source
            )
            protein.save()

            self.create_residues(protein, uniprot_data['sequence'], generic_numbers)
            self.create_protein_aliases(protein, uniprot_data['names'])
            self.create_genes(protein, uniprot_data['genes'], species)

            sig = uniprot_data.get("signal_peptide")
            if sig:
                start, end = sig["start"], sig["end"]
                generic_start = generic_numbers.get(start)
                generic_end = generic_numbers.get(end)
                SignalSequence.objects.create(
                    protein=protein,
                    start=start,
                    end=end,
                    generic_start=generic_start,
                    generic_end=generic_end,
                    sequence=sig["sequence"]
                )
            else:
                self.logger.info("No signal peptide for %s", accession)

        except IntegrityError as e:
            self.logger.error("Failed creating protein %s: %s", name, e)

    def create_residues(self, protein, sequence, generic_numbers):
        amino_acids = {
            'A': 'Ala','R': 'Arg','N': 'Asn','D': 'Asp','C': 'Cys',
            'E': 'Glu','Q': 'Gln','G': 'Gly','H': 'His','I': 'Ile',
            'L': 'Leu','K': 'Lys','M': 'Met','F': 'Phe','P': 'Pro',
            'S': 'Ser','T': 'Thr','W': 'Trp','Y': 'Tyr','V': 'Val'
        }
        residues = []
        for idx, aa in enumerate(sequence, start=1):
            label = generic_numbers.get(idx)
            try:
                as_int = int(label) if label is not None else None
            except ValueError:
                as_int = None

            residues.append(
                Residue(
                    protein=protein,
                    sequence_number=idx,
                    amino_acid=aa,
                    amino_acid_three_letter=amino_acids.get(aa, 'Xxx'),
                    #generic_number=as_int,
                    ccn_number=label
                )
            )
        Residue.objects.bulk_create(residues)

    def create_protein_aliases(self, protein, names):
        aliases = [
            ProteinAlias(protein=protein, name=name, position=i)
            for i, name in enumerate(names)
        ]
        ProteinAlias.objects.bulk_create(aliases)

    def create_genes(self, protein, genes, species):
        for pos, gene in enumerate(genes):
            gene_obj, _ = Gene.objects.get_or_create(
                name=gene, species=species, position=pos
            )
            gene_obj.proteins.add(protein)

    def parse_uniprot_file(self, accession):
        up = {'genes': [], 'names': [], 'accessions': [], 'sequence': '', 'signal_peptide': None}
        url = f"{self.remote_uniprot_dir}{accession}"
        try:
            with urlopen(url) as uf:
                data = json.load(uf)
        except Exception as e:
            self.logger.error("Error reading UniProt JSON for %s: %s", accession, e)
            return None

        up['entry_name'] = data.get("uniProtkbId", "").lower()
        up['source'] = "SWISSPROT" if data.get("entryType", "").lower() == "reviewed" else "TREMBL"
        up['accessions'] = [data.get("primaryAccession", "")] + data.get("secondaryAccessions", [])

        desc = data.get("proteinDescription", {}).get("recommendedName", {})
        full = desc.get("fullName", {}).get("value", "")
        if full:
            up['names'].append(full)
        up['full_name'] = full

        for gene in data.get("genes", []):
            name = gene.get("geneName", {}).get("value")
            if name:
                up['genes'].append(name)
            up['genes'] += [syn.get("value") for syn in gene.get("synonyms", []) if syn.get("value")]

        org = data.get("organism", {})
        up['species_latin_name'] = org.get("scientificName", "")
        up['species_common_name'] = org.get("commonName", up['species_latin_name'])

        seq_obj = data.get("sequence", {})
        up['sequence'] = seq_obj.get("value", "")

        for feature in data.get("features", []):
            if feature.get("type") == "Signal":
                loc = feature.get("location", {})
                start = loc.get("start", {}).get("value")
                end = loc.get("end", {}).get("value")
                if isinstance(start, int) and isinstance(end, int):
                    up['signal_peptide'] = {
                        "start": start,
                        "end": end,
                        "sequence": up['sequence'][start-1:end]
                    }
                    break

        return up

from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError
from django.utils.text import slugify

from protein.models import Protein, ChemokineType, Uniprot, Species, ProteinFamily, ProteinSegment, ProteinSequenceType, ProteinSource, ProteinAlias, Gene
from residue.models import Residue, ResidueNumberingScheme, ResidueGenericNumber
#from common.tools import test_model_updates
from common.models import WebResource, WebLink

import django.apps
import shlex
import os
from urllib.request import urlopen

import pandas as pd
import numpy  as np
import math

import csv


class Command(BaseBuild):
    help = 'Reads source data and creates protein families, proteins, and associated tables'
    
    protein_source_file = os.path.join(settings.DATA_DIR, 'pfam_protein_data.csv')
    remote_uniprot_dir = 'http://www.uniprot.org/uniprot/'
    
    #Setting the variables for the test tracking of the model upadates
    tracker = {}
    all_models = django.apps.apps.get_models()[6:]
    #test_model_updates(all_models, tracker, initialize=True)


    def add_arguments(self, parser):
        parser.add_argument("--debug", default=True, action="store_true", help="Debug mode")
        parser.add_argument("--test", default=False, action="store_true", help="Test mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")
        
        
    def handle(self, *args, **options):
        # test mode
        if options['test']:
            self.protein_source_file = os.sep.join([settings.DATA_DIR, 'protein_data',
                'proteins_and_families_test.txt'])
        
        if options['purge']:
            Protein.objects.all().delete()
            ChemokineType.objects.all().delete()
            Uniprot.objects.all().delete()
            Species.objects.all().delete()
            ProteinFamily.objects.all().delete()
            
        # create parent protein family, 000
        try:
            self.create_parent_protein_family()
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            # print(msg)
            self.logger.error(msg)
            
        # create proteins and families
        try:
            self.create_proteins_and_families()
            test_model_updates(self.all_models, self.tracker, check=True)
        except Exception as msg:
            print(msg)
            self.logger.error(msg)
            
            
    def create_parent_protein_family(self):
        pf = ProteinFamily.objects.get_or_create(slug='000', defaults={
            'name': 'Parent family'})
            
            
    def extract_aligned_sequences(self, file_path):
        """
        Extracts all aligned sequences from a CLUSTAL alignment file, correcting for accessions with extra identifiers.

        Args:
        file_path (str): Path to the CLUSTAL alignment file.

        Returns:
        dict: A dictionary where keys are UniProt accessions (cleaned of any additional identifiers) and values are the corresponding aligned sequences.
        """
        sequences = {}
        with open(file_path, 'r') as file:
            for line in file:
                if line.strip() and not line.startswith('CLUSTAL') and not line.startswith(' '):
                    parts = line.split()
                    if len(parts) > 1:
                        # Attempt to clean up the accession to remove any extra identifiers or ranges
                        cleaned_accession = parts[0].split('/')[0]
                        if cleaned_accession in sequences:
                            sequences[cleaned_accession] += parts[1]
                        else:
                            sequences[cleaned_accession] = parts[1]
        return sequences

    def map_residues_to_positions(self, sequence, alignment):
        """
        Maps each residue in a sequence to its position in the alignment.
        Adjusts the key to start from 1 and increment for each matched residue.
        :param sequence: The original protein sequence.
        :param alignment: The aligned sequence from the MSA file.
        :return: Dictionary with sequential keys starting from 1 and values as positions in the alignment.
        """
        generic_number_mapping = {}
        seq_index = 0  # Index for the original sequence (non-gap)
        position_count = 1  # Starting position for the output keys
        for align_index, residue in enumerate(alignment, start=1):
            if residue != '-':
                if seq_index < len(sequence) and residue == sequence[seq_index]:
                    generic_number_mapping[position_count] = align_index
                    position_count += 1
                    seq_index += 1
        return generic_number_mapping

    def create_residues(self, protein, sequence, generic_numbers):
        amino_acids = {
            'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'E': 'Glu',
            'Q': 'Gln', 'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys',
            'M': 'Met', 'F': 'Phe', 'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp',
            'Y': 'Tyr', 'V': 'Val'
        }
        print(sequence)
        for index, amino_acid in enumerate(sequence, start=1):
            amino_acid_three_letter = amino_acids.get(amino_acid, 'Xxx')
            generic_number = generic_numbers.get(index)
            residue = Residue(
                protein=protein,
                sequence_number=index,
                amino_acid=amino_acid,
                amino_acid_three_letter=amino_acid_three_letter,
                generic_number=generic_number
            )
            try:
                residue.save()
                self.logger.info(f'Created residue {amino_acid} at position {index} with generic number {generic_number}')
            except Exception as e:
                self.logger.error(f'Failed to create residue {amino_acid} at position {index}: {e}')

            
    def create_proteins_and_families(self):
        self.logger.info('CREATING PROTEINS')
        self.logger.info('Parsing file ' + self.protein_source_file)

        # get/create protein sequence type
        # Wild-type for all sequences from source file, isoforms handled separately
        try:
            sequence_type, created = ProteinSequenceType.objects.get_or_create(slug='wt',
                defaults={
                'slug': 'wt',
                'name': 'Wild-type',
                })
            if created:
                self.logger.info('Created protein sequence type Wild-type')
        except:
                self.logger.error('Failed creating protein sequence type Wild-type')
        
        # Read reference MSA and make dictionairy of accessions and aligned sequences
        aligned_sequences = self.extract_aligned_sequences(os.path.join(settings.DATA_DIR, 'chemokines_ref_MSA.aln'))
        
        with open(self.protein_source_file, newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=';')
            # Skip the header or you can use next(reader) if there's a header
            headers = next(reader)  # This reads the header row
            for row in reader:
                try:
                    # Generate a unique slug for the family
                    slug = slugify(row[0])
                    family, created = ProteinFamily.objects.get_or_create(slug=slug, defaults={'name': row[0]})
                    
                    species, created = Species.objects.get_or_create(latin_name=row[2])
                    protein_name = f"{row[1]}_{species.latin_name}"
                    protein_accession = row[3]
                    residue_numbering_scheme = False
                    
                    # parse uniprot file for this protein
                    self.logger.info('Parsing uniprot file for protein ' + protein_name)
                    up = self.parse_uniprot_file(protein_accession)
                    if not up:
                        self.logger.error('Failed parsing uniprot file for protein ' + protein_name + ', skipping')
                        continue

                    # Get aligned sequence for this protein
                    aligned_sequence = aligned_sequences.get(protein_accession)
                    print(aligned_sequence)
                    # Store wt protein residues in the database
                    self.create_protein(protein_name, family, sequence_type, residue_numbering_scheme,
                        protein_accession, up, aligned_sequence)

                except Exception as e:
                    self.logger.error(f"Error processing row {row}: {e}")            
            
            
            # GPCRDB method of storing protein objects
            
            ## family hierarchy is determined by indent
            #spaces_per_indent_level = 4
            #last_indent = 0
            #level_family_counter = [0]
            #parent_family = [0]
            #residue_numbering_scheme = False
#
            #for row in protein_file:
            #    # determine the level of indent
            #    indent = len(row) - len(row.lstrip(' '))
            #    indent = int(indent / spaces_per_indent_level)
#
            #    # has the indent changed
            #    if indent != last_indent:
            #        # did the level increase or decrease?
            #        if indent > last_indent:
            #            parent_family.append(0)
            #            level_family_counter.append(0)
            #        elif indent < last_indent:
            #            for j in range(last_indent-indent):
            #                parent_family.pop()
            #                level_family_counter.pop()
#
            #        last_indent = indent
#
            #    # process the line
            #    # is this a family or protein line?
#
            #    ###########
            #    # family/family type
            #    ###########
            #    create_protein = False
            #    if row.strip().startswith('"'): # protein row
            #        split_row = row.strip().split('","')
            #        family_name = split_row[4]
            #        create_protein = True
            #    else: # family row
            #        # check for residue numbering scheme
            #        split_row = row.strip().split('|')
            #        if len(split_row) > 1:
            #            try:
            #                rns = split_row[1].strip()
            #                residue_numbering_scheme = ResidueNumberingScheme.objects.get(slug=rns)
            #            except:
            #                # abort if residue numbering scheme is not found in db
            #                raise Exception('Residue numbering scheme ' + rns + ' not found, aborting')
            #        else:
            #            if not residue_numbering_scheme:
            #                # abort if no residue numbering scheme is specified in the protein source file
            #                raise Exception('No residue numbering scheme specified in source data, aborting')
#
            #        family_name = split_row[0].strip()
#
            #        # create the protein family
            #        created_family = self.create_protein_family(family_name, indent, parent_family,
            #            level_family_counter)
            #        if created_family:
            #            parent_family = created_family['parent_family']
            #            level_family_counter = created_family['level_family_counter']
            #        else:
            #            continue
#
#
            #    ###########
            #    # protein
            #    ###########
            #    if create_protein:
            #        protein_name = split_row[4]
#
            #        # accession codes for human, mouse and rat receptors (from IU-PHAR)
            #        # accessions = [split_row[15], split_row[31], split_row[23]]
            #        accessions = [split_row[15]]
#
            #        # create a family for this protein
            #        created_family = self.create_protein_family(protein_name, indent, parent_family,
            #            level_family_counter)
            #        if created_family:
            #            pf = created_family['pf']
            #            parent_family = created_family['parent_family']
            #            level_family_counter = created_family['level_family_counter']
            #        else:
            #            continue
#
            #        for protein_accession in accessions:
            #            # skip protein if there is no accession code
            #            if not protein_accession:
            #                self.logger.error('No accession code for protein ' + protein_name + ', skipping')
            #                continue
#
            #            # skip protein if accession code already exists
            #            if Protein.objects.filter(accession=protein_accession).count() > 0:
            #                self.logger.error('Protein accession ' + protein_accession + ' already exists, skipping')
            #                continue
#
            #            # parse uniprot file for this protein
            #            self.logger.info('Parsing uniprot file for protein ' + protein_name)
            #            up = self.parse_uniprot_file(protein_accession)
            #            if not up:
            #                self.logger.error('Failed parsing uniprot file for protein ' + protein_name + ', skipping')
            #                continue
#
            #            self.create_protein(protein_name, pf, sequence_type, residue_numbering_scheme,
            #                protein_accession, up)

        self.logger.info('COMPLETED CREATING PROTEINS')     
            

    def create_protein(self, name, family, sequence_type, residue_numbering_scheme, accession, uniprot, aligned_sequence):
        # get/create protein source
        try:
            source, created = ProteinSource.objects.get_or_create(name=uniprot['source'],
                defaults={'name': uniprot['source']})
            if created:
                self.logger.info('Created protein source ' + source.name)
        except IntegrityError:
            source = ProteinSource.objects.get(name=uniprot['source'])

        # get/create species
        try:
            species, created = Species.objects.get_or_create(latin_name=uniprot['species_latin_name'],
                defaults={
                'common_name': uniprot['species_common_name'],
                })
            if created:
                self.logger.info('Created species ' + species.latin_name)
        except IntegrityError:
            species = Species.objects.get(latin_name=uniprot['species_latin_name'])

        # create protein
        p = Protein()
        p.family = family
        p.species = species
        p.source = source
        p.residue_numbering_scheme = residue_numbering_scheme
        p.sequence_type = sequence_type
        if accession:
            p.accession = accession
        p.entry_name = uniprot['entry_name']
        p.name = name
        p.sequence = uniprot['sequence']

        try:
            p.save()
            self.logger.info('Created protein {}'.format(p.entry_name))
             
            generic_numbers = self.map_residues_to_positions(sequence=p.sequence, alignment=aligned_sequence)
            print(generic_numbers)
            
            self.create_residues(p, uniprot['sequence'], generic_numbers)  # Call to create residues
        except Exception as e:
            self.logger.error('Failed creating protein {} {}'.format(p.entry_name, str(e)))
            print('WARNING:', p.family, p.species, p.source, p.residue_numbering_scheme, p.sequence_type, p.accession, p.entry_name, p.name, p.sequence, e)

        # protein uniprot links
        if accession:
            for ac in uniprot['accessions'][1:]:
                resource = WebResource.objects.get(slug='uniprot')

                # create a link
                link, created = WebLink.objects.get_or_create(web_resource=resource, index=ac)

                # add the link to this protein
                p.web_links.add(link)

        # protein aliases
        for i, alias in enumerate(uniprot['names']):
            a = ProteinAlias()
            a.protein = p
            a.name = alias
            a.position = i

            try:
                a.save()
                self.logger.info('Created protein alias ' + a.name + ' for protein ' + p.name)
            except:
                self.logger.error('Failed creating protein alias ' + a.name + ' for protein ' + p.name)

        # genes
        for i, gene in enumerate(uniprot['genes']):
            g = False
            try:
                g, created = Gene.objects.get_or_create(name=gene, species=species, position=i)
                if created:
                    self.logger.info('Created gene ' + g.name + ' for protein ' + p.name)
            except IntegrityError:
                g = Gene.objects.get(name=gene, species=species, position=i)

            if g:
                g.proteins.add(p)

    def create_protein_family(self, family_name, indent, parent_family, level_family_counter):
        # find the parent family
        if indent == 0:
            try:
                ppf = ProteinFamily.objects.get(parent__isnull=True)
            except ProteinFamily.DoesNotExist:
                raise Exception('Family 000 not found, aborting')
        else:
            parent_family_id = parent_family[indent-1]
            try:
                ppf = ProteinFamily.objects.get(pk=parent_family_id)
            except ProteinFamily.DoesNotExist:
                self.logger.error('Parent family of ' + family_name + ' not found, skipping')
                return False

        # does this family already exists in db?
        try:
            pf = ProteinFamily.objects.get(name=family_name, parent=ppf)
        except ProteinFamily.DoesNotExist:
            # increment the family counter for the current indent level
            level_family_counter[indent] += 1

            # protein family slug
            family_slug = []
            for level in level_family_counter:
                family_slug.append(str(level).zfill(3))
            family_slug = '_'.join(family_slug)

            # create the protein family
            pf = ProteinFamily()
            pf.parent = ppf
            pf.slug = family_slug
            pf.name = family_name
            try:
                pf.save()
                parent_family[indent] = pf.id

                self.logger.info('Created protein family ' + family_name)
            except Exception as msg:
                self.logger.error('Failed creating protein family' + family_name,msg)

        return {
            'pf': pf,
            'parent_family': parent_family,
            'level_family_counter': level_family_counter,
        }

    def parse_uniprot_file(self, accession):
        filename = accession + '.txt'
        remote_file_path = self.remote_uniprot_dir + filename

        up = {}
        up['genes'] = []
        up['names'] = []
        up['accessions'] = []
        up['sequence'] = ''

        read_sequence = False

        # record whether organism has been read
        os_read = False

        try:
            uf = urlopen(remote_file_path)
            self.logger.info('Reading remote file ' + remote_file_path)

            for raw_line in uf:
                line = raw_line.decode('UTF-8')

                if line.startswith('//'):  # end of file
                    break

                if line.startswith('ID'):
                    split_id_line = line.split()
                    up['entry_name'] = split_id_line[1].lower()
                    review_status = split_id_line[2].strip(';')
                    up['source'] = 'SWISSPROT' if review_status == 'Reviewed' else 'TREMBL'

                elif line.startswith('OS') and not os_read:
                    species_full = line[2:].strip().strip('.')
                    species_split = species_full.split('(')
                    up['species_latin_name'] = species_split[0].strip()
                    up['species_common_name'] = species_split[1].strip().strip(')') if len(species_split) > 1 else up['species_latin_name']
                    os_read = True

                elif line.startswith('AC'):
                    sline = line.split()
                    up['accessions'].extend(ac.strip(';') for ac in sline)

                elif line.startswith('DE'):
                    split_de_line = line.split('=')
                    if len(split_de_line) > 1:
                        split_segment = split_de_line[1].split('{')
                        up['names'].append(split_segment[0].strip().strip(';'))

                elif line.startswith('GN'):
                    split_gn_line = line.split(';')
                    for segment in split_gn_line:
                        if '=' in segment:
                            split_segment = segment.split('=')
                            genes = split_segment[1].split(',')
                            for gene_name in genes:
                                split_gene_name = gene_name.split('{')
                                up['genes'].append(split_gene_name[0].strip())

                elif line.startswith('SQ'):
                    read_sequence = True

                elif read_sequence:
                    up['sequence'] += line.strip().replace(' ', '')

            uf.close()  # close the Uniprot file
        except Exception as e:
            self.logger.error('Error reading remote file: {}'.format(e))
            return False

        return up


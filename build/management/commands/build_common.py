from django.core.management.base import BaseCommand
from django.conf import settings
from django.db import IntegrityError
from common.models import WebResource
from residue.models import ResidueNumberingScheme
from protein.models import ProteinSegment

import os
import csv
import shlex
import logging


class Command(BaseCommand):
    help = "Import resources, protein segments, and residue numbering schemes."

    def __init__(self):
        super().__init__()
        # Define paths for data files
        self.resource_source_file = os.path.join(settings.DATA_DIR, 'common_data', 'resources.txt')
        self.segment_source_file = os.path.join(settings.DATA_DIR, 'common_data', 'protein_segments.csv')
        # Set up logging
        self.logger = logging.getLogger(__name__)

    def add_arguments(self, parser):
        parser.add_argument("--debug", action="store_true", help="Enable debug mode")
        parser.add_argument("--purge", action="store_true", help="Purge existing data before import")

    def handle(self, *args, **options):
        if options["debug"]:
            self.logger.setLevel(logging.DEBUG)
        
        if options["purge"]:
            self.logger.info("Purging existing WebResource data.")
            WebResource.objects.all().delete()

        # Create records from source files
        self.create_resources()
        self.create_protein_segments()
        self.create_ccn_positions()
        
    def create_resources(self):
        """Read and create WebResource entries from the resource source file."""
        self.logger.info("Creating resources from file: %s", self.resource_source_file)

        try:
            with open(self.resource_source_file, "r", encoding="UTF-8") as resource_file:
                for row in resource_file:
                    split_row = shlex.split(row)
                    if len(split_row) < 3:
                        self.logger.warning("Skipping malformed line: %s", row)
                        continue

                    slug, name, url = split_row[0], split_row[1], split_row[2]
                    defaults = {"name": name, "url": url}

                    try:
                        wr, created = WebResource.objects.get_or_create(slug=slug, defaults=defaults)
                        if created:
                            self.logger.info("Created resource with slug: %s", slug)
                    except IntegrityError:
                        self.logger.error("Failed to create resource with slug: %s", slug)
        except FileNotFoundError:
            self.logger.error("Resource file not found: %s", self.resource_source_file)
        except Exception as e:
            self.logger.error("An error occurred while creating resources: %s", e)

    def create_protein_segments(self):
        """Read and create ProteinSegment entries from the segment source file."""
        self.logger.info("Creating protein segments from file: %s", self.segment_source_file)

        ProteinSegment.objects.all().delete()

        try:
            with open(self.segment_source_file, newline="") as csvfile:
                ps_reader = csv.reader(csvfile, delimiter=",", quotechar="|")
                for row in ps_reader:
                    if len(row) < 4:
                        self.logger.warning("Skipping malformed line: %s", row)
                        continue

                    slug, category, start, end = row[0], row[1], int(row[2]), int(row[3])

                    try:
                        segment, created = ProteinSegment.objects.get_or_create(
                            slug=slug,
                            defaults={
                                "name": slug,
                                "category": category,
                                "proteinfamily": "CK",
                                "start": start,
                                "end": end
                            },
                        )
                        if created:
                            self.logger.info("Created protein segment: %s", slug)
                    except IntegrityError:
                        self.logger.error("Failed to create protein segment: %s", slug)
        except FileNotFoundError:
            self.logger.error("Segment file not found: %s", self.segment_source_file)
        except Exception as e:
            self.logger.error("An error occurred while creating protein segments: %s", e)

    def create_ccn_positions(self):
        """Read and create ResiduePosition entries from the ccn_numbers csv file."""
        ccn_csv_file = os.path.join(settings.DATA_DIR, 'ccn_numbers.csv')
        self.logger.info("Creating ccn positions from file: %s", ccn_csv_file)

        try:
            from common.models import ResiduePosition
            ResiduePosition.objects.all().delete()

            with open(ccn_csv_file, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    try:
                        position = int(row["number"])
                        ccn_number = row["ccn_number"]
                        ResiduePosition.objects.get_or_create(
                            position=position,
                            ccn_number=ccn_number
                        )
                        self.logger.debug("Created ResiduePosition: %s, %s", position, ccn_number)
                        print("Created ResiduePosition: %s, %s", position, ccn_number)
                    except Exception as e:
                        self.logger.error("Error adding row %s: %s", row, e)
            self.logger.info("Finished creating ccn positions.")
        except FileNotFoundError:
            self.logger.error("CCN positions file not found: %s", ccn_csv_file)
        except Exception as e:
            self.logger.error("An error occurred while creating ccn positions: %s", e)

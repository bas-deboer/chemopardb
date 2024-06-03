from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.db import IntegrityError

from common.models import WebResource
from residue.models import ResidueNumberingScheme
from protein.models import ProteinSegment

import os, csv
import django.apps
import logging
import shlex
import os


class Command(BaseBuild):

    resource_source_file = os.sep.join([settings.DATA_DIR, 'common_data', 'resources.txt'])
    segment_source_file = os.sep.join([settings.DATA_DIR, 'common_data', 'protein_segments.csv'])

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")

    def handle(self, *args, **options):
        if options['purge']:
            WebResource.objects.all().delete()

        self.create_resources()
        self.create_protein_segments()
        self.create_residue_numbering_schemes()
        
    def create_resources(self):
        self.logger.info('CREATING RESOURCES')
        self.logger.info('Parsing file ' + self.resource_source_file)

        with open(self.resource_source_file, "r", encoding='UTF-8') as resource_source_file:
            for row in resource_source_file:
                split_row = shlex.split(row)

                # create resource
                try:
                    defaults = {
                        'name': split_row[1],
                        'url': split_row[2]
                    }

                    wr, created = WebResource.objects.get_or_create(slug=split_row[0], defaults=defaults)

                    if created:
                        self.logger.info('Created resource ' + wr.slug)

                except:
                    self.logger.error('Failed creating resource ' + split_row[0])
                    continue

    def create_protein_segments(self):
        ProteinSegment.objects.all().delete()
        with open(self.segment_source_file, newline='') as csvfile:
            ps_reader = csv.reader(csvfile, delimiter=',', quotechar='|')
            for row in ps_reader:
                start = int(row[2])
                end = int(row[3])
                
                seg, created = ProteinSegment.objects.get_or_create(slug=row[0], name=row[0], category=row[1], proteinfamily='CK', start=start, end=end)
    
    def create_residue_numbering_schemes(self):
        ns, created = ResidueNumberingScheme.objects.get_or_create(parent=None, slug='gn', short_name='GN', name='ChemoPardb generic number')

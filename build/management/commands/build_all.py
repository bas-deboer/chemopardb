from django.core.management.base import BaseCommand, CommandError
from django.core.management import call_command

import datetime


class Command(BaseCommand):
    help = 'Runs all build functions'

    def add_arguments(self, parser):
        parser.add_argument('-p', '--proc',
                            type=int,
                            action='store',
                            dest='proc',
                            default=1,
                            help='Number of processes to run')
        parser.add_argument('-t', '--test',
                            action='store_true',
                            dest='test',
                            default=False,
                            help='Include only a subset of data for testing')
        parser.add_argument('--phase',
                            type=int,
                            action='store',
                            dest='phase',
                            default=None,
                            help='Specify build phase to run (1 or 2, default: None)')

    def handle(self, *args, **options):
        if options['test']:
            print('Running in test mode')

        if options['proc']>4:
            safe_proc_num = 4
        else:
            safe_proc_num = options['proc']

        phase1 = [
            ['clear_cache'],
            ['build_common'],         
            ['build_proteins'],
            ['build_structures'],
            ['build_models'],
            ['build_sequence_alignment'],
            ['build_structural_alignment'],
            ['build_split_models'],
            ['build_partners'],
            ['build_interactions'],
        ]
        phase2 = [
        ]

        if options['phase']:
            if options['phase']==1:
                commands = phase1
            elif options['phase']==2:
                commands = phase2
        else:
            commands = phase1+phase2

        for c in commands:
            print('{} Running {}'.format(
                datetime.datetime.strftime(datetime.datetime.now(), '%Y-%m-%d %H:%M:%S'), c[0]))
            if len(c) == 2:
                call_command(c[0], **c[1])
            elif len(c) == 3:
                call_command(c[0], *c[1], **c[2])
            else:
                call_command(c[0])

        print('{} Build completed'.format(datetime.datetime.strftime(
            datetime.datetime.now(), '%Y-%m-%d %H:%M:%S')))

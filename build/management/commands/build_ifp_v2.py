from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import connection
from django.utils.text import slugify
from django.db import IntegrityError
from django.db import transaction

from common.models import WebLink, WebResource, Publication
from partner.models import PartnerProteinStructure
from protein.models import Protein, ProteinConformation, ProteinState, ProteinSegment
from structure.models import Structure, PdbData, StructureType, Chain, Rotamer, Fragment, Entity, EntityType
from residue.models import Residue
from interaction.models import *
from django.core.management.base import BaseCommand
from django.db.models import Count
from interaction.models import ChemokinePartnerPair, ChemokinePartnerInteraction, ChemokinePartnerIFP

class Command(BaseCommand):
    help = 'Generates interaction fingerprint bitstrings for each chemokine partner pair based on generic residue numbers'

    def handle(self, *args, **kwargs):
        self.stdout.write("Generating interaction fingerprints...")
        ChemokinePartnerIFP.objects.all().delete()
        
        pairs = ChemokinePartnerPair.objects.all()

        for pair in pairs:
            interaction_data = ChemokinePartnerInteraction.objects.filter(
                chemokine_partner_pair=pair
            ).select_related('chemokine_residue')

            interaction_dict = self.organize_interactions_by_generic_number(interaction_data)
            ifp_string = self.generate_interaction_fingerprint(interaction_dict)
            
            ChemokinePartnerIFP.objects.create(
                chemokine_partner_pair=pair,
                ifp_string=ifp_string,
            )

        self.stdout.write("Interaction fingerprints generation completed.")

    def organize_interactions_by_generic_number(self, interactions):
        interaction_dict = {}

        for interaction in interactions:
            generic_number = interaction.chemokine_residue.generic_number
            if generic_number not in interaction_dict:
                interaction_dict[generic_number] = set()

            interaction_dict[generic_number].add(interaction.interaction_type)

        return interaction_dict

    def generate_interaction_fingerprint(self, interaction_dict):
        # 9 possible interactions
        possible_interactions = [
            "Hydrophobic", "Pistacking", "HBDonor", "HBAcceptor",
            "Anionic", "Cationic", "CationPi", "PiCation", "VdWContact"
        ]

        bitstring_parts = []

        # Loop through generic numbers 30 to 150
        for generic_number in range(30, 151):
            interactions = interaction_dict.get(generic_number, set())
            bitstring = ["0"] * len(possible_interactions)
            for interaction in interactions:
                if interaction in possible_interactions:
                    index = possible_interactions.index(interaction)
                    bitstring[index] = "1"
            bitstring_parts.append("".join(bitstring))

        # Concatenate all bitstrings for each generic number
        return "".join(bitstring_parts)

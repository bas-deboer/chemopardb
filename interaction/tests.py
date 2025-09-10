from django.test import TestCase, RequestFactory
from django.urls import reverse

from common.models import WebResource, WebLink, ResiduePosition
from structure.models import Structure, PdbData, Rotamer, ChemokineBindingPartner
from interaction.models import ChemokinePartnerInteraction
from interaction.views import ResidueSearchView
from residue.models import Residue


class ResidueSearchViewTests(TestCase):
    def setUp(self):
        ResiduePosition.objects.create(position=1, ccn_number="B1.10")
        ResiduePosition.objects.create(position=2, ccn_number="B2.20")

        wr = WebResource.objects.create(slug="pdb", name="PDB")
        wl = WebLink.objects.create(web_resource=wr, index="1ABC")
        self.structure = Structure.objects.create(pdb_code=wl)
        pdbdata = PdbData.objects.create(pdb="test")

        res1 = Residue.objects.create(
            protein=None,
            sequence_number=1,
            ccn_number="B1.10",
            amino_acid="A",
            amino_acid_three_letter="ALA",
        )
        rot1 = Rotamer.objects.create(
            residue=res1,
            structure=self.structure,
            pdbdata=pdbdata,
            chain="A",
            sequence_number=1,
            ccn_number="B1.10",
            amino_acid="A",
            amino_acid_three_letter="ALA",
        )
        res2 = Residue.objects.create(
            protein=None,
            sequence_number=2,
            ccn_number="B2.20",
            amino_acid="G",
            amino_acid_three_letter="GLY",
        )
        rot2 = Rotamer.objects.create(
            residue=res2,
            structure=self.structure,
            pdbdata=pdbdata,
            chain="A",
            sequence_number=2,
            ccn_number="B2.20",
            amino_acid="G",
            amino_acid_three_letter="GLY",
        )

        self.bp1 = ChemokineBindingPartner.objects.create(
            structure=self.structure,
            chemokine_chain="A",
            partner_chain="B",
            chemokine_name="CK1",
            partner_name="P1",
        )
        self.bp2 = ChemokineBindingPartner.objects.create(
            structure=self.structure,
            chemokine_chain="A",
            partner_chain="C",
            chemokine_name="CK1",
            partner_name="P2",
        )
        ChemokinePartnerInteraction.objects.create(
            chemokine_residue=rot1,
            partner_residue="X",
            partner_name="P1",
            partner_chain="B",
            interaction_type="contact",
            structure=self.structure,
            chemokine_binding_partner=self.bp1,
        )
        ChemokinePartnerInteraction.objects.create(
            chemokine_residue=rot2,
            partner_residue="Y",
            partner_name="P2",
            partner_chain="C",
            interaction_type="contact",
            structure=self.structure,
            chemokine_binding_partner=self.bp2,
        )

    def test_filter_by_ccn_number(self):
        factory = RequestFactory()
        url = reverse("interaction:residue_search")
        request = factory.get(url, {"ccn_numbers": ["B1.10"]})
        response = ResidueSearchView.as_view()(request)
        context = response.context_data
        partners = list(context.get("binding_partners"))
        self.assertEqual(partners, [self.bp1])

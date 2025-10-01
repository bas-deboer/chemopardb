from django.db import models


class ChemokinePartnerInteraction(models.Model):
    """Residue-level interactions for experimental chemokine–partner complexes."""
    
    chemokine_residue = models.ForeignKey('structure.Rotamer', on_delete=models.CASCADE, null=True)
    partner_residue = models.CharField(max_length=10, null=True)
    partner_name = models.CharField(max_length=100, null=True)
    partner_chain = models.CharField(max_length=3, null=True)
    interaction_type = models.CharField(max_length=100, null=True)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    chemokine_binding_partner = models.ForeignKey('structure.ChemokineBindingPartner', on_delete=models.CASCADE, null=True)

    class Meta:
        db_table = 'interaction_chemokine_partner_interaction'
        

class PredictedChemokinePartnerInteraction(models.Model):
    """Residue-level interactions for predicted chemokine–GPCR complexes."""

    predicted_complex = models.ForeignKey(
        'structure.PredictedComplex',
        on_delete=models.CASCADE,
        related_name='predicted_interactions',
    )
    chemokine_residue = models.ForeignKey(
        'residue.Residue',
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name='predicted_complex_interactions',
    )
    chemokine_chain = models.CharField(max_length=10, null=True, blank=True)
    chemokine_residue_number = models.IntegerField(null=True, blank=True)
    chemokine_residue_name = models.CharField(max_length=5, null=True, blank=True)
    partner_chain = models.CharField(max_length=10, null=True, blank=True)
    partner_residue_number = models.IntegerField(null=True, blank=True)
    partner_residue_name = models.CharField(max_length=5, null=True, blank=True)
    partner_generic_number = models.CharField(max_length=16, null=True, blank=True)
    interaction_type = models.CharField(max_length=100)
    created = models.DateTimeField(auto_now_add=True)

    class Meta:
        db_table = 'interaction_predicted_chemokine_partner'
        indexes = [
            models.Index(fields=['predicted_complex', 'interaction_type']),
            models.Index(fields=['predicted_complex', 'partner_generic_number']),
        ]
        unique_together = (
            (
                'predicted_complex',
                'chemokine_chain',
                'chemokine_residue_number',
                'partner_chain',
                'partner_residue_number',
                'interaction_type',
            ),
        )

    def __str__(self):
        chem_id = f"{self.chemokine_residue_name or '?'} {self.chemokine_residue_number or '?'}"
        partner_id = f"{self.partner_residue_name or '?'} {self.partner_residue_number or '?'}"
        generic = f" [{self.partner_generic_number}]" if self.partner_generic_number else ""
        return (
            f"Predicted interaction ({self.interaction_type})"
            f" {chem_id.strip()}->{partner_id.strip()}{generic}"
        )


class ResidueFragmentInteraction(models.Model):
    structure_partner_pair = models.ForeignKey('StructurePartnerInteraction', on_delete=models.CASCADE)
    rotamer = models.ForeignKey('structure.Rotamer', on_delete=models.CASCADE)
    fragment = models.ForeignKey('structure.Fragment', on_delete=models.CASCADE)
    #interaction_type = models.ForeignKey('ResidueFragmentInteractionType', on_delete=models.CASCADE)

    def __str__(self):
        if self.rotamer.residue.display_generic_number is not None:
            return "{!s} {!s} {!s}".format(self.structure_partner_pair.structure.pdb_code.index, self.rotamer.residue.display_generic_number.label, self.structure_partner_pair.partner.name)
        else:
            return "{!s} {!s} {!s}".format(self.structure_partner_pair.structure.pdb_code.index, self.rotamer.residue, self.structure_partner_pair.partner.name)
    class Meta():
        db_table = 'interaction_residue_fragment'

    def get_pdbdata(self):
        return "{!s}\n{!s}".format(self.rotamer.pdbdata, self.fragment.pdbdata)


    def generate_filename(self):

        if self.rotamer.residue.display_generic_number is not None:
            generic_num = self.rotamer.residue.display_generic_number.label
        else:
            generic_num = self.rotamer.residue.sequence_number
        res_name = self.rotamer.residue.amino_acid
        prot_entry_name = str(self.structure_partner_pair.structure.protein_conformation.protein.parent.entry_name)
        pdb_code = self.structure_partner_pair.structure.pdb_code.index
        interaction = self.interaction_type.slug

        return "{}_{}_{}_{}_{}.pdb".format(generic_num.replace('.','_'), res_name, prot_entry_name, pdb_code, interaction)

        
class ChemokinePartnerIFP(models.Model):
    structure = models.ForeignKey('structure.Structure', null=True, on_delete=models.CASCADE)
    ifp_string = models.CharField(max_length=2000, null=True)
    binding_pair = models.ForeignKey('structure.ChemokineBindingPartner', null=True, on_delete=models.CASCADE)

    class Meta:
        db_table = 'interaction_chemokine_partner_IFP'
    
        
class StructurePartnerInteraction(models.Model):
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE, null=True)
    partner = models.ForeignKey('partner.Partner', on_delete=models.CASCADE)
    pdb_reference = models.CharField(max_length=3, null=True)
    pdb_file = models.ForeignKey('structure.PdbData', null=True, on_delete=models.CASCADE)
    annotated = models.BooleanField(default=False)

    def __str__(self):
        return "{} {}".format(self.structure.pdb_code, self.partner.name)

    class Meta():
        db_table = 'interaction_structure_partner'
        
        
class LigandNetworkHTML(models.Model):
    ligand_name = models.CharField(max_length=255)
    html_content = models.TextField()
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE, related_name='ligand_network_htmls')

    class Meta:
        db_table = 'ligand_network_html'

    def __str__(self):
        return self.ligand_name


class ChemokineGpcrGtPBinding(models.Model):
    """
    Stores curated bindings between chemokines and GPCRs based on Guide to Pharmacology (GtP).
    Minimal schema linking a chemokine Protein to a GPCR Partner with optional GtP metadata.
    """
    chemokine = models.ForeignKey('protein.Protein', on_delete=models.CASCADE, related_name='gtp_gpcr_bindings')
    gpcr_partner = models.ForeignKey('partner.Partner', on_delete=models.CASCADE, related_name='gtp_chemokine_bindings')

    chemokine_gtp_id = models.CharField(max_length=32, null=True, blank=True)
    gpcr_gtp_target_id = models.CharField(max_length=32, null=True, blank=True)

    type = models.CharField(max_length=64, null=True, blank=True, help_text='e.g., agonist/antagonist/modulator')
    action = models.CharField(max_length=64, null=True, blank=True, help_text='e.g., agonist/antagonist/modulator')
    notes = models.TextField(null=True, blank=True)
    source_url = models.URLField(null=True, blank=True)

    # Affinity-related fields (as provided by GtP payload)
    affinity = models.CharField(max_length=128, null=True, blank=True)
    affinity_parameter = models.CharField(max_length=32, null=True, blank=True)
    original_affinity = models.CharField(max_length=128, null=True, blank=True)
    original_affinity_type = models.CharField(max_length=64, null=True, blank=True)
    original_affinity_relation = models.CharField(max_length=8, null=True, blank=True)
    concentration_range = models.CharField(max_length=64, null=True, blank=True)

    class Meta:
        db_table = 'interaction_chemokine_gpcr_gtp_binding'
        unique_together = ('chemokine', 'gpcr_partner')
        verbose_name = 'Chemokine–GPCR GtP Binding'
        verbose_name_plural = 'Chemokine–GPCR GtP Bindings'

    def __str__(self):
        return f"{self.chemokine.name} ↔ {self.gpcr_partner.name} (GtP)"

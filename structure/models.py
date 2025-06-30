from django.db import models
from django.db.models import JSONField
from common.diagrams_chemokine import *
from django.utils.html import escape


class Structure(models.Model):
    STATE_CHOICES = [
        ('monomer', 'Monomer'),
        ('dimer', 'Dimer'),
        ('tetramer', 'Tetramer'),
        ('polymer', 'Polymer'),
        (None, 'None'),
    ]
    
    structure_type = models.ForeignKey('StructureType', null=True, on_delete=models.CASCADE)
    pdb_code = models.ForeignKey('common.WebLink', on_delete=models.CASCADE, null=True)
    resolution =  models.CharField(max_length=50, null=True)
    protein = models.ForeignKey('protein.Protein', on_delete=models.SET_NULL, related_name='structures', null=True)
    state = models.CharField(max_length=50, choices=STATE_CHOICES, null=True, blank=True)
    publication = models.ForeignKey('common.Publication', null=True, on_delete=models.CASCADE)
    publication_date = models.DateField(null=True)
    pdb_data = models.ForeignKey('PdbData', null=True, on_delete=models.CASCADE)
    chain_id = models.CharField(max_length=100, null=True, blank=True)  # New field to store chain ID
    
    def __str__(self):
        if self.pdb_code:
            return f"{self.pdb_code.index}"

    class Meta():
        db_table = 'structure'
        
    def get_snake_plot(self, chain=None, nobuttons=False):
        try:
            if chain:
                rotamer_list = Rotamer.objects.filter(structure=self, chain=chain)
            else:
                rotamer_list = Rotamer.objects.filter(structure=self)

            if not rotamer_list.exists():
                return "<p>No snake plot available: no rotamer data for this structure.</p>"

            snake_plot = DrawArrestinPlot(rotamer_list, str(self), nobuttons=nobuttons)

            return str(snake_plot)

        except Exception as e:
            # Log the actual error somewhere in production
            return f"<p>Error generating snake plot</p>"


class Entity(models.Model):
    structure = models.ForeignKey('Structure', null=True, on_delete=models.CASCADE)
    entity_type = models.ForeignKey('EntityType', null=True, on_delete=models.CASCADE)
    name = models.CharField(max_length=1000, null=True)
    chain = models.CharField(max_length=200, null=True)
    residues = JSONField(null=True, blank=True)
    organism = models.CharField(max_length=100, null=True)
    unp_accession = models.CharField(max_length=50, null=True, blank=True)
    pfam_accession = models.CharField(max_length=50, null=True, blank=True)
    pubchem_id = models.CharField(max_length=255, null=True, blank=True)
    comp_id = models.CharField(max_length=50, null=True, blank=True)
    chembl_id = models.CharField(max_length=255, null=True, blank=True)
    smiles = models.CharField(max_length=500, null=True, blank=True)
    inchikey = models.CharField(max_length=500, null=True, blank=True)

    def __str__(self):
        return f"{self.name} ({self.chain})"

    class Meta:
        db_table = 'structure_entity'


class ChemokineBindingPartner(models.Model):
    """
    Captures a binding pair between a chemokine and its partner within a given structure.
    """
    
    structure = models.ForeignKey('Structure', on_delete=models.CASCADE, related_name='chemokine_binding_partners', help_text="The structure in which the binding pair is observed.", null=True)
    partner = models.ForeignKey('partner.Partner', null=True, blank=True, on_delete=models.SET_NULL,
                                related_name="binding_instances",
                                help_text="Link to the partner object if defined in the Partner table.")
    chemokine_entity = models.ForeignKey('Entity', on_delete=models.CASCADE, related_name='chemokine_bindings', help_text="The chemokine entity.", null=True)
    chemokine_name = models.CharField(max_length=100, null=True)
    partner_entity = models.ForeignKey('Entity', on_delete=models.CASCADE, related_name='binding_partner_entities', help_text="The binding partner entity.", null=True)
    chemokine_chain = models.CharField(max_length=10, help_text="Chain identifier for the chemokine within the structure.", null=True)
    partner_name = models.CharField(max_length=100, null=True)
    partner_type = models.CharField(max_length=100, null=True)
    partner_chain = models.CharField(max_length=10, help_text="Chain identifier for the binding partner within the structure.", null=True)
    partner_residues = models.JSONField(null=True, blank=True)
    pdb_data = models.ForeignKey('PdbData', null=True, on_delete=models.CASCADE)

    def __str__(self):
        return (f"{self.chemokine_entity} (chain {self.chemokine_chain}) "
                f"– {self.partner_entity} (chain {self.partner_chain}) in {self.structure}")

    class Meta:
        verbose_name = "Chemokine Binding Partner Pair"
        verbose_name_plural = "Chemokine Binding Partner Pairs"
        unique_together = (
            'structure',
            'chemokine_entity',
            'partner_entity',
            'chemokine_chain',
            'partner_chain'
        )


class EntityInstance(models.Model):
    structure = models.ForeignKey(Structure, on_delete=models.CASCADE)
    chain_id = models.CharField(max_length=1)  # Chain ID of the sample structure
    rmsd = models.FloatField()
    rotation_matrix = models.JSONField()  # JSON field to store the rotation matrix
    translation_vector = models.JSONField()  # JSON field to store the translation vector

    def __str__(self):
        return f"EntityInstance: {self.structure.pdb_code.index} Chain {self.chain_id} (RMSD: {self.rmsd})"


class EntityType(models.Model):
    slug = models.SlugField(max_length=100)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'structure_entity_type' 


class StructureType(models.Model):
    slug = models.SlugField(max_length=100)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'structure_type' 


class PdbData(models.Model):
    pdb = models.TextField(max_length=1000000000)

    def __str__(self):
        return self.pdb

    class Meta():
        db_table = "structure_pdb_data"


class Rotamer(models.Model):
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE, null=True)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    pdbdata = models.ForeignKey('PdbData', on_delete=models.CASCADE)
    pdbseq_number = models.IntegerField(null=True)
    chain = models.CharField(max_length=1, null=True)
    
    # Residue information fields
    sequence_number = models.SmallIntegerField(null=True)
    generic_number = models.SmallIntegerField(null=True)
    ccn_number = models.CharField(max_length=10, null=True, blank=True)
    amino_acid = models.CharField(max_length=1, null=True)
    amino_acid_three_letter = models.CharField(max_length=3, null=True)
    residue_type = models.CharField(max_length=10, choices=[('canonical', 'Canonical'), ('mutated', 'Mutated')], null=True)
    segment = models.CharField(max_length=20, null=True)

    def __str__(self):
        return '{} {}{} (Chain: {})'.format(self.structure.pdb_code.index, self.residue.amino_acid, self.residue.sequence_number, self.chain)

    class Meta:
        db_table = "structure_rotamer"



class Fragment(models.Model):
    residue = models.ForeignKey('residue.Residue', on_delete=models.CASCADE)
    #ligand = models.ForeignKey('ligand.Ligand', on_delete=models.CASCADE)
    structure = models.ForeignKey('structure.Structure', on_delete=models.CASCADE)
    pdbdata = models.ForeignKey('PdbData', on_delete=models.CASCADE)
    # TODO
    # Make Ligand models
    def __str__(self):
        return '{} {}{} {}'.format(self.structure.pdb_code.index, self.residue.amino_acid,
            self.residue.sequence_number, self.partner.name)

    class Meta():
        db_table = "structure_fragment"
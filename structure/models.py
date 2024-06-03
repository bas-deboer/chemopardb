from django.db import models
from common.diagrams_chemokine import *


class Structure(models.Model):
    structure_type = models.ForeignKey('StructureType', null=True, on_delete=models.CASCADE)
    pdb_code = models.ForeignKey('common.WebLink', on_delete=models.CASCADE, null=True)
    resolution =  models.CharField(max_length=50, null=True)
    protein = models.ForeignKey('protein.Protein', on_delete=models.SET_NULL, related_name='structures', null=True)
    partner = models.ManyToManyField('partner.Partner')
    state = models.ForeignKey('protein.ProteinState', on_delete=models.CASCADE, null=True)
    publication = models.ForeignKey('common.Publication', null=True, on_delete=models.CASCADE)
    publication_date = models.DateField(null=True)
    pdb_data = models.ForeignKey('PdbData', null=True, on_delete=models.CASCADE)

    def __str__(self):
        if self.pdb_code:
            return f"{self.pdb_code.index}"

    class Meta():
        db_table = 'structure'
        
    def get_snake_plot(self):
        residuelist = Residue.objects.filter(protein=self.protein)
        print(len(residuelist))
        return DrawArrestinPlot(residuelist, str(self))


class Entity(models.Model):
    structure = models.ForeignKey('Structure', null=True, on_delete=models.CASCADE)
    entity_type = models.ForeignKey('EntityType', null=True, on_delete=models.CASCADE)
    name = models.CharField(max_length=100, null=True)
    chain = models.CharField(max_length=5, null=True)
    organism = models.CharField(max_length=100, null=True)
    unp_accession = models.CharField(max_length=50, null=True)
    pfam_accession = models.CharField(max_length=50, null=True)

    def __str__(self):
        return f"{self.name} ({self.chain})"

    class Meta:
        db_table = 'structure_entity'


class EntityType(models.Model):
    slug = models.SlugField()
    name = models.CharField(max_length=20)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'structure_entity_type' 


class Chain(models.Model):
    entity = models.ForeignKey(Entity, null=True, on_delete=models.CASCADE, related_name='chains')
    chain_id = models.CharField(max_length=2, null=True)
    sequence = models.TextField(null=True, blank=True)
    length = models.IntegerField(null=True, blank=True)

    def __str__(self):
        return f"{self.entity.structure.pdb_id} - {self.chain_id} Chain"


#class Chain(models.Model):
#    structure = models.ForeignKey('Structure', on_delete=models.CASCADE)
#    chain = models.CharField(max_length=100, null=True)
#    type = models.CharField(max_length=100, null=True)
#    pdbx_accession = models.CharField(max_length=100, null=True)
#    pfam_accession = models.CharField(max_length=100, null=True)
#    name = models.CharField(max_length=100, null=True)
#    is_chemokine = models.BooleanField(null=True)
#    contacts_chemokine = models.BooleanField(null=True)
#    
#    def __str__(self):
#        return self.chain_ID
#
#    class Meta():
#        db_table = 'structure_chain'


class StructureType(models.Model):
    slug = models.SlugField()
    name = models.CharField(max_length=20)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'structure_type' 


class PdbData(models.Model):
    pdb = models.TextField()

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
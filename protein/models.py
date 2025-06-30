from common.diagrams_chemokine import *
from django.db import models
from residue.models import Residue


class Protein(models.Model):
    name = models.CharField(max_length=50)
    gene_name = models.CharField(max_length=50)
    subfamily = models.CharField(max_length=50)
    type = models.CharField(max_length=50)
    species = models.CharField(max_length=200)
    full_name = models.CharField(max_length=200)
    uniprot_id = models.CharField(max_length=100, db_index=True, null=True)
    iuphar = models.CharField(max_length=50)
    sequence = models.CharField(max_length=1000, null=True)
    sequence_type = models.ForeignKey('ProteinSequenceType', on_delete=models.CASCADE, null=True)
    source = models.ForeignKey('ProteinSource', on_delete=models.CASCADE, null=True)
    web_links = models.ManyToManyField('common.WebLink')
    
    def __str__(self):
    	return self.name

    def get_snake_plot(self):
        residuelist = Residue.objects.filter(protein=self)

        # Access the signal sequence directly via the one-to-one relation.
        try:
            signal_seq = self.signal_sequence
            signal_start = signal_seq.start
            signal_end = signal_seq.end
        except SignalSequence.DoesNotExist:
            print("Signal sequence does not exist for this protein.")
            signal_start = None
            signal_end = None

        return DrawArrestinPlot(residuelist, str(self), signal_start, signal_end, nobuttons='chemokine')


    class Meta():
        db_table = 'protein'


class ProteinSegment(models.Model):
    slug = models.SlugField(max_length=100)
    name = models.CharField(max_length=50)
    category = models.CharField(max_length=50)
    start = models.IntegerField(null=True)
    end = models.IntegerField(null=True)
    proteinfamily = models.CharField(max_length=20)

    def __str__(self):
        return self.slug

    class Meta():
        ordering = ('id', )
        db_table = 'protein_segment'
        unique_together = ('slug', 'proteinfamily')


class Species(models.Model):
    latin_name = models.CharField(max_length=200, unique=False)
    common_name = models.CharField(max_length=200)

    def __str__(self):
        return self.latin_name

    class Meta():
        db_table = 'species'


class ProteinFamily(models.Model):
    parent = models.ForeignKey('self', null=True, on_delete=models.CASCADE)
    slug = models.SlugField(max_length=100, unique=True)
    name = models.CharField(max_length=200)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_family'


class SignalSequence(models.Model):
    protein = models.OneToOneField(
        'Protein',
        related_name='signal_sequence',
        on_delete=models.CASCADE
    )
    start = models.PositiveIntegerField()
    end = models.PositiveIntegerField()
    generic_start = models.CharField(
        max_length=20, 
        blank=True, 
        null=True,
        help_text="Generic numbering for the signal sequence start (e.g., 'SS1')"
    )
    generic_end = models.CharField(
        max_length=20, 
        blank=True, 
        null=True,
        help_text="Generic numbering for the signal sequence end (e.g., 'SS21')"
    )
    sequence = models.TextField()

    class Meta:
        verbose_name = "Signal Sequence"
        verbose_name_plural = "Signal Sequences"
        ordering = ['protein']

    def __str__(self):
        return f"{self.protein.name} signal: {self.start}-{self.end}"


 
class ProteinSequenceType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.slug

    class Meta():
        db_table = 'protein_sequence_type'
        
        
class ProteinSource(models.Model):
    name = models.CharField(max_length=20, unique=True)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'protein_source'
        
        
class ProteinAlias(models.Model):
    protein = models.ForeignKey('Protein', on_delete=models.CASCADE)
    name = models.CharField(max_length=200)
    position = models.SmallIntegerField()

    def __str__(self):
        return self.name

    class Meta():
        ordering = ('position', )
        db_table = 'protein_alias'


class Gene(models.Model):
    proteins = models.ManyToManyField('Protein', related_name='genes')
    species = models.ForeignKey('Species', on_delete=models.CASCADE)
    name = models.CharField(max_length=100)
    position = models.SmallIntegerField()

    def __str__(self):
        return self.name

    class Meta():
        ordering = ('position', )
        unique_together = ('name', 'species','position')
        db_table = 'protein_gene'

        
class GPCRProtein(models.Model):
    name = models.CharField(max_length=200)
    gpcr_class = models.CharField(max_length=50, null=True, blank=True)  # e.g., Class A, B, C
    family = models.CharField(max_length=100, null=True, blank=True)  # Specific GPCR family (e.g., Chemokine)
    endogenous_ligand = models.CharField(max_length=200, null=True, blank=True)  # Endogenous ligand if known
    g_protein_coupling = models.CharField(max_length=100, null=True, blank=True)  # Type of G-protein coupling
    signaling_pathways = models.JSONField(null=True, blank=True)  # JSON field to store pathways information
    residues = models.JSONField(null=True, blank=True)  # JSON field to store GPCR-specific residues
    unp_accession = models.CharField(max_length=100, db_index=True, null=True)

    def __str__(self):
        return f"{self.name} - GPCR Details"

    class Meta:
        db_table = 'gpcr_protein'
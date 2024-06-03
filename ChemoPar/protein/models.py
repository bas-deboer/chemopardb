from common.diagrams_chemokine import *
from django.db import models


class Protein(models.Model):
    name = models.CharField(max_length=200)
    uniprot_id = models.CharField(max_length=100, db_index=True, null=True)
    family = models.CharField(max_length=200)
    species = models.CharField(max_length=200)
    chemokinetype = models.CharField(max_length=200)
    sequence = models.CharField(max_length=1000, null=True)
    sequence_type = models.ForeignKey('ProteinSequenceType', on_delete=models.CASCADE, null=True)
    source = models.ForeignKey('ProteinSource', on_delete=models.CASCADE, null=True)
    web_links = models.ManyToManyField('common.WebLink')
    accession = models.CharField(max_length=100, db_index=True, null=True)
    
    def __str__(self):
    	return self.name

    class Meta():
        db_table = 'protein'
        
    def get_snake_plot(self):
        residuelist = Residue.objects.filter(protein=self)
        print(len(residuelist))
        return DrawArrestinPlot(residuelist, str(self))


class ProteinConformation(models.Model):
    protein = models.ForeignKey('Protein', on_delete=models.CASCADE)
    state = models.ForeignKey('ProteinState', on_delete=models.CASCADE)
    #protein_anomalies = models.ManyToManyField('protein.ProteinAnomaly')

    def __str__(self):
        return '<{}-{}>'.format(self.domain.isoform.protein.entry_name, self.state.slug)

    class Meta():
        ordering = ('id', )
        db_table = "protein_conformation"


class ProteinState(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.slug

    class Meta():
        db_table = "protein_state"


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


class Uniprot(models.Model):
    uniprot_id = models.CharField(max_length=200)
    
    def __str__(self):
    	return self.name


class ChemokineType(models.Model):
    chemokinetype = models.CharField(max_length=200)
    
    def __str__(self):
    	return self.name
 
 
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
        db_table = 'gene'
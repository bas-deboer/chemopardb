from django.db import models

class Model(models.Model):
    name = models.CharField(max_length=200, null=True)
    pdb_id = models.CharField(max_length=4, null=True)
    model_id = models.CharField(max_length=4, null=True)
    entity_ref_id = models.CharField(max_length=10, null=True)
    db_name = models.CharField(max_length=100, null=True)
    db_code = models.CharField(max_length=100, null=True)
    sequence = models.TextField(null=True)
    mmcif_text = models.TextField(null=True)
    chemokine_similarity_score = models.CharField(null=True, max_length=100)
    
    def __str__(self):
    	return self.name
 
    class Meta():
        db_table = 'model'

class Partner(models.Model):
    name = models.CharField(max_length=50, null=True)
    pdb_id = models.CharField(max_length=4, null=True)

class Model_PDB(models.Model):
    pdbdata = models.ForeignKey('PDBData', on_delete=models.CASCADE)
    name = models.SlugField(max_length=100, unique=True, null=True)
    #sequence = models.ForeignKey('protein.Sequence', on_delete=models.CASCADE)

    def __str__(self):
        return '{} {}'.format(self.domain, self.chain)

    class Meta():
        db_table = 'Model_PDB'    


class PDBData(models.Model):
    pdb = models.TextField()

    def __str__(self):
        return ('pdb data object')

    class Meta():
        db_table = 'Model_pdbdata'
from django.db import models

class Partner(models.Model):
    name = models.CharField(max_length=255, unique=True)
    description = models.TextField(null=True, blank=True)
    type = models.CharField(max_length=100, null=True, blank=True)
    structures = models.ManyToManyField('structure.Structure', related_name="partners", blank=True)

    def __str__(self):
        return self.name

    class Meta:
        db_table = 'partner'
        verbose_name = 'Partner'
        verbose_name_plural = 'Partners'


class PartnerEntity(models.Model):
    partner = models.ForeignKey(Partner, on_delete=models.CASCADE, related_name='entities')
    entity = models.ForeignKey('structure.Entity', on_delete=models.CASCADE, related_name='partners')
    unp_accession = models.CharField(max_length=50, null=True, blank=True)
    pfam_accession = models.CharField(max_length=50, null=True, blank=True)
    comp_id = models.CharField(max_length=50, null=True, blank=True)
    chembl_id = models.CharField(max_length=255, null=True, blank=True)

    def __str__(self):
        return f"{self.partner.name} - {self.entity.name}"

    class Meta:
        db_table = 'partner_entity'

        
class Partner_PDB(models.Model):
    pdbdata = models.ForeignKey('PDBData', on_delete=models.CASCADE)
    name = models.SlugField(max_length=100, unique=True, null=True)

    def __str__(self):
        return '{} {}'.format(self.domain, self.chain)

    class Meta():
        db_table = 'partner_pdb'    


class PDBData(models.Model):
    pdb = models.TextField()

    def __str__(self):
        return ('pdb data object')

    class Meta():
        db_table = 'partner_pdbdata'
from django.db import models


class Partner(models.Model):
    name = models.CharField(max_length=255, unique=True, null=True)
    partner_type = models.CharField(max_length=255, null=True, blank=True)
    
    # structure definition
    smiles = models.TextField(null=True)
    inchikey = models.CharField(max_length=27, null=True, unique=True)
    clean_inchikey = models.CharField(max_length=27, null=True)
    sequence = models.CharField(max_length=1000, null=True)

    # Ligand properties
    mw = models.DecimalField(max_digits=15, decimal_places=3, null=True)
    rotatable_bonds = models.SmallIntegerField(null=True)
    hacc = models.SmallIntegerField(null=True)
    hdon = models.SmallIntegerField(null=True)
    logp = models.DecimalField(max_digits=10, decimal_places=3, null=True)    

    def __str__(self):
        if self.name:
            return self.name

    class Meta():
        db_table = 'partner'
        
class PartnerType(models.Model):
    slug = models.SlugField(max_length=20, unique=True)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'ligand_type'
        
        
class PartnerProteinStructure(models.Model):
    structure = models.ForeignKey(
        'structure.Structure', on_delete=models.CASCADE, null=True)
    partner = models.ForeignKey('partner.Partner', on_delete=models.CASCADE)
    chain = models.CharField(max_length=20)
    #model = models.ForeignKey(
        #'structure.StructureModel', on_delete=models.CASCADE, null=True)

    def __str__(self):
        return '<ProteinPartner: {} {} {}>'.format(self.structure, self.ligand, self.chain)

    class Meta():
        db_table = "partner_protein_structure"    
        
        
class Partner_PDB(models.Model):
    pdbdata = models.ForeignKey('PDBData', on_delete=models.CASCADE)
    name = models.SlugField(max_length=100, unique=True, null=True)
    #sequence = models.ForeignKey('protein.Sequence', on_delete=models.CASCADE)

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
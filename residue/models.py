from django.db import models


class Residue(models.Model):
    protein = models.ForeignKey('protein.Protein', null=True, on_delete=models.CASCADE)
    #protein_segment = models.ForeignKey('protein.ProteinSegment', null=True, on_delete=models.CASCADE)
    sequence_number = models.SmallIntegerField(null=True)
    ccn_number = models.CharField(max_length=10, null=True, blank=True)
    amino_acid = models.CharField(max_length=1)
    amino_acid_three_letter = models.CharField(max_length=3, null=True)

    def __str__(self):
        return self.amino_acid + str(self.sequence_number)
    
    @property
    def segment(self):
        if self.ccn_number:
            ccn = self.ccn_number
            if ccn.startswith('NTc.Cm'):
                return 'N-term'
            if ccn.startswith('CX.'):
                return 'CX'
            if ccn.startswith('cxb1.'):
                return 'N-loop'
            if ccn.startswith('B1.'):
                return 'B1'
            if ccn.startswith('b1b2.'):
                return '30s-loop'
            if ccn.startswith('B2.'):
                return 'B2'
            if ccn.startswith('b2b3.'):
                return '40s-loop'
            if ccn.startswith('B3.'):
                return 'B3'
            if ccn.startswith('b3h.'):
                return '50s-loop'
            if ccn.startswith('H.'):
                return 'Helix'
            if ccn.startswith('CT.'):
                return 'C-term'
        return None

    class Meta():
        db_table = 'residue'
        ordering = ['sequence_number']


class ResidueGenericNumber(models.Model):
    scheme = models.ForeignKey('ResidueNumberingScheme', on_delete=models.CASCADE)
    protein_segment = models.ForeignKey('protein.ProteinSegment', related_name='generic_numbers', null=True, on_delete=models.CASCADE)
    label = models.CharField(db_index=True, max_length=12)

    def __str__(self):
        return self.label

    class Meta():
        db_table = 'residue_generic_number'
        unique_together = ('scheme', 'label')


class ResidueNumberingScheme(models.Model):
    parent = models.ForeignKey('self', null=True, on_delete=models.CASCADE)
    slug = models.SlugField(max_length=20)
    short_name = models.CharField(max_length=20)
    name = models.CharField(max_length=100)

    def __str__(self):
        return self.name

    class Meta():
        db_table = 'residue_generic_numbering_scheme'
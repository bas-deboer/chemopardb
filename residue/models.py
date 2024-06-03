from django.db import models


class Residue(models.Model):
    protein = models.ForeignKey('protein.Protein', null=True, on_delete=models.CASCADE)
    #protein_segment = models.ForeignKey('protein.ProteinSegment', null=True, on_delete=models.CASCADE)
    sequence_number = models.SmallIntegerField(null=True)
    generic_number = models.SmallIntegerField(null=True)
    amino_acid = models.CharField(max_length=1)
    amino_acid_three_letter = models.CharField(max_length=3, null=True)

    def __str__(self):
        return self.amino_acid + str(self.sequence_number)

    @property
    def segment(self):
        if self.generic_number is not None:
            if 1 <= self.generic_number <= 62:
                return 'N-term'
            if 62 <= self.generic_number <= 67:
                return 'CX'
            if 67 <= self.generic_number <= 82:
                return 'N-loop'
            if 82 <= self.generic_number <= 89:
                return 'B1'
            if 89 <= self.generic_number <= 100:
                return '30s-loop'
            if 100 <= self.generic_number <= 105:
                return 'B2'
            if 105 <= self.generic_number <= 114:
                return '40s-loop'
            if 114 <= self.generic_number <= 118:
                return 'B3'
            if 118 <= self.generic_number <= 122:
                return '50s-loop'
            if 122 <= self.generic_number <= 135:
                return 'Helix'
            elif 135 <= self.generic_number <= 500:
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
import django.db.models.deletion
from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('partner', '0005_partner_uniprot_id'),
        ('structure', '0018_statstext_predictedstructure_predictedcomplex'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='predictedcomplex',
            name='gpcr_protein',
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='gpcr_partner',
            field=models.ForeignKey(blank=True, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='predicted_complex_models', to='partner.partner'),
        ),
    ]

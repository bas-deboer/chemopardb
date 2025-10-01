from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0014_predictedchemokinepartnerinteraction'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictedchemokinepartnerinteraction',
            name='partner_generic_number',
            field=models.CharField(blank=True, max_length=16, null=True),
        ),
        migrations.AddIndex(
            model_name='predictedchemokinepartnerinteraction',
            index=models.Index(
                fields=['predicted_complex', 'partner_generic_number'],
                name='interaction_predicted_generic_idx',
            ),
        ),
    ]

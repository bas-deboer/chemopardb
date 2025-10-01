from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('residue', '0001_initial'),
        ('interaction', '0015_predictedchemokinepartnerinteraction_generic'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictedchemokinepartnerinteraction',
            name='chemokine_residue',
            field=models.ForeignKey(blank=True, null=True, on_delete=models.SET_NULL, related_name='predicted_complex_interactions', to='residue.residue'),
        ),
    ]

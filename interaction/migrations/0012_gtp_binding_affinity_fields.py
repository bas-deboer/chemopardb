from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interaction', '0011_chemokinegpcrgtpbinding'),
    ]

    operations = [
        migrations.AddField(
            model_name='chemokinegpcrgtpbinding',
            name='affinity',
            field=models.CharField(blank=True, max_length=128, null=True),
        ),
        migrations.AddField(
            model_name='chemokinegpcrgtpbinding',
            name='affinity_parameter',
            field=models.CharField(blank=True, max_length=32, null=True),
        ),
        migrations.AddField(
            model_name='chemokinegpcrgtpbinding',
            name='concentration_range',
            field=models.CharField(blank=True, max_length=64, null=True),
        ),
        migrations.AddField(
            model_name='chemokinegpcrgtpbinding',
            name='original_affinity',
            field=models.CharField(blank=True, max_length=128, null=True),
        ),
        migrations.AddField(
            model_name='chemokinegpcrgtpbinding',
            name='original_affinity_relation',
            field=models.CharField(blank=True, max_length=8, null=True),
        ),
        migrations.AddField(
            model_name='chemokinegpcrgtpbinding',
            name='original_affinity_type',
            field=models.CharField(blank=True, max_length=64, null=True),
        ),
    ]


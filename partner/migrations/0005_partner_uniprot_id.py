from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('partner', '0004_partner_gpcrdb_url'),
    ]

    operations = [
        migrations.AddField(
            model_name='partner',
            name='uniprot_id',
            field=models.CharField(max_length=20, null=True, blank=True),
        ),
    ]


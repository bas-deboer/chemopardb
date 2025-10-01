from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('partner', '0003_alter_partner_options_partner_structures_and_more'),
    ]

    operations = [
        migrations.AddField(
            model_name='partner',
            name='gpcrdb_url',
            field=models.URLField(blank=True, null=True),
        ),
    ]


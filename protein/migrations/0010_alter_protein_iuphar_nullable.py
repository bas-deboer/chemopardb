from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0009_alter_protein_web_links'),
    ]

    operations = [
        migrations.AlterField(
            model_name='protein',
            name='iuphar',
            field=models.CharField(max_length=50, null=True, blank=True),
        ),
    ]


# Generated by Django 5.1.2 on 2025-04-03 08:37

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('protein', '0007_alter_protein_web_links'),
    ]

    operations = [
        migrations.DeleteModel(
            name='ProteinState',
        ),
        migrations.AlterModelTable(
            name='gene',
            table='protein_gene',
        ),
        migrations.DeleteModel(
            name='Target',
        ),
    ]

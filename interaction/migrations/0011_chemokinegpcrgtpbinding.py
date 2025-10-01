from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('partner', '0005_partner_uniprot_id'),
        ('protein', '0010_alter_protein_iuphar_nullable'),
        ('interaction', '0010_remove_residuefragmentinteraction_interaction_type_and_more'),
    ]

    operations = [
        migrations.CreateModel(
            name='ChemokineGpcrGtPBinding',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('chemokine_gtp_id', models.CharField(blank=True, max_length=32, null=True)),
                ('gpcr_gtp_target_id', models.CharField(blank=True, max_length=32, null=True)),
                ('action', models.CharField(blank=True, help_text='e.g., agonist/antagonist/modulator', max_length=64, null=True)),
                ('notes', models.TextField(blank=True, null=True)),
                ('source_url', models.URLField(blank=True, null=True)),
                ('chemokine', models.ForeignKey(on_delete=models.deletion.CASCADE, related_name='gtp_gpcr_bindings', to='protein.protein')),
                ('gpcr_partner', models.ForeignKey(on_delete=models.deletion.CASCADE, related_name='gtp_chemokine_bindings', to='partner.partner')),
            ],
            options={
                'verbose_name': 'Chemokine–GPCR GtP Binding',
                'verbose_name_plural': 'Chemokine–GPCR GtP Bindings',
                'db_table': 'interaction_chemokine_gpcr_gtp_binding',
                'unique_together': {('chemokine', 'gpcr_partner')},
            },
        ),
    ]


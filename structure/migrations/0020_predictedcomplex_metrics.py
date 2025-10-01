import re

from django.db import migrations, models


METRIC_LABELS = [
    ('ipTM', 'iptm'),
    ('pTM', 'ptm'),
    ('ranking confidence', 'ranking_confidence'),
    ('mean pLDDT', 'mean_plddt'),
    ('mean pLDDT (chain A)', 'mean_plddt_A'),
    ('mean pLDDT (chain B)', 'mean_plddt_B'),
    ('mean PAE', 'mean_pae'),
    ('mean interface PAE', 'mean_ipae'),
]

METRIC_PATTERNS = {
    field: re.compile(rf"{re.escape(label)}:\s*([0-9]+(?:\.[0-9]+)?)", re.IGNORECASE)
    for label, field in METRIC_LABELS
}

METRIC_PREFIXES = tuple(f"- {label.lower()}" for label, _ in METRIC_LABELS)


def _migrate_predicted_complex_metrics(apps, schema_editor):
    PredictedComplex = apps.get_model('structure', 'PredictedComplex')

    for complex_obj in PredictedComplex.objects.select_related('stats_text').all():
        stats = getattr(complex_obj, 'stats_text', None)
        raw_text = getattr(stats, 'stats_text', '') if stats else ''

        updates = {}
        if raw_text:
            for field_name, pattern in METRIC_PATTERNS.items():
                match = pattern.search(raw_text)
                if not match:
                    continue
                try:
                    updates[field_name] = float(match.group(1))
                except (TypeError, ValueError):
                    continue

        if updates:
            for field_name, value in updates.items():
                setattr(complex_obj, field_name, value)
            complex_obj.save(update_fields=list(updates.keys()))

        if stats and raw_text:
            lines = raw_text.splitlines()
            cleaned_lines = []
            skip_metrics = False

            for line in lines:
                stripped = line.strip()
                if stripped == 'Metrics:':
                    skip_metrics = True
                    continue
                if skip_metrics:
                    lower = stripped.lower()
                    if lower.startswith('- '):
                        if any(lower.startswith(prefix) for prefix in METRIC_PREFIXES):
                            continue
                    skip_metrics = False
                cleaned_lines.append(line)

            cleaned_text = '\n'.join(cleaned_lines)
            if cleaned_text and not cleaned_text.endswith('\n'):
                cleaned_text += '\n'

            if cleaned_text != raw_text:
                stats.stats_text = cleaned_text
                stats.save(update_fields=['stats_text'])


class Migration(migrations.Migration):

    dependencies = [
        ('structure', '0019_predictedcomplex_gpcr_partner'),
    ]

    operations = [
        migrations.AddField(
            model_name='predictedcomplex',
            name='iptm',
            field=models.FloatField(blank=True, help_text='AlphaFold ipTM score', null=True),
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='ptm',
            field=models.FloatField(blank=True, help_text='AlphaFold pTM score', null=True),
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='ranking_confidence',
            field=models.FloatField(blank=True, help_text='AlphaFold ranking confidence', null=True),
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='mean_plddt',
            field=models.FloatField(blank=True, help_text='Mean predicted LDDT across all chains', null=True),
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='mean_plddt_A',
            field=models.FloatField(blank=True, help_text='Mean predicted LDDT for chain A', null=True),
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='mean_plddt_B',
            field=models.FloatField(blank=True, help_text='Mean predicted LDDT for chain B', null=True),
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='mean_pae',
            field=models.FloatField(blank=True, help_text='Mean predicted aligned error', null=True),
        ),
        migrations.AddField(
            model_name='predictedcomplex',
            name='mean_ipae',
            field=models.FloatField(blank=True, help_text='Mean interface predicted aligned error', null=True),
        ),
        migrations.RunPython(_migrate_predicted_complex_metrics, migrations.RunPython.noop),
    ]

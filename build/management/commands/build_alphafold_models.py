import logging
import os
from datetime import date

import requests
from django.core.management.base import BaseCommand
from django.db import transaction

from Bio.PDB import PDBParser

from protein.models import Protein
from structure.models import PredictedStructure, PdbData, StatsText


ALPHAFOLD_URL_TEMPLATES = [
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v4.pdb",
    "https://alphafold.ebi.ac.uk/files/AF-{acc}-F1-model_v3.pdb",
]


class Command(BaseCommand):
    help = "Download AlphaFold models for human chemokines and store as PredictedStructure entries."

    def add_arguments(self, parser):
        parser.add_argument(
            "--accessions",
            nargs="*",
            default=None,
            help="Optional UniProt accessions to process (defaults to all human chemokines)",
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="Overwrite existing PredictedStructure entries for the protein",
        )
        parser.add_argument(
            "--limit",
            type=int,
            default=None,
            help="Process at most N proteins",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Do not write to the database; just report what would happen",
        )

    def handle(self, *args, **options):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)

        accessions = options.get("accessions")
        force = options.get("force")
        limit = options.get("limit")
        dry_run = options.get("dry_run")

        qs = self._select_target_proteins(accessions)
        if limit is not None:
            qs = qs[:limit]

        count = 0
        for protein in qs:
            try:
                created = self._process_protein(protein, force=force, dry_run=dry_run)
                if created:
                    count += 1
            except Exception as e:
                self.logger.error("Failed processing %s (%s): %s", protein.name, protein.uniprot_id, e)

        msg = f"Completed. New PredictedStructure entries: {count}"
        if dry_run:
            msg += " (dry-run)"
        self.stdout.write(self.style.SUCCESS(msg))

    def _select_target_proteins(self, accessions):
        if accessions:
            return Protein.objects.filter(uniprot_id__in=accessions)

        # Heuristic for chemokines based on our import pipeline:
        # - subfamily (ProteinFamily name) is one of common chemokine families
        # - species uses UniProt 5-letter code; 'HUMAN' for Homo sapiens
        families = ["CCL", "CXCL", "XCL", "CX3CL"]
        return Protein.objects.filter(species="Homo sapiens", subfamily__in=families).order_by("name")

    def _process_protein(self, protein: Protein, force: bool, dry_run: bool) -> bool:
        acc = (protein.uniprot_id or '').strip()
        if not acc:
            self.logger.warning("Skipping %s: missing UniProt accession", protein)
            return False

        existing = PredictedStructure.objects.filter(protein=protein, method__icontains="AlphaFold")
        if existing.exists() and not force:
            self.logger.info("Skipping %s (%s): AlphaFold model already present", protein.name, acc)
            return False

        pdb_text, url_used = self._fetch_alphafold_pdb(acc)
        if not pdb_text:
            self.logger.warning("No AlphaFold PDB for %s (%s)", protein.name, acc)
            return False

        # Remove signal peptide residues if annotated
        pdb_text = self._strip_signal_sequence(pdb_text, protein)

        stats_str = self._compute_plddt_stats(pdb_text)
        method = f"AlphaFold DB ({os.path.basename(url_used)})"

        if dry_run:
            self.logger.info("[dry-run] Would create PredictedStructure for %s (%s)", protein.name, acc)
            self.logger.info("[dry-run] Stats: %s", stats_str)
            return True

        with transaction.atomic():
            pdb_obj = PdbData.objects.create(pdb=pdb_text)
            stats_obj = StatsText.objects.create(stats_text=stats_str)

            if existing.exists():
                existing.delete()

            PredictedStructure.objects.create(
                protein=protein,
                state=None,
                pdb_data=pdb_obj,
                stats_text=stats_obj,
                date_generated=date.today(),
                method=method,
            )

        self.logger.info("Created AlphaFold PredictedStructure for %s (%s)", protein.name, acc)
        return True

    def _fetch_alphafold_pdb(self, accession: str):
        for tmpl in ALPHAFOLD_URL_TEMPLATES:
            url = tmpl.format(acc=accession)
            try:
                r = requests.get(url, timeout=30)
            except Exception:
                continue
            if r.status_code == 200 and r.text and r.text.startswith("HEADER"):
                return r.text, url
        return None, None

    def _compute_plddt_stats(self, pdb_text: str) -> str:
        # AlphaFold encodes pLDDT in B-factor field; compute simple summary
        plddts = []
        for line in pdb_text.splitlines():
            if not line.startswith("ATOM"):
                continue
            try:
                # Columns 61-66 store B-factor in PDB format (1-based), Python slicing [60:66]
                b = float(line[60:66].strip())
                plddts.append(b)
            except Exception:
                continue

        if not plddts:
            return "No pLDDT values parsed from PDB."

        n = len(plddts)
        mean_plddt = sum(plddts) / n
        high = sum(1 for v in plddts if v >= 90)
        conf = sum(1 for v in plddts if v >= 70)
        low = sum(1 for v in plddts if v < 50)

        return (
            f"AlphaFold pLDDT summary: mean={mean_plddt:.1f}, N={n}; "
            f">=90: {high}, >=70: {conf}, <50: {low}"
        )

    def _strip_signal_sequence(self, pdb_text: str, protein: Protein) -> str:
        """
        Remove residues within the signal peptide [start, end] from the PDB text.
        Keeps all non-ATOM/HETATM records as-is. Works for typical AFDB single-chain PDBs.
        """
        try:
            sig = getattr(protein, 'signal_sequence', None)
        except Exception:
            sig = None

        if not sig or not sig.start or not sig.end:
            return pdb_text

        start_idx = int(sig.start)
        end_idx = int(sig.end)
        if end_idx < start_idx:
            return pdb_text

        kept_lines = []
        for line in pdb_text.splitlines(keepends=True):
            rec = line[:6]
            if rec.startswith('ATOM') or rec.startswith('HETATM'):
                # PDB resSeq is columns 23-26 (1-based). Python slice [22:26]
                try:
                    resseq_str = line[22:26].strip()
                    resseq = int(resseq_str) if resseq_str else None
                except Exception:
                    resseq = None

                if resseq is not None and start_idx <= resseq <= end_idx:
                    # Skip atoms in the signal peptide range
                    continue
            kept_lines.append(line)

        result = ''.join(kept_lines)
        # Ensure trailing newline for PDB compliance
        if not result.endswith('\n'):
            result += '\n'
        return result

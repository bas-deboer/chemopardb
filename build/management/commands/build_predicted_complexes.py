import csv
import logging
import re
from datetime import date
from pathlib import Path
from types import SimpleNamespace
from typing import Dict, Iterable, List, Optional, Tuple, Union

from django.conf import settings
from django.core.management import CommandError
from django.core.management.base import BaseCommand
from django.db import transaction
from django.db.models import Q

from protein.models import GPCRProtein, Protein, ProteinAlias
from residue.models import Residue
from partner.models import Partner
from structure.models import PredictedComplex, PdbData, StatsText

from interaction.calculation2024 import prolif_calculation
from interaction.models import PredictedChemokinePartnerInteraction

InteractionRecord = Dict[str, Optional[Union[str, int]]]


DEFAULT_METHOD = "AlphaFold-Multimer"


class Command(BaseCommand):
    """Load local AlphaFold-Multimer chemokine-GPCR complexes into PredictedComplex."""

    def __init__(self):
        super().__init__()
        self.logger = logging.getLogger(__name__)
        base_dir = Path(settings.BASE_DIR)
        self.data_root = base_dir / settings.DATA_DIR / "structure_data" / "predicted_structures"
        self.default_csv = self.data_root / "panel_sample.csv"
        self.default_models_dir = self.data_root / "models"
        self.interactions_dir = self.data_root / "interactions"

    def add_arguments(self, parser):
        parser.add_argument(
            "--csv",
            default=str(self.default_csv),
            help="Path to the CSV file with complex metadata (default: %(default)s)",
        )
        parser.add_argument(
            "--models-dir",
            default=str(self.default_models_dir),
            help="Directory containing PDB models (default: %(default)s)",
        )
        parser.add_argument(
            "--only",
            nargs="*",
            help="Optional list of complex identifiers to import (matches the 'complex' column).",
        )
        parser.add_argument(
            "--limit",
            type=int,
            help="Process at most N complexes after filtering.",
        )
        parser.add_argument(
            "--force",
            action="store_true",
            help="Overwrite existing PredictedComplex entries for the same pair and method.",
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
            help="Validate the import without writing to the database.",
        )
        parser.add_argument(
            "--method",
            default=DEFAULT_METHOD,
            help="Method label stored on PredictedComplex.method (default: %(default)s)",
        )

    def handle(self, *args, **options):
        self.logger.setLevel(logging.INFO)

        csv_path = self._resolve_path(options["csv"])
        models_dir = self._resolve_path(options["models_dir"])
        limit = options.get("limit")
        force = options.get("force")
        dry_run = options.get("dry_run")
        method = (options.get("method") or DEFAULT_METHOD).strip()
        only_tokens = {
            token.strip(): token.strip()
            for token in (options.get("only") or [])
            if token and token.strip()
        }

        if not csv_path.exists():
            raise CommandError(f"Metadata CSV not found at {csv_path}")
        if not models_dir.exists() or not models_dir.is_dir():
            raise CommandError(f"Models directory not found: {models_dir}")

        if not dry_run:
            self.interactions_dir.mkdir(parents=True, exist_ok=True)

        created = 0
        updated = 0
        skipped = 0
        selected = 0

        for row in self._iter_rows(csv_path):
            complex_id = (row.get("complex") or "").strip()
            if not complex_id:
                self.logger.warning("Skipping row without 'complex' identifier: %s", row)
                skipped += 1
                continue

            if only_tokens and complex_id not in only_tokens:
                continue

            selected += 1
            if limit and selected > limit:
                break

            result = self._process_row(row, models_dir, method, dry_run, force)
            if result == "created":
                created += 1
            elif result == "updated":
                updated += 1
            else:
                skipped += 1

        summary = (
            f"Completed import: created={created}, updated={updated}, skipped={skipped}."
        )
        if dry_run:
            summary += " (dry-run)"
        self.stdout.write(self.style.SUCCESS(summary))




    def _iter_rows(self, csv_path: Path) -> Iterable[Dict[str, str]]:
        with csv_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            for row in reader:
                yield row

    def _process_row(
        self,
        row: Dict[str, str],
        models_dir: Path,
        method: str,
        dry_run: bool,
        force: bool,
    ) -> str:
        complex_id = (row.get("complex") or "").strip()
        gpcr_symbol = (row.get("id1") or "").strip()
        chemokine_symbol = (row.get("id2") or "").strip()

        if not gpcr_symbol or not chemokine_symbol:
            self.logger.warning("Skipping %s: missing id1/id2 values", complex_id)
            return "skipped"

        gpcr_partner = self._resolve_partner(gpcr_symbol)
        if gpcr_partner is None:
            self.logger.warning("Skipping %s: GPCR '%s' not found in Partner table", complex_id, gpcr_symbol)
            return "skipped"

        chemokine = self._resolve_protein(chemokine_symbol, role="chemokine")
        if chemokine is None:
            self.logger.warning(
                "Skipping %s: chemokine '%s' not found in Protein table",
                complex_id,
                chemokine_symbol,
            )
            return "skipped"

        pdb_filename = complex_id if complex_id.lower().endswith(".pdb") else f"{complex_id}.pdb"
        pdb_path = models_dir / pdb_filename
        if not pdb_path.exists():
            self.logger.warning("Skipping %s: PDB file not found (%s)", complex_id, pdb_path)
            return "skipped"

        original_pdb_text = self._read_pdb_text(pdb_path)

        extra_stats: List[str] = []
        interaction_records: List[InteractionRecord] = []
        pdb_text_to_store = original_pdb_text

        if dry_run:
            extra_stats.append(" - Interaction fingerprint not generated (dry-run).")
        else:
            _, interaction_lines, interaction_records = self._calculate_interactions(
                complex_id=complex_id,
                pdb_text=pdb_text_to_store,
                row=row,
                force=force,
            )
            extra_stats.extend(interaction_lines)

        metrics = self._extract_metrics(row)
        stats_text = self._build_stats_text(row, extra_stats)

        existing = PredictedComplex.objects.filter(
            chemokine_protein=chemokine,
            gpcr_partner=gpcr_partner,
            method=method,
        ).first()

        if existing and not force:
            self.logger.info(
                "Skipping %s: PredictedComplex already exists (id=%s). Use --force to overwrite.",
                complex_id,
                existing.id,
            )
            return "skipped"

        if dry_run:
            action = "update" if existing else "create"
            self.logger.info(
                "[dry-run] Would %s PredictedComplex %s (GPCR=%s, chemokine=%s)",
                action,
                complex_id,
                gpcr_symbol,
                chemokine_symbol,
            )
            return "updated" if existing else "created"

        with transaction.atomic():
            if existing:
                self._update_existing(existing, pdb_text_to_store, stats_text, method, metrics, gpcr_partner)
                self.logger.info(
                    "Updated PredictedComplex %s (id=%s)",
                    complex_id,
                    existing.id,
                )
                if not dry_run:
                    self._persist_predicted_interactions(existing, interaction_records)
                return "updated"

            new_obj = self._create_complex(chemokine, gpcr_partner, pdb_text_to_store, stats_text, method, metrics)
            self.logger.info(
                "Created PredictedComplex %s (id=%s)",
                complex_id,
                new_obj.id,
            )
            if not dry_run:
                self._persist_predicted_interactions(new_obj, interaction_records)
            return "created"

    def _resolve_path(self, value: str) -> Path:
        path = Path(value)
        if not path.is_absolute():
            path = Path(settings.BASE_DIR) / path
        return path

    def _resolve_protein(self, symbol: str, role: str) -> Optional[Protein]:
        token = (symbol or "").strip()
        if not token:
            return None

        matches = list(
            Protein.objects.filter(
                Q(gene_name__iexact=token)
                | Q(name__iexact=token)
                | Q(uniprot_id__iexact=token)
            )
        )

        if not matches:
            prefixed = f"{token}_"
            matches = list(Protein.objects.filter(name__istartswith=prefixed))

        if not matches:
            alias = (
                ProteinAlias.objects.filter(name__iexact=token)
                .select_related("protein")
                .first()
            )
            if alias:
                matches = [alias.protein]

        if matches:
            preferred = [p for p in matches if (p.species or "").lower() == "homo sapiens"]
            chosen = preferred[0] if preferred else matches[0]

            if len(matches) > 1:
                self.logger.debug(
                    "Multiple matches for %s '%s'; using %s (id=%s)",
                    role,
                    token,
                    chosen.name,
                    chosen.id,
                )

            return chosen

        return None

    def _resolve_partner(self, symbol: str) -> Optional[Partner]:
        token = (symbol or "").strip()
        if not token:
            return None

        partner = Partner.objects.filter(
            Q(name__iexact=token) | Q(uniprot_id__iexact=token)
        ).first()
        if partner:
            return partner

        simplified = token.replace('-', '').replace(' ', '')
        if simplified != token:
            partner = Partner.objects.filter(name__iexact=simplified).first()
            if partner:
                return partner

        return None


    def _read_pdb_text(self, pdb_path: Path) -> str:
        try:
            text = pdb_path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = pdb_path.read_text(encoding="latin-1")
        if not text.endswith("\n"):
            text += "\n"
        return text

    def _build_stats_text(
        self,
        row: Dict[str, str],
        extra_lines: Optional[Iterable[str]] = None,
    ) -> str:
        complex_id = (row.get("complex") or "").strip()
        gpcr_symbol = (row.get("id1") or "").strip()
        chemokine_symbol = (row.get("id2") or "").strip()

        lines = []
        header = f"Complex: {complex_id}" if complex_id else "Complex metadata"
        lines.append(header)
        if gpcr_symbol or chemokine_symbol:
            lines.append(
                f"Pair: GPCR={gpcr_symbol or '?'}; chemokine={chemokine_symbol or '?'}"
            )

        for chain_id in ("A", "B"):
            seq = (row.get(chain_id) or "").strip()
            if seq:
                lines.append(f"Chain {chain_id} length: {len(seq)} aa")

        notes: List[str] = []
        if extra_lines:
            notes = [line for line in extra_lines if line]
            if notes:
                lines.append("")
                lines.append("Processing notes:")
                lines.extend(notes)

        stats_text = "\n".join(lines)
        if not stats_text.endswith("\n"):
            stats_text += "\n"
        return stats_text

    def _extract_metrics(
        self,
        row: Dict[str, str],
    ) -> Dict[str, Optional[float]]:
        metric_keys = [
            'iptm',
            'ptm',
            'ranking_confidence',
            'mean_plddt',
            'mean_plddt_A',
            'mean_plddt_B',
            'mean_pae',
            'mean_ipae',
        ]
        metrics: Dict[str, Optional[float]] = {}
        for key in metric_keys:
            metrics[key] = self._safe_float(row.get(key))
        return metrics

    def _calculate_interactions(
        self,
        complex_id: str,
        pdb_text: str,
        row: Dict[str, str],
        force: bool,
    ) -> Tuple[Optional[Path], List[str], List[InteractionRecord]]:
        lines: List[str] = []
        records: List[InteractionRecord] = []

        if not pdb_text.strip():
            lines.append(" - Interaction fingerprint not generated (empty PDB text).")
            return None, lines, records

        target_path = self.interactions_dir / f"{complex_id}_ifp.pkl"

        try:
            self.interactions_dir.mkdir(parents=True, exist_ok=True)

            chemokine_chain = (row.get("chemokine_chain") or "B").strip() or "B"
            partner_chain = (row.get("gpcr_chain") or "A").strip() or "A"

            chemokine_residues = self._parse_residue_tokens(row.get("chemokine_residues"))
            partner_residues = self._parse_residue_tokens(row.get("partner_residues"))

            if chemokine_residues is None:
                chem_seq = (row.get("chemokine_sequence") or row.get("B") or "").strip()
                if chem_seq:
                    chemokine_residues = list(range(1, len(chem_seq) + 1))

            if partner_residues is None:
                partner_seq = (row.get("gpcr_sequence") or row.get("A") or "").strip()
                if partner_seq:
                    partner_residues = list(range(1, len(partner_seq) + 1))

            temp_complex = SimpleNamespace(pdb_data=SimpleNamespace(pdb=pdb_text))

            fp = prolif_calculation(
                chemokine_partner_complex=temp_complex,
                chemokine_chain=chemokine_chain,
                chemokine_residues=chemokine_residues,
                partner_chains=[partner_chain],
                partner_residues=partner_residues,
            )

            fp.to_pickle(str(target_path))
            lines.append(f" - Interaction fingerprint saved to {target_path.name}.")
            lines.extend(self._summarize_interactions(fp.ifp))
            records = self._flatten_interactions(fp.ifp)
            return target_path, lines, records
        except Exception as exc:  # pragma: no cover - depends on ProLIF runtime
            self.logger.error(
                "Failed to calculate interactions for %s: %s",
                complex_id,
                exc,
                exc_info=True,
            )
            lines.append(f" - Interaction calculation failed: {exc}")
            return None, lines, records

    def _parse_residue_tokens(self, raw: Optional[str]) -> Optional[List[int]]:
        if not raw:
            return None

        residues: List[int] = []
        for token in re.split(r"[;,\s]+", raw.strip()):
            if not token:
                continue
            if "-" in token:
                start, end = token.split("-", 1)
                try:
                    start_idx = int(start)
                    end_idx = int(end)
                except ValueError:
                    continue
                if end_idx < start_idx:
                    start_idx, end_idx = end_idx, start_idx
                residues.extend(range(start_idx, end_idx + 1))
            else:
                try:
                    residues.append(int(token))
                except ValueError:
                    continue

        return residues or None

    def _summarize_interactions(self, interaction_map) -> List[str]:
        counts: Dict[str, int] = {}
        contact_pairs = 0

        for residues in interaction_map.values():
            for interactions in residues.values():
                if interactions:
                    contact_pairs += 1
                for interaction_type in interactions.keys():
                    counts[interaction_type] = counts.get(interaction_type, 0) + 1

        lines: List[str] = []
        if contact_pairs:
            lines.append(f" - Interaction residue pairs assessed: {contact_pairs}")

        if counts:
            for interaction_type, count in sorted(counts.items(), key=lambda kv: kv[0].lower()):
                lines.append(f" - {interaction_type}: {count}")
        elif not contact_pairs:
            lines.append(" - No interactions detected by ProLIF.")

        return lines

    def _flatten_interactions(self, interaction_map) -> List[InteractionRecord]:
        rows: List[InteractionRecord] = []
        for residues in interaction_map.values():
            for residue_pair, interactions in residues.items():
                if not interactions:
                    continue

                chem_res = residue_pair[0]
                partner_res = residue_pair[1]
                chem_chain = getattr(chem_res, "chain", None)
                partner_chain = getattr(partner_res, "chain", None)
                chem_number = self._safe_int(getattr(chem_res, "number", None))
                partner_number = self._safe_int(getattr(partner_res, "number", None))
                chem_name = getattr(chem_res, "name", None)
                partner_name = getattr(partner_res, "name", None)

                for interaction_type in interactions.keys():
                    rows.append(
                        {
                            "chemokine_chain": chem_chain,
                            "chemokine_residue_number": chem_number,
                            "chemokine_residue_name": chem_name,
                            "partner_chain": partner_chain,
                            "partner_residue_number": partner_number,
                            "partner_residue_name": partner_name,
                            "partner_generic_number": None,
                            "interaction_type": interaction_type,
                        }
                    )

        return rows

    def _persist_predicted_interactions(
        self,
        predicted_complex: PredictedComplex,
        records: List[InteractionRecord],
    ) -> None:
        if predicted_complex is None:
            return

        PredictedChemokinePartnerInteraction.objects.filter(
            predicted_complex=predicted_complex
        ).delete()

        if not records:
            return

        generic_map = self._partner_generic_map(predicted_complex.gpcr_partner)
        chemokine_map = self._chemokine_residue_map(predicted_complex.chemokine_protein)

        objects = [
            PredictedChemokinePartnerInteraction(
                predicted_complex=predicted_complex,
                chemokine_residue=chemokine_map.get(
                    self._safe_int(record.get("chemokine_residue_number"))
                ),
                chemokine_chain=record.get("chemokine_chain"),
                chemokine_residue_number=record.get("chemokine_residue_number"),
                chemokine_residue_name=record.get("chemokine_residue_name"),
                partner_chain=record.get("partner_chain"),
                partner_residue_number=record.get("partner_residue_number"),
                partner_residue_name=record.get("partner_residue_name"),
                partner_generic_number=self._lookup_generic_number(
                    generic_map,
                    record.get("partner_residue_number"),
                ),
                interaction_type=record.get("interaction_type") or "Unknown",
            )
            for record in records
        ]

        PredictedChemokinePartnerInteraction.objects.bulk_create(objects, ignore_conflicts=True)

    def _lookup_generic_number(
        self,
        generic_map: Dict[int, str],
        sequence_number: Optional[Union[str, int]],
    ) -> Optional[str]:
        if not generic_map:
            return None
        if sequence_number in (None, ""):
            return None
        try:
            seq_int = int(sequence_number)
        except (TypeError, ValueError):
            seq_int = self._safe_int(str(sequence_number))
        if seq_int is None:
            return None
        return generic_map.get(seq_int)

    def _chemokine_residue_map(
        self,
        chemokine: Optional[Protein],
    ) -> Dict[int, Residue]:
        if chemokine is None:
            return {}

        residues = Residue.objects.filter(
            protein=chemokine,
            sequence_number__isnull=False,
        ).only("id", "sequence_number")

        return {res.sequence_number: res for res in residues}

    def _safe_float(self, value: Optional[str]) -> Optional[float]:
        if value in (None, ""):
            return None
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    def _safe_int(self, value: Optional[str]) -> Optional[int]:
        if value in (None, ""):
            return None
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    def _partner_generic_map(self, partner: Optional[Partner]) -> Dict[int, str]:
        if partner is None:
            return {}

        lookup_candidates = []
        partner_uniprot = (partner.uniprot_id or "").strip()
        if partner_uniprot:
            lookup_candidates.append({"unp_accession__iexact": partner_uniprot})
        lookup_candidates.append({"name__iexact": (partner.name or "").strip()})

        gpcr_obj = None
        for filters in lookup_candidates:
            if not filters.get(next(iter(filters))):
                continue
            gpcr_obj = GPCRProtein.objects.filter(**filters).first()
            if gpcr_obj:
                break

        if not gpcr_obj or not gpcr_obj.residues:
            return {}

        mapping: Dict[int, str] = {}
        residues_data = gpcr_obj.residues

        if isinstance(residues_data, dict):
            iterable = residues_data.items()
        else:
            # Unexpected format (e.g., list); convert to iterable of tuples
            iterable = enumerate(residues_data)

        for key, value in iterable:
            try:
                seq_number = int(key)
            except (TypeError, ValueError):
                seq_number = None
                if isinstance(value, dict):
                    seq_candidate = value.get("sequence_number")
                    try:
                        seq_number = int(seq_candidate)
                    except (TypeError, ValueError):
                        seq_number = None
            if seq_number is None:
                continue

            generic = None
            if isinstance(value, dict):
                generic = value.get("generic_number") or value.get("display_generic_number")
            elif isinstance(value, str):
                generic = value

            if generic:
                mapping[seq_number] = generic

        return mapping

    def _create_complex(
        self,
        chemokine: Protein,
        gpcr_partner: Partner,
        pdb_text: str,
        stats_text: str,
        method: str,
        metrics: Dict[str, Optional[float]],
    ) -> PredictedComplex:
        pdb_obj = PdbData.objects.create(pdb=pdb_text)
        stats_obj = StatsText.objects.create(stats_text=stats_text)
        metric_values = {k: v for k, v in (metrics or {}).items() if v is not None}
        return PredictedComplex.objects.create(
            chemokine_protein=chemokine,
            gpcr_partner=gpcr_partner,
            state=None,
            pdb_data=pdb_obj,
            stats_text=stats_obj,
            date_generated=date.today(),
            method=method,
            **metric_values,
        )

    def _update_existing(
        self,
        instance: PredictedComplex,
        pdb_text: str,
        stats_text: str,
        method: str,
        metrics: Dict[str, Optional[float]],
        gpcr_partner: Partner,
    ) -> None:
        if instance.pdb_data:
            instance.pdb_data.pdb = pdb_text
            instance.pdb_data.save(update_fields=['pdb'])
        else:
            instance.pdb_data = PdbData.objects.create(pdb=pdb_text)

        if instance.stats_text:
            instance.stats_text.stats_text = stats_text
            instance.stats_text.save(update_fields=['stats_text'])
        else:
            instance.stats_text = StatsText.objects.create(stats_text=stats_text)

        metric_field_names = []
        for field_name, value in (metrics or {}).items():
            setattr(instance, field_name, value)
            metric_field_names.append(field_name)

        instance.gpcr_partner = gpcr_partner
        instance.method = method
        instance.date_generated = date.today()
        update_fields = ['pdb_data', 'stats_text', 'gpcr_partner', 'method', 'date_generated'] + metric_field_names
        instance.save(update_fields=update_fields)


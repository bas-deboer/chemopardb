import time
import logging
import requests

from django.core.management.base import BaseCommand

from protein.models import Protein
from partner.models import Partner
from interaction.models import ChemokineGpcrGtPBinding


GTP_BASE = "https://www.guidetopharmacology.org/services"


class Command(BaseCommand):
    help = "Fetches chemokine–GPCR endogenous interactions from GtoP and stores them in ChemokineGpcrGtPBinding"

    def add_arguments(self, parser):
        parser.add_argument("--ligand", type=int, help="Only fetch for a specific GtoP ligand id")
        parser.add_argument("--sleep", type=float, default=0.2, help="Sleep seconds between API calls to be polite")
        parser.add_argument("--refresh", action="store_true", help="Delete existing bindings before loading")

    def handle(self, *args, **opts):
        logger = logging.getLogger(__name__)

        if opts.get("refresh"):
            ChemokineGpcrGtPBinding.objects.all().delete()
            self.stdout.write(self.style.WARNING("Cleared existing GtP bindings"))

        # Preload GPCR partners by name (case-insensitive)
        partners = {p.name.lower(): p for p in Partner.objects.filter(type__iexact="GPCR")}
        if not partners:
            self.stdout.write(self.style.WARNING("No GPCR partners found. Run build_proteins_gpcr first."))

        # Pick chemokines by presence of iuphar id
        qs = Protein.objects.exclude(iuphar__isnull=True).exclude(iuphar="")
        if opts.get("ligand"):
            qs = qs.filter(iuphar=str(opts["ligand"]))

        count = 0
        for chem in qs.iterator():
            ligand_id = str(chem.iuphar).strip()
            try:
                interactions = self.fetch_interactions(ligand_id)
            except Exception as e:
                logger.warning(f"Failed fetching interactions for ligand {ligand_id}: {e}")
                continue

            # Filter: human targets and endogenous interactions; target type GPCR is implicit in our partner list
            for it in interactions:
                if not it.get("endogenous", False):
                    continue
                if (it.get("targetSpecies") or "").lower() != "human":
                    continue

                target_name = (it.get("targetName") or "").strip()
                if not target_name:
                    continue
                partner = partners.get(target_name.lower())
                if not partner:
                    # Could not map by name; skip quietly
                    continue

                action = it.get("action")
                type = it.get("type")
                target_id = it.get("targetId")
                try:
                    ChemokineGpcrGtPBinding.objects.update_or_create(
                        chemokine=chem,
                        gpcr_partner=partner,
                        defaults={
                            "chemokine_gtp_id": ligand_id,
                            "gpcr_gtp_target_id": str(target_id) if target_id is not None else None,
                            "type": type,
                            "action": action,
                            "source_url": f"https://www.guidetopharmacology.org/GRAC/LigandDisplayForward?ligandId={ligand_id}",
                            # Affinity fields
                            "affinity": it.get("affinity") or None,
                            "affinity_parameter": it.get("affinityParameter") or None,
                            "original_affinity": it.get("originalAffinity") or None,
                            "original_affinity_type": it.get("originalAffinityType") or None,
                            "original_affinity_relation": it.get("originalAffinityRelation") or None,
                            "concentration_range": it.get("concentrationRange") or None,
                        },
                    )
                    count += 1
                except Exception as e:
                    logger.warning(f"Upsert failed for {chem.name} ↔ {partner.name}: {e}")

            # polite pause
            time.sleep(float(opts.get("sleep") or 0))

        self.stdout.write(self.style.SUCCESS(f"Stored/updated {count} chemokine–GPCR GtP bindings"))

    def fetch_interactions(self, ligand_id: str):
        url = f"{GTP_BASE}/ligands/{ligand_id}/interactions?format=json"
        r = requests.get(url, timeout=30)
        r.raise_for_status()
        return r.json() or []

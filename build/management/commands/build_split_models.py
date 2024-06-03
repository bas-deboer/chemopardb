from django.core.management.base import BaseCommand, CommandError
from build.management.commands.base_build import Command as BaseBuild
from django.conf import settings
from django.db import transaction

from protein.models import Protein
from structure.models import Structure
from model.models import Model

import os
import io
import sys
import pandas as pd
from pathlib import Path

from Bio import PDB
from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBExceptions import PDBConstructionWarning


class Command(BaseBuild):

    def add_arguments(self, parser):
        parser.add_argument("--debug", default=False, action="store_true", help="Debug mode")
        parser.add_argument("--purge", default=False, action="store_true", help="Purge data")

    def handle(self, *args, **options):

        # Split models
        self.split_structure()

    def split_structure(self):
        
        protein	= ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
        
        alternative_aminoacids = ['000', '004', '005', '00B', '00C', '01W', '0A0', '0A1', '0A2', '0A8', '0A9', '0AA', '0AB', '0AC', '0AF', '0AG', '0AH', '0AK', 
                               '0AR', '0AY', '0AZ', '0BN', '0CS', '0FL', '0NC', '0YG', '143', '175', '192', '193', '1AC', '1LU', '1PA', '1PI', '1TQ', '1TY', '1ZN', 
                               '200', '23F', '23S', '26B', '2AD', '2AG', '2AO', '2CO', '2FM', '2HF', '2LU', '2ML', '2MR', '2MT', '2P0', '2PI', '2PP', '2RA', '2TL',
                               '2TY', '2VA', '2XA', '32S', '32T', '3AH', '3AR', '3CF', '3FG', '3GA', '3GL', '3MD', '3MM', '3MY', '3NF', '3O3', '3PA', '3QN', '3TY',
                               '3XH', '4AF', '4BA', '4CF', '4CY', '4DB', '4DP', '4F3', '4FB', '4FW', '4HT', '4IN', '4MM', '4PH', '4U7', '5AB', '5CS', '5HP', '5OH',
                               '5PG', '5ZA', '6CL', '6CW', '6HC', '6HG', '6HN', '6HT', '7HA', '7JA', '7MN', '9DN', '9DS', '9NE', '9NF', '9NR', '9NV', 'A5N', 'A66',
                               'A9D', 'AA3', 'AA4', 'AA6', 'AAR', 'AB7', 'ABA', 'ABU', 'AC5', 'ACA', 'ACB', 'ACE', 'ACL', 'ACY', 'AEA', 'AEI', 'AFA', 'AGD', 'AGM',
                               'AGQ', 'AGT', 'AHB', 'AHO', 'AHP', 'AIB', 'AKL', 'ALC', 'ALM', 'ALN', 'ALO', 'ALS', 'ALY', 'AME', 'AN6', 'APH', 'API', 'APK', 'APM',
                               'APN', 'APO', 'APP', 'AR2', 'AR4', 'ARM', 'ARO', 'ARV', 'AS2', 'AS9', 'ASA', 'ASB', 'ASE', 'ASI', 'ASK', 'ASL', 'ASM', 'ASQ', 'ASX',
                               'AVN', 'AYA', 'AYG', 'AZK', 'AZS', 'AZY', 'B1F', 'B27', 'B2A', 'B2F', 'B2I', 'B2V', 'B3A', 'B3D', 'B3E', 'B3K', 'B3L', 'B3M', 'B3Q',
                               'B3S', 'B3T', 'B3U', 'B3X', 'B3Y', 'BAL', 'BB6', 'BB7', 'BB8', 'BB9', 'BBC', 'BCC', 'BCS', 'BCX', 'BFD', 'BG1', 'BH2', 'BHD', 'BIF',
                               'BIL', 'BIU', 'BLE', 'BLY', 'BMT', 'BNN', 'BNO', 'BOR', 'BPE', 'BSE', 'BTA', 'BTC', 'BTK', 'BTR', 'BUC', 'BUG', 'C12', 'C1T', 'C1X',
                               'C2N', 'C3Y', 'C4R', 'C5C', 'C66', 'C6C', 'C99', 'CAB', 'CAF', 'CAL', 'CAS', 'CAV', 'CAY', 'CBX', 'CCL', 'CCS', 'CCY', 'CDE', 'CDV',
                               'CEA', 'CEG', 'CFY', 'CGA', 'CGN', 'CGU', 'CH2', 'CH3', 'CH6', 'CH7', 'CHF', 'CHG', 'CHP', 'CHS', 'CIR', 'CJO', 'CLB', 'CLD', 'CLE', 
                               'CLG', 'CLH', 'CLV', 'CME', 'CMH', 'CML', 'CMT', 'CPC', 'CPI', 'CPN', 'CR0', 'CR2', 'CR5', 'CR7', 'CR8', 'CRF', 'CRG', 'CRK', 'CRO',
                               'CRQ', 'CRU', 'CRW', 'CRX', 'CS0', 'CS1', 'CS3', 'CS4', 'CSA', 'CSB', 'CSD', 'CSE', 'CSH', 'CSJ', 'CSO', 'CSP', 'CSR', 'CSS', 'CSU',
                               'CSW', 'CSX', 'CSY', 'CSZ', 'CTE', 'CTH', 'CUC', 'CUD', 'CWD', 'CWR', 'CXM', 'CY0', 'CY1', 'CY3', 'CY4', 'CYA', 'CYD', 'CYF', 'CYG',
                               'CYJ', 'CYM', 'CYQ', 'CYR', 'CYW', 'CZ2', 'CZO', 'CZZ', 'D11', 'D3P', 'D4P', 'DA2', 'DAB', 'DAH', 'DAL', 'DAM', 'DAR', 'DAS', 'DBB',
                               'DBS', 'DBU', 'DBY', 'DBZ', 'DC2', 'DCY', 'DDE', 'DDZ', 'DFI', 'DFO', 'DGH', 'DGL', 'DGN', 'DHA', 'DHI', 'DHL', 'DHN', 'DHP', 'DHV',
                               'DIL', 'DIR', 'DIV', 'DLE', 'DLS', 'DLY', 'DM0', 'DMH', 'DMK', 'DMT', 'DNE', 'DNG', 'DNL', 'DNM', 'DNP', 'DNS', 'DO2', 'DOA', 'DOH',
                               'DON', 'DPL', 'DPN', 'DPP', 'DPQ', 'DPR', 'DSE', 'DSG', 'DSN', 'DSP', 'DTH', 'DTR', 'DTY', 'DVA', 'DYG', 'DYS', 'ECX', 'EFC', 'EHP',
                               'ELY', 'ESB', 'ETA', 'EXY', 'EYG', 'F2F', 'FB5', 'FB6', 'FBE', 'FC0', 'FCL', 'FCY', 'FFM', 'FGA', 'FGL', 'FGP', 'FH7', 'FHL', 'FHO',
                               'FLA', 'FLE', 'FMA', 'FME', 'FOE', 'FOR', 'FP9', 'FPR', 'FRD', 'FT6', 'FTR', 'FTY', 'FVA', 'FZN', 'G01', 'GAL', 'GAU', 'GCU', 'GGL',
                               'GHC', 'GHG', 'GHP', 'GHW', 'GL3', 'GLC', 'GLH', 'GLJ', 'GLM', 'GLQ', 'GLX', 'GLZ', 'GMA', 'GME', 'GND', 'GPL', 'GPN', 'GSC', 'GSU',
                               'GT9', 'GTH', 'GVL', 'GYC', 'GYS', 'H5M', 'HAC', 'HAR', 'HBN', 'HCS', 'HFA', 'HGL', 'HHI', 'HIA', 'HIC', 'HIP', 'HIQ', 'HL2', 'HLU',
                               'HMF', 'HMR', 'HNC', 'HOX', 'HPC', 'HPE', 'HPQ', 'HQA', 'HRG', 'HRP', 'HS8', 'HS9', 'HSE', 'HSL', 'HSO', 'HTI', 'HTN', 'HTR', 'HV5',
                               'HVA', 'HY3', 'HYP', 'HZP', 'I2M', 'I58', 'IAM', 'IAR', 'IAS', 'ICY', 'IEL', 'IEY', 'IGL', 'IIC', 'IIL', 'ILG', 'ILX', 'IML', 'IOY',
                               'IT1', 'IYR', 'IYT', 'JJJ', 'JJK', 'JJL', 'K1R', 'KBE', 'KCX', 'KGC', 'KOR', 'KPI', 'KST', 'KYN', 'KYQ', 'L2A', 'L3O', 'LA2', 'LAA',
                               'LAL', 'LBY', 'LCK', 'LCX', 'LDH', 'LE1', 'LED', 'LEF', 'LEH', 'LEI', 'LEM', 'LEN', 'LEP', 'LET', 'LGY', 'LHC', 'LLO', 'LLP', 'LLY',
                               'LLZ', 'LME', 'LMQ', 'LNT', 'LP6', 'LPD', 'LPG', 'LPS', 'LSO', 'LTA', 'LTR', 'LVG', 'LVN', 'LWY', 'LYM', 'LYN', 'LYR', 'LYU', 'LYX',
                               'LYZ', 'M0H', 'M2L', 'M2S', 'M3L', 'M3R', 'MA', 'MAA', 'MAI', 'MAN', 'MBQ', 'MC1', 'MCB', 'MCG', 'MCL', 'MCS', 'MDH', 'MDO', 'ME0', 
                               'MEA', 'MED', 'MEG', 'MEN', 'MEQ',  'MEU', 'MF3', 'MFC', 'MGG', 'MGN', 'MGY', 'MH1', 'MH6', 'MHE', 'MHL', 'MHO', 'MHS', 'MHT', 'MHU',
                               'MHV', 'MHW', 'MIR', 'MIS', 'MK8', 'ML3', 'MLE', 'MLL', 'MLU', 'MLY', 'MLZ', 'MME', 'MMO', 'MN1', 'MN2', 'MN7', 'MN8', 'MND', 'MNL',
                               'MNV', 'MOD', 'MP8', 'MPH', 'MPJ', 'MPQ', 'MPR', 'MSA', 'MSE', 'MSL', 'MSO', 'MSP', 'MT2', 'MVA', 'MYN', 'MYR', 'N10', 'N2C', 'N7P',
                               'NA8', 'NAL', 'NAM', 'NB8', 'NBQ', 'NC1', 'NCB', 'NCY', 'NDF', 'NEM', 'NEP', 'NFA', 'NH2', 'NHL', 'NIT', 'NIY', 'NKS', 'NLE', 'NLN',
                               'NLO', 'NLP', 'NLQ', 'NMC', 'NME', 'NMM', 'NNH', 'NOT', 'NPH', 'NPI', 'NRP', 'NRQ', 'NSK', 'NTY', 'NVA', 'NWD', 'NYB', 'NYC', 'NYG',
                               'NYS', 'NZH', 'O12', 'OAS', 'OBF', 'OBS', 'OCS', 'OCY', 'ODA', 'ODS', 'OEM', 'OHI', 'OHS', 'OIC', 'OLE', 'OLZ', 'OMH', 'OMT', 'OMX',
                               'OMY', 'OMZ', 'ONH', 'ONL', 'OPR', 'ORD', 'ORN', 'ORQ', 'OSE', 'OTH', 'OTY', 'OXX', 'OZT', 'P1L', 'P2Q', 'P2Y', 'P3Q', 'PAA', 'PAQ',
                               'PAS', 'PAT', 'PAU', 'PBB', 'PBF', 'PBI', 'PCA', 'PCC', 'PCE', 'PCS', 'PDD', 'PDL', 'PEC', 'PF5', 'PFF', 'PFX', 'PG1', 'PG9', 'PGA',
                               'PGL', 'PGY', 'PHA', 'PHD', 'PHI', 'PHL', 'PHM', 'PIA', 'PLE', 'PM3', 'PN2', 'PNL', 'PNZ', 'POM', 'PPN', 'PR3', 'PR4', 'PR7', 'PR9',
                               'PRJ', 'PRK', 'PRQ', 'PRR', 'PRS', 'PRV', 'PSA', 'PSH', 'PSS', 'PSW', 'PTA', 'PTH', 'PTL', 'PTM', 'PTR', 'PVH', 'PVL', 'PXZ', 'PYA',
                               'PYT', 'PYX', 'QDS', 'QIL', 'QLG', 'QMM', 'QPH', 'QUA', 'Q000', 'Q004', 'Q005', 'Q00B', 'Q00C', 'Q01W', 'Q0A0', 'Q0A1', 'Q0A2', 'Q0A8',
                               'Q0A9', 'Q0AA', 'Q0AB', 'Q0AC', 'Q0AF', 'Q0AG', 'Q0AH', 'Q0AK', 'Q0AR', 'Q0AY', 'Q0AZ', 'Q0BN', 'Q0CS', 'Q0FL', 'Q0NC', 'Q0YG', 'Q143',
                               'Q175', 'Q192', 'Q193', 'Q1AC', 'Q1LU', 'Q1PA', 'Q1PI', 'Q1TQ', 'Q1TY', 'Q1ZN', 'Q200', 'Q23F', 'Q23S', 'Q26B', 'Q2AD', 'Q2AG', 'Q2AO',
                               'Q2CO', 'Q2FM', 'Q2HF', 'Q2LU', 'Q2ML', 'Q2MR', 'Q2MT', 'Q2P0', 'Q2PI', 'Q2PP', 'Q2RA', 'Q2TL', 'Q2TY', 'Q2VA', 'Q2XA', 'Q32S', 'Q32T',
                               'Q3AH', 'Q3AR', 'Q3CF', 'Q3FG', 'Q3GA', 'Q3GL', 'Q3MD', 'Q3MM', 'Q3MY', 'Q3NF', 'Q3O3', 'Q3PA', 'Q3QN', 'Q3TY', 'Q3XH', 'Q4AF', 'Q4BA',
                               'Q4CF', 'Q4CY', 'Q4DB', 'Q4DP', 'Q4F3', 'Q4FB', 'Q4FW', 'Q4HT', 'Q4IN', 'Q4MM', 'Q4PH', 'Q4U7', 'Q5AB', 'Q5CS', 'Q5HP', 'Q5OH', 'Q5PG',
                               'Q5ZA', 'Q6CL', 'Q6CW', 'Q6HC', 'Q6HG', 'Q6HN', 'Q6HT', 'Q7HA', 'Q7JA', 'Q7MN', 'Q9DN', 'Q9DS', 'Q9NE', 'Q9NF', 'Q9NR', 'Q9NV', 'R1A',
                               'R1B', 'R1F', 'R2P', 'R4K', 'R7A', 'RC7', 'RCY', 'RE0', 'RE3', 'RGL', 'ROP', 'S1H', 'S2C', 'S2D', 'S2P', 'SAC', 'SAH', 'SAR', 'SBD',
                               'SBL', 'SCH', 'SCS', 'SCY', 'SD2', 'SDP', 'SEB', 'SEC', 'SEE', 'SEG', 'SEL', 'SEM', 'SEN', 'SEP', 'SET', 'SGB', 'SHC', 'SHP', 'SHR',
                               'SIB', 'SIC', 'SIN', 'SLA', 'SLR', 'SLZ', 'SMC', 'SME', 'SMF', 'SNC', 'SNN', 'SOC', 'SOY', 'STA', 'STY', 'SUB', 'SUI', 'SUN', 'SVA',
                               'SVV', 'SVW', 'SVX', 'SVZ', 'SXE', 'SYS', 'T0I', 'T11', 'T16', 'T29', 'T66', 'TA4', 'TAV', 'TBG', 'TBM', 'TCQ', 'TCR', 'TDD', 'TEF',
                               'TFQ', 'TH5', 'TH6', 'THC', 'THO', 'THZ', 'TIH', 'TIS', 'TLX', 'TMB', 'TMD', 'TNB', 'TNR', 'TNY', 'TOX', 'TPH', 'TPL', 'TPN', 'TPO',
                               'TPQ', 'TQI', 'TQQ', 'TRF', 'TRG', 'TRN', 'TRO', 'TRQ', 'TRW', 'TRX', 'TS9', 'TSI', 'TST', 'TTQ', 'TTS', 'TY2', 'TY3', 'TY5', 'TY8',
                               'TY9', 'TYB', 'TYI', 'TYJ', 'TYN', 'TYO', 'TYQ', 'TYS', 'TYT', 'TYX', 'TYY', 'TZB', 'TZO', 'UAL', 'UDS', 'UM1', 'UM2', 'UMA', 'UN1',
                               'UN2', 'V1A', 'VAD', 'VAF', 'VAH', 'VAS', 'VB1', 'VDL', 'VLL', 'VLM', 'VME', 'VMS', 'VOL', 'VR0', 'WFP', 'X2W', 'X9Q', 'XAA', 'XBB',
                               'XCP', 'XPC', 'XPR', 'XSN', 'XX1', 'XXA', 'XXY', 'XYG', 'Y28', 'YCM', 'YCP', 'YNM', 'YOF', 'YYA', 'Z3E', 'Z70', 'ZAE', 'ZAL', 'ZBZ', 
                               'ZFB', 'ZUK', 'ZYJ', 'ZYK', 'ZZD', 'ZZJ', 'ZZU']
        
        cofactors = ['12P', '144', '15P', '16D', '1BO', '1CP', '1PE', '1PS', '2OS', '2PE', '7PE', 'ACN', 'ACO', 'AGC', 'B12', 'B3P', 'B7G', 'BCB', 'BCL', 'BCN', 'BE7',
                     'BEQ', 'BGC', 'BMA', 'BME', 'BNG', 'BOG', 'BPB', 'BPH', 'BTB', 'BU1', 'BU2', 'BU3', 'BUD', 'C10', 'C15', 'C8E', 'CAA', 'CAO', 'CBM', 'CCN', 'CDP',
                     'CE1', 'CHL', 'CIT', 'CL1', 'CL2', 'CLA', 'CLN', 'CM', 'CM5', 'CN', 'CNC', 'COA', 'COB', 'COH', 'COY', 'CP3', 'CPS', 'CRY', 'CTP', 'CXE', 'DDH', 'DDQ',
                     'DEU', 'DHD', 'DHE', 'DIA', 'DIO', 'DMF', 'DMS', 'DMU', 'DMX', 'DOX', 'DR6', 'DTT', 'DTV', 'DXG', 'EDO', 'EEE', 'EGL', 'EOH', 'EPE', 'ETF', 'FAD',
                     'FLC', 'FMN', 'FMT', 'FRU', 'FUC', 'GBL', 'GCD', 'GDP', 'GLO', 'GMP', 'GNP', 'GOL', 'GTP', 'H4B', 'HAS', 'HCO', 'HDD', 'HDM', 'HEA', 'HEB', 'HEC',
                     'HEG', 'HEM', 'HEO', 'HES', 'HEZ', 'HIF', 'HNI', 'HTG', 'HTO', 'ICI', 'ICT', 'IDT', 'IMD', 'IOH', 'IPA', 'IPH', 'JEF', 'LAK', 'LAT', 'LBT', 'LDA',
                     'LMT', 'MA4', 'MCA', 'MLA', 'MES', 'MG8', 'MHA', 'MMP', 'MOH', 'MP1', 'MPD', 'MPO', 'MRD', 'MRY', 'MYR', 'MTL', 'N8E', 'NAD', 'NAG', 'NAI', 'NAP',
                     'NDG', 'NDO', 'NDP', 'NHE', 'OCT' ,'OTE', 'P33', 'P4C', 'P6G', 'PC3', 'PCU', 'PDO', 'PE4', 'PE7', 'PE8', 'PEG', 'PG4', 'PG5', 'PG6', 'PGE', 'PGO',
                     'PGQ', 'PGR', 'PIG', 'PIN', 'PNI', 'PSE', 'POL', 'POR', 'PP9', 'Q12P', 'Q144', 'Q15P', 'Q16D', 'Q1BO', 'Q1CP', 'Q1PE', 'Q1PS', 'Q2OS', 'Q2PE', 'Q7PE',
                     'SAL', 'SBT', 'SCA', 'SDS', 'SGM', 'SOR', 'SPD', 'SPK', 'SPM', 'SRM', 'SRT', 'SUC', 'TAM', 'TAR', 'TAU', 'TBU', 'TLA', 'tmp', 'TRE', 'TRS', 'TRT',
                     'U10', 'UMP', 'UMQ', 'UNK', 'URE', 'XPE', 'ZEM', 'ZNH']
        
        ions = ['1CU', '2HP', '2OF', '3CO', '3MT', '3NI', '4MO', '543', '6MO', 'ACT', 'AG', 'AL', 'ATH', 'AU', 'AU3', 'AZI', 'BA', 'BCT', 'BF4', 'BO4', 'BR', 'CA', 'CD',
                'CD1', 'CD3', 'CD5', 'CE', 'CHT', 'CL', 'CLO', 'CO', 'CO3', 'CO5', 'CON', 'CR', 'CS', 'CU', 'CU1', 'CUA', 'CUZ', 'CYN', 'DMI', 'E4N', 'EU', 'EU3', 'F', 'FE',
                'FE2', 'FLO', 'FPO', 'GA', 'GD3', 'HAI', 'HG', 'HO', 'IN', 'IR', 'IR3', 'IRI', 'K', 'KO4', 'LA', 'LCO', 'LI', 'LU', 'MG', 'MH2', 'MH3', 'MLI', 'MLT', 'MN',
                'MN3', 'NA', 'NC', 'NET', 'NH4', 'NI', 'NO2', 'NO3', 'OAA', 'OH', 'OS', 'OXL', 'PB', 'PD', 'PER', 'PI', 'PO3', 'PO4', 'PR', 'PT', 'PT4', 'PTN', 'Q1CU', 'Q2HP',
                'Q2OF', 'Q3CO', 'Q3MT', 'Q3NI', 'Q4MO', 'Q543', 'Q6MO', 'RB', 'RHD', 'RU', 'SB', 'SCN', 'SM', 'SMO', 'SO3', 'SO4', 'SOH', 'SUL', 'SR', 'TB', 'TCN', 'TEA',
                'THE', 'TL', 'TMA', 'TRA', 'V', 'W', 'Y1', 'YB', 'YT3', 'ZN', 'ZN2', 'ZN3', 'ZNO']
        
        organometallics = ['1CM', '2BM', '2MO', '6HE', '6WO', '749', '7HE', 'A48', 'A71', 'A72', 'AF3', 'AIO', 'ALF', 'AMS', 'APW', 'ARS', 'AST', 'AUC', 'B69', 'B70',
                           'BEF', 'BF2', 'BH1', 'BPT', 'BRO', 'C1O', 'C2C', 'CAC', 'CAD', 'CFM', 'CFN', 'CFO', 'CFQ', 'CLF', 'CLP', 'CM6', 'CN1', 'CNB', 'CNF', 'CUB',
                           'CUM', 'CUN', 'CUO', 'DP4', 'DRN', 'DW2', 'DWC', 'EMC', 'EMT', 'F3S', 'FCO', 'FDC', 'FEO', 'FES', 'FNE', 'FS1', 'FS2', 'FS3', 'FS4', 'FSO',
                           'HC0', 'HC1', 'HF3', 'HF5', 'HGB', 'HGC', 'HGI', 'HO3', 'IAP', 'IBZ', 'IDO', 'IOD', 'ISU', 'IUM', 'JM1', 'KEG', 'KYS', 'LCP', 'LPT', 'MAC',
                           'MAP', 'MBO', 'MGF', 'MGO', 'MM4', 'MMC', 'MN5', 'MN6', 'MNH', 'MNR', 'MO1', 'MO2', 'MO3', 'MO4', 'MO5', 'MO6', 'MO7', 'MOO', 'MOS', 'MTD',
                           'MW1', 'MW2', 'MW3', 'NA2', 'NA5', 'NA6', 'NAO', 'NAW', 'NBF', 'NCO', 'ND4', 'NFS', 'NI1', 'NI2', 'NI3', 'NPB', 'O4M', 'OC1', 'OC2', 'OC3',
                           'OC4', 'OC5', 'OC6', 'OC7', 'OC8', 'OCL', 'OCM', 'OCN', 'OCO', 'OF1', 'OF2', 'OF3', 'PBA', 'PBC', 'PBM', 'PCL', 'PEU', 'PHF', 'PHG', 'PMB',
                           'PPB', 'Q1CM', 'Q2BM', 'Q2MO', 'Q6HE', 'Q6WO', 'Q749', 'Q7HE', 'R4A', 'R6A', 'REO', 'RHX', 'RTC', 'SBO', 'SE4', 'SES', 'SF3', 'SF4', 'SIF',
                           'SM2', 'SM3', 'SM4', 'SRB', 'SRD', 'SUF', 'T42', 'TAS', 'UNX', 'V35', 'V36', 'V7O', 'VA1', 'VO4', 'WCC', 'WO4', 'WO5', 'XCC', 'YBT', 'ZO3']
        
        water = ['HOH', 'WAT', 'TIP', 'SOL', 'OH2', 'DOD', 'D20']

        def is_small_molecule(residue_name):
            # Check if the residue name is not in any of the other residue_names sets
            return (
                residue_name not in water and
                residue_name not in protein and
                residue_name not in ions and
                residue_name not in cofactors and
                residue_name not in organometallics
            )

        def extract_and_save(input_pdb_file, output_dir):
            # Get the base filename without the extension
            base_filename = os.path.splitext(os.path.basename(input_pdb_file))[0]

            # Create a new output directory for this PDB file
            pdb_output_dir = os.path.join(output_dir, base_filename)

            # Create the output directory if it doesn't exist
            if not os.path.exists(pdb_output_dir):
                os.makedirs(pdb_output_dir)

            # Lists to store lines for each category
            solvent_lines = []
            protein_lines = []
            ion_lines = []
            cofactors_lines = []
            organometallics_lines = []
            small_molecule_lines = []

            with open(input_pdb_file, 'r') as input_file:
                for line in input_file:
                    if line.startswith('ATOM') or line.startswith('HETATM'):
                        residue_name = line[17:20].strip()
                        if residue_name in water:
                            solvent_lines.append(line)
                        elif residue_name in protein:
                            protein_lines.append(line)
                        elif residue_name in alternative_aminoacids:
                            protein_lines.append(line)
                        elif residue_name in ions:
                            ion_lines.append(line)
                        elif residue_name in cofactors:
                            cofactors_lines.append(line)
                        elif residue_name in organometallics:
                            organometallics_lines.append(line)
                        elif is_small_molecule(residue_name):
                            small_molecule_lines.append(line)

            # Write lines to respective output files if they are not empty
            if solvent_lines:
                with open(os.path.join(pdb_output_dir, f'{base_filename}_solvent.pdb'), 'w') as f:
                    f.writelines(solvent_lines)
            if protein_lines:
                with open(os.path.join(pdb_output_dir, f'{base_filename}_protein.pdb'), 'w') as f:
                    f.writelines(protein_lines)
            if ion_lines:
                with open(os.path.join(pdb_output_dir, f'{base_filename}_ion.pdb'), 'w') as f:
                    f.writelines(ion_lines)
            if cofactors_lines:
                with open(os.path.join(pdb_output_dir, f'{base_filename}_cofactors.pdb'), 'w') as f:
                    f.writelines(cofactors_lines)
            if organometallics_lines:
                with open(os.path.join(pdb_output_dir, f'{base_filename}_organometallics.pdb'), 'w') as f:
                    f.writelines(organometallics_lines)
            if small_molecule_lines:
                with open(os.path.join(pdb_output_dir, f'{base_filename}_small_molecule.pdb'), 'w') as f:
                    f.writelines(small_molecule_lines)

        def process_pdb_files_in_directory(directory, output_directory):
            pdb_files = [file for file in os.listdir(directory) if file.endswith('.pdb')]
            for pdb_file in pdb_files:
                pdb_file_path = os.path.join(directory, pdb_file)
                extract_and_save(pdb_file_path, output_directory)

        pdb_directory = os.sep.join([settings.DATA_DIR, 'protonated_models'])
        output_directory = os.sep.join([settings.DATA_DIR, 'protonated_models_split'])
        if not os.path.exists(output_directory):
            os.makedirs(output_directory)

        process_pdb_files_in_directory(pdb_directory, output_directory)
        
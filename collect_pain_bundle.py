#!/usr/bin/env python3
"""
Pain Radar OS · Cloud Multi-Dimensional Collector (Zero External Dependencies)
==============================================================================
Designed to run cleanly on GitHub Actions (free US Azure runner) or any machine
with direct access to official biomedical APIs:

  1. ClinicalTrials.gov API v2    -> Global & China Pain/Analgesia clinical trials
  2. openFDA Drugs@FDA API        -> US FDA benchmark analgesia approvals
  3. openFDA FAERS API            -> Real-world adverse events (2.2类改良 pain points)
  4. Open Targets GraphQL API     -> Pain target genetics, druggability & clinical drugs
  5. RCSB PDB Search API v2       -> Experimental 3D protein structure counts per target
  6. ChEMBL REST API (EMBL-EBI)   -> Molecule pharmacology, max phase & mechanisms
  7. Europe PMC REST API          -> High-impact Pain literature + linked NCT IDs

Outputs:
  - data/pain_radar_bundle.json     (plain JSON bundle for local Pain Radar OS import)
  - data/pain_radar_bundle.json.gz  (compressed bundle for fast GitHub Raw pull)
  - data/sync_summary.md            (human-readable run summary)

Usage:
  python collect_pain_bundle.py --mode full
  python collect_pain_bundle.py --mode incremental --lookback-days 14
  python collect_pain_bundle.py --self-test
"""

import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timedelta, timezone

BUNDLE_VERSION = "2.0"
USER_AGENT = "PainRadarOS-CloudSync/2.0 (Open Biomedical Intelligence Research)"
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_OUT_DIR = os.path.join(BASE_DIR, "data")


def utc_now_iso():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")


def utc_today_str():
    return datetime.now(timezone.utc).strftime("%Y-%m-%d")


# ==================== HTTP HELPERS WITH RETRY ====================

def http_request_json(url, method="GET", payload=None, headers=None, timeout=25, retries=3):
    req_headers = {
        "User-Agent": USER_AGENT,
        "Accept": "application/json",
    }
    if headers:
        req_headers.update(headers)
    data = None
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        req_headers.setdefault("Content-Type", "application/json")

    last_err = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers=req_headers, method=method)
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                if resp.status == 204:
                    return None
                raw = resp.read()
                if not raw:
                    return None
                return json.loads(raw.decode("utf-8"))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            last_err = e
            if e.code in (429, 500, 502, 503, 504) and attempt < retries:
                time.sleep(1.5 * attempt)
                continue
            raise
        except Exception as e:
            last_err = e
            if attempt < retries:
                time.sleep(1.5 * attempt)
                continue
            raise last_err
    return None


# ==================== PAIN DOMAIN DEFINITIONS ====================

PAIN_TARGET_RULES = [
    ("NaV1.8 / 钠通道", r"nav\s*1\.?8|nav\s*1\.?7|nav\s*1\.?9|scn10a|scn9a|scn11a|suzetrigine|vx-?548|vx-?993|ltg-?001|stc-?004|hrs-?4800|hrs-?2129|jmkx000623|journavx|钠离子通道|钠通道"),
    ("外周KOR / 偏向性MOR", r"oliceridine|trv-?130|奥赛利定|泰吉利定|shr-?8554|安瑞克芬|anrikefon|hsk21542|difelikefalin|oprk1|oprm1|偏向性|外周.*阿片|κ阿片|kappa.*opioid|biased.*mor"),
    ("长效局麻 / 2类改良NSAID", r"bupivacaine|ropivacaine|meloxicam|ketorolac|celecoxib|布比卡因|罗哌卡因|美洛昔康|脂质体|纳米晶|微球|缓释|exparel|htx-?011|zynrelef|普坦宁|ptgs2|cox-?2|505\(b\)\(2\)|2\.2类"),
    ("CGRP / 偏头痛", r"cgrp|calcrl|ramp1|rimegepant|ubrogepant|atogepant|erenumab|galcanezumab|fremanezumab|eptinezumab|zavegepant|瑞美吉泮|依瑞奈尤|瑞玛奈珠|偏头痛|migraine"),
    ("TRPV1 / NGF / P2X3 / 新靶点", r"trpv1|trpa1|ngf|ntrk1|tanezumab|fasinumab|p2x3|p2rx3|kv7|kcnq2|sigma-?1|nk1|adra2a|dexmedetomidine|右美托咪定|resiniferatoxin|capsaicin|辣椒素|肉毒毒素"),
    ("阿片改良 / 镇痛复方", r"oxycodone|hydromorphone|tapentadol|tramadol|fentanyl|morphine|sufentanil|pregabalin|gabapentin|duloxetine|羟考酮|氢吗啡酮|他喷他多|曲马多|芬太尼|吗啡|舒芬太尼|纳布啡|地佐辛|透皮贴|鼻喷"),
]

# Non-drug / off-target noise filter for ClinicalTrials.gov
NOISE_TRIAL_PATTERN = re.compile(
    r"\b(acupuncture|massage|yoga|mindfulness|cognitive behavioral|virtual reality|music therapy|physiotherapy|physical therapy|chiropractic|transcranial magnetic|tens\b|spinal cord stimulat|nerve stimulat|radiofrequency ablation|cryoablation)\b",
    re.IGNORECASE,
)


def classify_pain_target(text):
    hay = (text or "").lower()
    for label, pattern in PAIN_TARGET_RULES:
        if re.search(pattern, hay, re.IGNORECASE):
            return label
    return "其他镇痛机制"


# Core Pain Targets (for Open Targets + RCSB PDB)
PAIN_TARGETS = [
    {"gene_symbol": "SCN10A", "ensembl_id": "ENSG00000185313", "target_group": "NaV1.8 / 钠通道"},
    {"gene_symbol": "SCN9A",  "ensembl_id": "ENSG00000169432", "target_group": "NaV1.8 / 钠通道"},
    {"gene_symbol": "SCN11A", "ensembl_id": "ENSG00000168356", "target_group": "NaV1.8 / 钠通道"},
    {"gene_symbol": "OPRK1",  "ensembl_id": "ENSG00000082556", "target_group": "外周KOR / 偏向性MOR"},
    {"gene_symbol": "OPRM1",  "ensembl_id": "ENSG00000112038", "target_group": "外周KOR / 偏向性MOR"},
    {"gene_symbol": "CALCRL", "ensembl_id": "ENSG00000064989", "target_group": "CGRP / 偏头痛"},
    {"gene_symbol": "RAMP1",  "ensembl_id": "ENSG00000132329", "target_group": "CGRP / 偏头痛"},
    {"gene_symbol": "TRPV1",  "ensembl_id": "ENSG00000196689", "target_group": "TRPV1 / NGF / P2X3 / 新靶点"},
    {"gene_symbol": "TRPA1",  "ensembl_id": "ENSG00000104321", "target_group": "TRPV1 / NGF / P2X3 / 新靶点"},
    {"gene_symbol": "NTRK1",  "ensembl_id": "ENSG00000198400", "target_group": "TRPV1 / NGF / P2X3 / 新靶点"},
    {"gene_symbol": "P2RX3",  "ensembl_id": "ENSG00000109991", "target_group": "TRPV1 / NGF / P2X3 / 新靶点"},
    {"gene_symbol": "KCNQ2",  "ensembl_id": "ENSG00000075043", "target_group": "TRPV1 / NGF / P2X3 / 新靶点"},
    {"gene_symbol": "PTGS2",  "ensembl_id": "ENSG00000073756", "target_group": "长效局麻 / 2类改良NSAID"},
    {"gene_symbol": "ADRA2A", "ensembl_id": "ENSG00000150594", "target_group": "TRPV1 / NGF / P2X3 / 新靶点"},
]

# ClinicalTrials.gov curated multi-slice queries
CTGOV_SLICES = [
    {
        "name": "nav18_sodium",
        "cond": "Pain OR Neuropathic Pain OR Postoperative Pain OR Radiculopathy",
        "intr": "NaV1.8 OR NaV1.7 OR Suzetrigine OR VX-548 OR VX-993 OR LTG-001 OR STC-004 OR HRS4800 OR HRS-2129 OR JMKX000623",
        "locn": "",
        "extra_adv": "",
        "max_full": 300,
        "max_inc": 120,
    },
    {
        "name": "kor_biased_mor",
        "cond": "Pain OR Analgesia OR Postoperative Pain OR Pruritus",
        "intr": "Oliceridine OR TRV130 OR Tegileridine OR SHR8554 OR Anrikefon OR HSK21542 OR Difelikefalin OR CR845",
        "locn": "",
        "extra_adv": "",
        "max_full": 200,
        "max_inc": 80,
    },
    {
        "name": "postop_long_acting",
        "cond": "Postoperative Pain OR Acute Pain OR Surgical Pain",
        "intr": "Bupivacaine OR Ropivacaine OR Meloxicam OR Ketorolac OR Celecoxib OR HTX-011 OR EXPAREL",
        "locn": "",
        "extra_adv": "",
        "max_full": 300,
        "max_inc": 120,
    },
    {
        "name": "cgrp_migraine",
        "cond": "Migraine",
        "intr": "Rimegepant OR Ubrogepant OR Atogepant OR Zavegepant OR Erenumab OR Galcanezumab OR Fremanezumab OR Eptinezumab",
        "locn": "",
        "extra_adv": "",
        "max_full": 250,
        "max_inc": 100,
    },
    {
        "name": "pain_china",
        "cond": "Pain OR Analgesia OR Postoperative Pain OR Neuropathic Pain OR Migraine",
        "intr": "",
        "locn": "China",
        "extra_adv": "",
        "max_full": 350,
        "max_inc": 150,
    },
    {
        "name": "pain_terminated_watchlist",
        "cond": "Pain OR Analgesia OR Postoperative Pain OR Neuropathic Pain OR Migraine",
        "intr": "",
        "locn": "",
        "extra_adv": "AREA[OverallStatus](TERMINATED OR SUSPENDED)",
        "max_full": 250,
        "max_inc": 100,
    },
]

# Benchmark Pain Ingredients for openFDA Drugs@FDA
FDA_PAIN_INGREDIENTS = [
    "SUZETRIGINE", "OLICERIDINE", "DIFELIKEFALIN",
    "BUPIVACAINE", "ROPIVACAINE", "MELOXICAM", "KETOROLAC", "CELECOXIB",
    "RIMEGEPANT", "UBROGEPANT", "ATOGEPANT", "ZAVEGEPANT",
    "ERENUMAB", "GALCANEZUMAB", "FREMANEZUMAB", "EPTINEZUMAB",
    "TAPENTADOL", "OXYCODONE", "HYDROMORPHONE", "SUFENTANIL",
    "DEXMEDETOMIDINE", "CAPSAICIN", "PREGABALIN", "GABAPENTIN", "DULOXETINE",
]

# Benchmark Pain Drugs for openFDA FAERS (2.2类改良痛点挖掘)
FAERS_PAIN_DRUGS = [
    "MELOXICAM", "BUPIVACAINE", "ROPIVACAINE", "KETOROLAC", "CELECOXIB",
    "OXYCODONE", "TRAMADOL", "FENTANYL", "OLICERIDINE", "SUZETRIGINE",
    "RIMEGEPANT", "ATOGEPANT", "PREGABALIN", "GABAPENTIN",
]

# Benchmark Molecules for ChEMBL
CHEMBL_SEED_NAMES = [
    "SUZETRIGINE", "OLICERIDINE", "DIFELIKEFALIN", "MELOXICAM", "BUPIVACAINE",
    "ROPIVACAINE", "KETOROLAC", "CELECOXIB", "RIMEGEPANT", "UBROGEPANT",
    "ATOGEPANT", "ZAVEGEPANT", "ERENUMAB", "GALCANEZUMAB", "FREMANEZUMAB",
    "EPTINEZUMAB", "TAPENTADOL", "PREGABALIN", "GABAPENTIN", "DULOXETINE",
    "CAPSAICIN", "DEXMEDETOMIDINE", "TANEZUMAB", "GEFAPIXANT",
]

# Curated High-Impact Pain Literature Queries (Europe PMC)
EPMC_QUERIES = [
    '("NaV1.8" OR "SCN10A" OR "suzetrigine" OR "VX-548" OR "LTG-001") AND (pain OR analgesia) AND PUB_YEAR:[2022 TO 2027]',
    '("kappa opioid" OR "anrikefon" OR "HSK21542" OR "difelikefalin" OR "oliceridine" OR "tegileridine" OR "biased mu-opioid") AND (pain OR analgesia) AND PUB_YEAR:[2022 TO 2027]',
    '("liposomal bupivacaine" OR "meloxicam" OR "postoperative pain") AND (JOURNAL:"N Engl J Med" OR JOURNAL:"Lancet" OR JOURNAL:"JAMA" OR JOURNAL:"Anesthesiology" OR JOURNAL:"Pain") AND PUB_YEAR:[2022 TO 2027]',
    '("CGRP" OR "rimegepant" OR "atogepant" OR "ubrogepant") AND migraine AND (JOURNAL:"N Engl J Med" OR JOURNAL:"Lancet" OR JOURNAL:"JAMA" OR JOURNAL:"Lancet Neurol" OR JOURNAL:"Cephalalgia") AND PUB_YEAR:[2023 TO 2027]',
]


# ==================== NORMALIZERS (TESTABLE WITHOUT NETWORK) ====================

def normalize_ctgov_study(st):
    ps = st.get("protocolSection") or {}
    ident = ps.get("identificationModule") or {}
    status = ps.get("statusModule") or {}
    design = ps.get("designModule") or {}
    cond = ps.get("conditionsModule") or {}
    arms = ps.get("armsInterventionsModule") or {}
    spon = ps.get("sponsorCollaboratorsModule") or {}
    locs = (ps.get("contactsLocationsModule") or {}).get("locations") or []
    outs = (ps.get("outcomesModule") or {}).get("primaryOutcomes") or []

    trial_id = (ident.get("nctId") or "").strip().upper()
    title = (ident.get("briefTitle") or ident.get("officialTitle") or "").strip()
    if not trial_id or not title:
        return None

    iv_list = arms.get("interventions") or []
    drug_ivs = [
        iv for iv in iv_list
        if (iv.get("type") or "").upper() in ("DRUG", "BIOLOGICAL", "COMBINATION_PRODUCT", "")
    ]
    interventions = [
        f"{iv.get('name', '')} ({iv.get('type', '')})".strip()
        for iv in (drug_ivs or iv_list)
        if iv.get("name")
    ]
    conditions_str = "; ".join(cond.get("conditions") or [])
    interventions_str = "; ".join(interventions[:6])

    # Filter out obvious non-drug device/physical-therapy noise if no recognized pain drug keyword
    combined = f"{title} {conditions_str} {interventions_str}"
    if NOISE_TRIAL_PATTERN.search(combined) and classify_pain_target(combined) == "其他镇痛机制":
        return None

    countries = sorted({(loc.get("country") or "").strip() for loc in locs if loc.get("country")})
    primary_outs = [o.get("measure", "") for o in outs if o.get("measure")]
    phases = design.get("phases") or []

    return {
        "trial_id": trial_id,
        "registry": "ctgov",
        "brief_title": title,
        "conditions": conditions_str,
        "interventions": interventions_str,
        "phase": "/".join(phases) if phases else "NA",
        "overall_status": status.get("overallStatus") or "",
        "why_stopped": status.get("whyStopped") or "",
        "sponsor": (spon.get("leadSponsor") or {}).get("name") or "",
        "countries": ", ".join(countries),
        "enrollment": int((design.get("enrollmentInfo") or {}).get("count") or 0),
        "primary_outcomes": "; ".join(primary_outs[:3]),
        "start_date": (status.get("startDateStruct") or {}).get("date") or "",
        "primary_completion_date": (status.get("primaryCompletionDateStruct") or {}).get("date") or "",
        "has_results": 1 if st.get("hasResults") else 0,
        "last_update_date": (status.get("lastUpdatePostDateStruct") or {}).get("date") or utc_today_str(),
        "source_tag": "cloud_ctgov_v2",
    }


def normalize_openfda_approval(item):
    app_num = (item.get("application_number") or "").strip().upper()
    if not app_num:
        return None
    sponsor = (item.get("sponsor_name") or "").strip()
    prods = item.get("products") or []
    p0 = prods[0] if prods else {}
    brand = (p0.get("brand_name") or "").strip()
    dosage = (p0.get("dosage_form") or "").strip()
    route = (p0.get("route") or "").strip()
    active_list = p0.get("active_ingredients") or []
    ing = ", ".join((a.get("name") or "").strip() for a in active_list if a.get("name"))
    strength = "; ".join((a.get("strength") or "").strip() for a in active_list if a.get("strength"))[:80]

    subs = item.get("submissions") or []
    orig_date = ""
    latest_date = ""
    for s in subs:
        if (s.get("submission_status") or "").upper() != "AP":
            continue
        d = (s.get("submission_status_date") or "").strip()
        if len(d) == 8 and d.isdigit():
            fmt = f"{d[:4]}-{d[4:6]}-{d[6:8]}"
            if fmt > latest_date:
                latest_date = fmt
            if (s.get("submission_type") or "").upper() == "ORIG":
                if not orig_date or fmt < orig_date:
                    orig_date = fmt

    reg_prefix = re.sub(r"\d+$", "", app_num) or app_num[:3]
    drug_label = f"{ing} ({brand})" if (brand and ing) else (brand or ing or app_num)
    dosage_full = f"{dosage} / {route}".strip(" /")

    return {
        "approval_no": f"FDA-{app_num}",
        "drug_name_cn": drug_label,
        "drug_name_en": ing or brand,
        "brand_name": brand,
        "drug_category": "FDA批准药品",
        "dosage_form": dosage_full,
        "strength": strength,
        "mah_holder": sponsor,
        "manufacturer": sponsor,
        "approval_date": orig_date or latest_date,
        "is_imported": 1,
        "reg_class": f"FDA {reg_prefix}",
        "indication_detail": f"openFDA Drugs@FDA: {app_num} · 首次批准 {orig_date or 'N/A'} · 最近批准 {latest_date or 'N/A'} · 持有人 {sponsor}",
        "tags": f"openFDA,{reg_prefix},云端同步",
        "source_tag": "cloud_openfda",
    }


def normalize_epmc_article(art):
    pmid = str(art.get("pmid") or "").strip()
    pmcid = str(art.get("pmcid") or "").strip()
    pub_id = f"PMID:{pmid}" if pmid else (f"PMCID:{pmcid}" if pmcid else "")
    title = (art.get("title") or "").strip().rstrip(".")
    if not pub_id or not title:
        return None

    abstract = re.sub(r"<[^>]+>", " ", art.get("abstractText") or "")
    abstract = re.sub(r"\s+", " ", abstract).strip()
    nct_hits = sorted(set(re.findall(r"NCT\d{8}", f"{title} {abstract}", re.IGNORECASE)))
    nct_ids = ", ".join(x.upper() for x in nct_hits[:8])

    journal = (art.get("journalTitle") or "").strip()
    if not journal:
        ji = (art.get("journalInfo") or {}).get("journal") or {}
        journal = (ji.get("medlineAbbreviation") or ji.get("title") or "").strip()

    doi = (art.get("doi") or "").strip()
    url = f"https://doi.org/{doi}" if doi else (
        f"https://europepmc.org/article/MED/{pmid}" if pmid else f"https://europepmc.org/article/PMC/{pmcid}"
    )
    pub_year = int(art.get("pubYear") or 0) if str(art.get("pubYear") or "").isdigit() else 0

    return {
        "pub_id": pub_id,
        "title": title,
        "journal": journal,
        "pub_year": pub_year,
        "pub_date": (art.get("firstPublicationDate") or "").strip(),
        "authors": (art.get("authorString") or "")[:180],
        "cited_by_count": int(art.get("citedByCount") or 0),
        "doi": doi,
        "nct_ids": nct_ids,
        "abstract_snippet": abstract[:420],
        "url": url,
        "source_tag": "cloud_europepmc",
    }


# ==================== COLLECTORS ====================

def collect_ctgov(mode="full", lookback_days=14):
    trials_by_id = {}
    since_date = (datetime.now(timezone.utc) - timedelta(days=lookback_days)).strftime("%Y-%m-%d")
    slice_stats = []

    for sl in CTGOV_SLICES:
        limit = sl["max_full"] if mode == "full" else sl["max_inc"]
        adv_parts = ["AREA[InterventionType](DRUG OR BIOLOGICAL)"]
        if sl["extra_adv"]:
            adv_parts.append(sl["extra_adv"])
        if mode == "incremental":
            adv_parts.append(f"AREA[LastUpdatePostDate]RANGE[{since_date},MAX]")

        params = {
            "query.cond": sl["cond"],
            "filter.advanced": " AND ".join(adv_parts),
            "pageSize": str(min(200, limit)),
            "sort": "LastUpdatePostDate:desc",
            "countTotal": "true",
        }
        if sl["intr"]:
            params["query.intr"] = sl["intr"]
        if sl["locn"]:
            params["query.locn"] = sl["locn"]

        fetched = 0
        next_token = None
        total_remote = 0
        while fetched < limit:
            req_p = dict(params)
            if next_token:
                req_p["pageToken"] = next_token
            url = "https://clinicaltrials.gov/api/v2/studies?" + urllib.parse.urlencode(req_p)
            payload = http_request_json(url, timeout=25)
            if not payload:
                break
            if total_remote == 0:
                total_remote = int(payload.get("totalCount") or 0)
            studies = payload.get("studies") or []
            if not studies:
                break
            for st in studies:
                rec = normalize_ctgov_study(st)
                if rec:
                    trials_by_id[rec["trial_id"]] = rec
                fetched += 1
                if fetched >= limit:
                    break
            next_token = payload.get("nextPageToken")
            if not next_token:
                break
            time.sleep(0.25)

        slice_stats.append(f"{sl['name']}:{fetched}/{total_remote}")

    return list(trials_by_id.values()), f"{len(trials_by_id)} unique trials ({', '.join(slice_stats)})"


def collect_openfda_approvals():
    api_key = os.environ.get("OPENFDA_API_KEY", "").strip()
    approvals_by_no = {}

    # Batch ingredients in groups of 4 for clean URL length & high recall
    for i in range(0, len(FDA_PAIN_INGREDIENTS), 4):
        batch = FDA_PAIN_INGREDIENTS[i : i + 4]
        q_terms = " OR ".join(f'"{ing}"' for ing in batch)
        params = {
            "search": f"products.active_ingredients.name:({q_terms})",
            "limit": "100",
        }
        if api_key:
            params["api_key"] = api_key
        url = "https://api.fda.gov/drug/drugsfda.json?" + urllib.parse.urlencode(params)
        data = http_request_json(url, timeout=25)
        if not data:
            continue
        for item in data.get("results") or []:
            rec = normalize_openfda_approval(item)
            if rec:
                # Prioritize NDA / BLA approvals over duplicate generic ANDA noise (keep up to 2 ANDAs per ingredient if needed, always keep all NDA/BLA)
                is_nda_bla = rec["approval_no"].startswith(("FDA-NDA", "FDA-BLA"))
                if is_nda_bla or rec["approval_no"] not in approvals_by_no:
                    if is_nda_bla:
                        approvals_by_no[rec["approval_no"]] = rec
                    elif len([k for k in approvals_by_no if k.startswith("FDA-ANDA")]) < 35:
                        approvals_by_no[rec["approval_no"]] = rec
        time.sleep(0.3)

    return list(approvals_by_no.values()), f"{len(approvals_by_no)} FDA approvals"


def collect_openfda_faers():
    api_key = os.environ.get("OPENFDA_API_KEY", "").strip()
    signals = []

    for drug in FAERS_PAIN_DRUGS:
        base_search = f'patient.drug.openfda.generic_name:"{drug}"'
        p_total = {"search": base_search, "limit": "1"}
        if api_key:
            p_total["api_key"] = api_key
        meta_resp = http_request_json("https://api.fda.gov/drug/event.json?" + urllib.parse.urlencode(p_total), timeout=25)

        if not meta_resp:
            base_search = f'patient.drug.medicinalproduct:"{drug}"'
            p_total["search"] = base_search
            meta_resp = http_request_json("https://api.fda.gov/drug/event.json?" + urllib.parse.urlencode(p_total), timeout=25)

        if not meta_resp:
            continue

        total_reports = int(((meta_resp.get("meta") or {}).get("results") or {}).get("total") or 0)
        if total_reports <= 0:
            continue

        p_count = {
            "search": base_search,
            "count": "patient.reaction.reactionmeddrapt.exact",
        }
        if api_key:
            p_count["api_key"] = api_key
        cnt_resp = http_request_json("https://api.fda.gov/drug/event.json?" + urllib.parse.urlencode(p_count), timeout=25)
        if not cnt_resp:
            continue

        t_group = classify_pain_target(drug)
        for row in (cnt_resp.get("results") or [])[:12]:
            term = (row.get("term") or "").strip()
            cnt = int(row.get("count") or 0)
            if not term or cnt <= 0:
                continue
            share = round((cnt / total_reports) * 100.0, 2)
            signals.append({
                "signal_key": f"{drug.upper()}::{term.upper()}",
                "drug_name": drug.title(),
                "target_group": t_group,
                "reaction_pt": term,
                "report_count": cnt,
                "total_drug_reports": total_reports,
                "report_share_pct": share,
                "source_tag": "cloud_openfda_faers",
            })
        time.sleep(0.35)

    return signals, f"{len(signals)} FAERS signals across {len(FAERS_PAIN_DRUGS)} pain drugs"


def _query_rcsb_pdb_count(gene_symbol):
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.rcsb_gene_name.value",
                        "operator": "exact_match",
                        "value": gene_symbol,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_entity_source_organism.ncbi_taxonomy_id",
                        "operator": "exact_match",
                        "value": 9606,
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {"return_counts": True},
    }
    try:
        res = http_request_json("https://search.rcsb.org/rcsbsearch/v2/query", method="POST", payload=query, timeout=20)
        return int((res or {}).get("total_count") or 0)
    except Exception:
        return 0


def collect_opentargets_and_pdb():
    gql_full = """
    query TargetPainInfo($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        approvedName
        tractability {
          label
          modality
          value
        }
        safetyLiabilities {
          event
        }
        associatedDiseases(page: {index: 0, size: 15}) {
          rows {
            score
            disease {
              id
              name
            }
          }
        }
        knownDrugs(size: 100) {
          count
          rows {
            drugId
            prefName
            drugType
            phase
            status
            mechanismOfAction
            disease {
              name
            }
          }
        }
      }
    }
    """
    gql_fallback = """
    query TargetPainBasic($ensemblId: String!) {
      target(ensemblId: $ensemblId) {
        id
        approvedSymbol
        approvedName
        tractability {
          label
          modality
          value
        }
        associatedDiseases(page: {index: 0, size: 15}) {
          rows {
            score
            disease {
              id
              name
            }
          }
        }
        knownDrugs(size: 80) {
          count
          rows {
            drugId
            prefName
            drugType
            phase
            status
            mechanismOfAction
          }
        }
      }
    }
    """

    targets_out = []
    discovered_chembl_ids = set()

    for item in PAIN_TARGETS:
        ens_id = item["ensembl_id"]
        gene = item["gene_symbol"]
        t_group = item["target_group"]

        resp = None
        try:
            resp = http_request_json(
                "https://api.platform.opentargets.org/api/v4/graphql",
                method="POST",
                payload={"query": gql_full, "variables": {"ensemblId": ens_id}},
                timeout=25,
            )
            if not resp or resp.get("errors") or not (resp.get("data") or {}).get("target"):
                resp = http_request_json(
                    "https://api.platform.opentargets.org/api/v4/graphql",
                    method="POST",
                    payload={"query": gql_fallback, "variables": {"ensemblId": ens_id}},
                    timeout=25,
                )
        except Exception:
            resp = None

        tgt = ((resp or {}).get("data") or {}).get("target") or {}
        tracts = tgt.get("tractability") or []
        sm_pos = [t.get("label") for t in tracts if t.get("modality") == "SM" and t.get("value")]
        ab_pos = [t.get("label") for t in tracts if t.get("modality") == "AB" and t.get("value")]

        safeties = sorted({(s.get("event") or "").strip() for s in (tgt.get("safetyLiabilities") or []) if s.get("event")})
        assoc_rows = (tgt.get("associatedDiseases") or {}).get("rows") or []
        pain_scores = [
            float(r.get("score") or 0)
            for r in assoc_rows
            if re.search(r"pain|migraine|neuropath|analges|headache|pruritus|osteoarthritis", (r.get("disease") or {}).get("name") or "", re.I)
        ]
        top_dis = [
            f"{(r.get('disease') or {}).get('name', '')} ({round(float(r.get('score') or 0), 2)})"
            for r in assoc_rows[:5]
            if (r.get("disease") or {}).get("name")
        ]

        kd = tgt.get("knownDrugs") or {}
        kd_rows = kd.get("rows") or []
        for r in kd_rows:
            did = (r.get("drugId") or "").strip().upper()
            if did.startswith("CHEMBL"):
                discovered_chembl_ids.add(did)

        pdb_cnt = _query_rcsb_pdb_count(gene)

        targets_out.append({
            "gene_symbol": tgt.get("approvedSymbol") or gene,
            "ensembl_id": ens_id,
            "target_name": tgt.get("approvedName") or gene,
            "target_group": t_group,
            "tractability_sm": ", ".join(sm_pos[:4]) if sm_pos else "无小分子成药标记",
            "tractability_ab": ", ".join(ab_pos[:4]) if ab_pos else "无抗体成药标记",
            "pdb_structure_count": pdb_cnt,
            "known_drugs_count": int(kd.get("count") or len(kd_rows)),
            "safety_liabilities": "; ".join(safeties[:6]) if safeties else "Open Targets 暂无已知严重警示",
            "pain_association_score": round(max(pain_scores), 3) if pain_scores else (round(float(assoc_rows[0].get("score") or 0), 3) if assoc_rows else 0.0),
            "top_diseases": "; ".join(top_dis),
            "source_tag": "cloud_opentargets_pdb",
        })
        time.sleep(0.2)

    return targets_out, sorted(discovered_chembl_ids), f"{len(targets_out)} targets (Open Targets + RCSB PDB), {len(discovered_chembl_ids)} ChEMBL IDs discovered"


def collect_chembl_molecules(discovered_ids=None):
    mols_by_id = {}

    def _ingest_mol(m):
        cid = (m.get("molecule_chembl_id") or "").strip().upper()
        pref = (m.get("pref_name") or "").strip()
        if not cid:
            return
        props = m.get("molecule_properties") or {}
        structs = m.get("molecule_structures") or {}
        syns = sorted({(s.get("molecule_synonym") or "").strip() for s in (m.get("molecule_synonyms") or []) if s.get("molecule_synonym")})
        raw_phase = m.get("max_phase")
        try:
            max_phase = int(float(raw_phase)) if raw_phase is not None else 0
        except Exception:
            max_phase = 0
        try:
            mw = float(props.get("full_mwt")) if props.get("full_mwt") is not None else 0.0
        except Exception:
            mw = 0.0
        try:
            alogp = float(props.get("alogp")) if props.get("alogp") is not None else 0.0
        except Exception:
            alogp = 0.0
        try:
            ro5 = int(props.get("num_ro5_violations")) if props.get("num_ro5_violations") is not None else -1
        except Exception:
            ro5 = -1

        mols_by_id[cid] = {
            "chembl_id": cid,
            "pref_name": pref or cid,
            "synonyms": "; ".join(syns[:8]),
            "molecule_type": m.get("molecule_type") or "",
            "max_phase": max_phase,
            "first_approval": int(m.get("first_approval") or 0),
            "mechanism_of_action": "",
            "target_pref_name": "",
            "target_group": classify_pain_target(f"{pref} {' '.join(syns[:8])}"),
            "molecular_weight": round(mw, 2),
            "alogp": round(alogp, 2),
            "ro5_violations": ro5,
            "smiles": (structs.get("canonical_smiles") or "")[:240],
            "source_tag": "cloud_chembl",
        }

    # 1. Fetch curated benchmark pain drugs by preferred name
    for i in range(0, len(CHEMBL_SEED_NAMES), 8):
        batch = CHEMBL_SEED_NAMES[i : i + 8]
        url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?" + urllib.parse.urlencode({
            "pref_name__in": ",".join(batch),
            "limit": "50",
            "format": "json",
        })
        data = http_request_json(url, timeout=25)
        for m in (data or {}).get("molecules") or []:
            _ingest_mol(m)
        time.sleep(0.2)

    # 2. Fetch clinical molecules discovered from Open Targets (cap at 60 to keep bundle crisp)
    extra_ids = [cid for cid in (discovered_ids or []) if cid not in mols_by_id][:60]
    for i in range(0, len(extra_ids), 20):
        batch = extra_ids[i : i + 20]
        url = "https://www.ebi.ac.uk/chembl/api/data/molecule.json?" + urllib.parse.urlencode({
            "molecule_chembl_id__in": ",".join(batch),
            "limit": "50",
            "format": "json",
        })
        data = http_request_json(url, timeout=25)
        for m in (data or {}).get("molecules") or []:
            _ingest_mol(m)
        time.sleep(0.2)

    # 3. Attach Mechanism of Action from ChEMBL /mechanism endpoint
    all_ids = list(mols_by_id.keys())
    for i in range(0, len(all_ids), 20):
        batch = all_ids[i : i + 20]
        url = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json?" + urllib.parse.urlencode({
            "molecule_chembl_id__in": ",".join(batch),
            "limit": "100",
            "format": "json",
        })
        data = http_request_json(url, timeout=25)
        for mech in (data or {}).get("mechanisms") or []:
            cid = (mech.get("molecule_chembl_id") or "").upper()
            if cid in mols_by_id:
                moa = (mech.get("mechanism_of_action") or "").strip()
                tname = (mech.get("target_chembl_id") or "").strip()
                if moa and moa not in mols_by_id[cid]["mechanism_of_action"]:
                    mols_by_id[cid]["mechanism_of_action"] = (
                        (mols_by_id[cid]["mechanism_of_action"] + "; " + moa).strip("; ")
                    )
                if tname and not mols_by_id[cid]["target_pref_name"]:
                    mols_by_id[cid]["target_pref_name"] = tname
                # Re-classify target group with MoA context
                if mols_by_id[cid]["target_group"] == "其他镇痛机制":
                    mols_by_id[cid]["target_group"] = classify_pain_target(
                        f"{mols_by_id[cid]['pref_name']} {mols_by_id[cid]['mechanism_of_action']}"
                    )
        time.sleep(0.2)

    return list(mols_by_id.values()), f"{len(mols_by_id)} ChEMBL molecules"


def collect_europepmc_literature():
    arts_by_id = {}
    for q in EPMC_QUERIES:
        url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({
            "query": q,
            "format": "json",
            "resultType": "core",
            "pageSize": "35",
            "sort": "P_PDATE_D desc",
        })
        data = http_request_json(url, timeout=25)
        for art in ((data or {}).get("resultList") or {}).get("result") or []:
            rec = normalize_epmc_article(art)
            if rec:
                arts_by_id[rec["pub_id"]] = rec
        time.sleep(0.25)

    return list(arts_by_id.values()), f"{len(arts_by_id)} Europe PMC publications"


# ==================== ORCHESTRATOR ====================

def run_collection(mode="full", lookback_days=14, out_dir=DEFAULT_OUT_DIR):
    os.makedirs(out_dir, exist_ok=True)
    started = utc_now_iso()
    bundle = {
        "bundle_version": BUNDLE_VERSION,
        "generated_at": started,
        "generated_date": utc_today_str(),
        "mode": mode,
        "meta": {"sources": {}},
        "clinical_trials": [],
        "fda_approvals": [],
        "safety_signals": [],
        "target_science": [],
        "molecule_profiles": [],
        "literature_items": [],
    }

    def _safe_run(key, fn):
        t0 = time.time()
        try:
            res = fn()
            elapsed = round(time.time() - t0, 1)
            return res, {"status": "ok", "elapsed_s": elapsed}
        except Exception as err:
            elapsed = round(time.time() - t0, 1)
            print(f"[WARN] {key} failed: {err}", file=sys.stderr)
            return None, {"status": "error", "error": str(err), "elapsed_s": elapsed}

    # 1. ClinicalTrials.gov API v2
    r, st = _safe_run("ctgov", lambda: collect_ctgov(mode=mode, lookback_days=lookback_days))
    if r:
        bundle["clinical_trials"], st["summary"] = r
    bundle["meta"]["sources"]["ctgov"] = st

    # 2. openFDA Drugs@FDA
    r, st = _safe_run("openfda_approvals", collect_openfda_approvals)
    if r:
        bundle["fda_approvals"], st["summary"] = r
    bundle["meta"]["sources"]["openfda_approvals"] = st

    # 3. openFDA FAERS
    r, st = _safe_run("openfda_faers", collect_openfda_faers)
    if r:
        bundle["safety_signals"], st["summary"] = r
    bundle["meta"]["sources"]["openfda_faers"] = st

    # 4. Open Targets + RCSB PDB
    discovered_ids = []
    r, st = _safe_run("opentargets_pdb", collect_opentargets_and_pdb)
    if r:
        bundle["target_science"], discovered_ids, st["summary"] = r
    bundle["meta"]["sources"]["opentargets_pdb"] = st

    # 5. ChEMBL
    r, st = _safe_run("chembl", lambda: collect_chembl_molecules(discovered_ids))
    if r:
        bundle["molecule_profiles"], st["summary"] = r
    bundle["meta"]["sources"]["chembl"] = st

    # 6. Europe PMC
    r, st = _safe_run("europepmc", collect_europepmc_literature)
    if r:
        bundle["literature_items"], st["summary"] = r
    bundle["meta"]["sources"]["europepmc"] = st

    bundle["finished_at"] = utc_now_iso()
    bundle["counts"] = {
        "clinical_trials": len(bundle["clinical_trials"]),
        "fda_approvals": len(bundle["fda_approvals"]),
        "safety_signals": len(bundle["safety_signals"]),
        "target_science": len(bundle["target_science"]),
        "molecule_profiles": len(bundle["molecule_profiles"]),
        "literature_items": len(bundle["literature_items"]),
    }

    json_path = os.path.join(out_dir, "pain_radar_bundle.json")
    gz_path = os.path.join(out_dir, "pain_radar_bundle.json.gz")
    summary_path = os.path.join(out_dir, "sync_summary.md")

    raw_bytes = json.dumps(bundle, ensure_ascii=False, indent=2).encode("utf-8")
    with open(json_path, "wb") as f:
        f.write(raw_bytes)
    with gzip.open(gz_path, "wb", compresslevel=6) as f:
        f.write(raw_bytes)

    summary_lines = [
        f"# Pain Radar OS Cloud Sync Summary ({bundle['generated_date']})",
        "",
        f"- **Mode**: `{mode}`",
        f"- **Started**: {started}",
        f"- **Finished**: {bundle['finished_at']}",
        f"- **Bundle Size**: {round(len(raw_bytes) / 1024, 1)} KB (plain) / {round(os.path.getsize(gz_path) / 1024, 1)} KB (gzip)",
        "",
        "## Dimensions Collected",
        f"- **ClinicalTrials.gov**: {bundle['counts']['clinical_trials']} trials",
        f"- **openFDA Approvals**: {bundle['counts']['fda_approvals']} approvals",
        f"- **openFDA FAERS Signals**: {bundle['counts']['safety_signals']} AE pain-point signals",
        f"- **Target Science (Open Targets + RCSB PDB)**: {bundle['counts']['target_science']} targets",
        f"- **Molecule Profiles (ChEMBL)**: {bundle['counts']['molecule_profiles']} molecules",
        f"- **Literature (Europe PMC)**: {bundle['counts']['literature_items']} publications",
        "",
        "## Source Health",
    ]
    for k, info in bundle["meta"]["sources"].items():
        summary_lines.append(
            f"- `{k}`: **{info['status']}** ({info.get('elapsed_s', 0)}s) — {info.get('summary') or info.get('error', '')}"
        )

    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines) + "\n")

    print("\n".join(summary_lines))
    return bundle


def run_self_test():
    sample_study = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT05558410", "briefTitle": "Study of Suzetrigine (VX-548) After Abdominoplasty"},
            "statusModule": {
                "overallStatus": "COMPLETED",
                "startDateStruct": {"date": "2022-09"},
                "primaryCompletionDateStruct": {"date": "2023-12"},
                "lastUpdatePostDateStruct": {"date": "2025-02-01"},
            },
            "designModule": {"phases": ["PHASE3"], "enrollmentInfo": {"count": 1118}},
            "conditionsModule": {"conditions": ["Acute Pain", "Postoperative Pain"]},
            "armsInterventionsModule": {"interventions": [{"name": "VX-548", "type": "DRUG"}]},
            "sponsorCollaboratorsModule": {"leadSponsor": {"name": "Vertex Pharmaceuticals"}},
            "contactsLocationsModule": {"locations": [{"country": "United States"}]},
            "outcomesModule": {"primaryOutcomes": [{"measure": "SPID48 vs Placebo"}]},
        },
        "hasResults": True,
    }
    rec = normalize_ctgov_study(sample_study)
    assert rec and rec["trial_id"] == "NCT05558410" and rec["has_results"] == 1

    noise_study = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT00000001", "briefTitle": "Acupuncture and Yoga for Chronic Back Pain"},
            "statusModule": {"overallStatus": "RECRUITING"},
            "conditionsModule": {"conditions": ["Low Back Pain"]},
            "armsInterventionsModule": {"interventions": [{"name": "Acupuncture", "type": "OTHER"}]},
        }
    }
    assert normalize_ctgov_study(noise_study) is None

    sample_fda = {
        "application_number": "NDA219209",
        "sponsor_name": "VERTEX PHARMACEUTICALS INC",
        "products": [{"brand_name": "JOURNAVX", "dosage_form": "TABLET", "route": "ORAL", "active_ingredients": [{"name": "SUZETRIGINE", "strength": "50MG"}]}],
        "submissions": [{"submission_type": "ORIG", "submission_status": "AP", "submission_status_date": "20250130"}],
    }
    f_rec = normalize_openfda_approval(sample_fda)
    assert f_rec and f_rec["approval_no"] == "FDA-NDA219209" and f_rec["approval_date"] == "2025-01-30"

    sample_art = {
        "pmid": "39999999",
        "title": "NaV1.8 Inhibition With Suzetrigine in Acute Pain (NCT05558410).",
        "journalTitle": "N Engl J Med",
        "pubYear": "2025",
        "firstPublicationDate": "2025-04-01",
        "authorString": "Jones J, et al.",
        "citedByCount": 42,
        "doi": "10.1056/NEJMoa000001",
        "abstractText": "<p>Evaluated in NCT05558410 and NCT05553366.</p>",
    }
    a_rec = normalize_epmc_article(sample_art)
    assert a_rec and "NCT05558410" in a_rec["nct_ids"] and a_rec["pub_id"] == "PMID:39999999"

    print("SELF-TEST PASSED: All normalizers and noise filters verified.")


def main():
    parser = argparse.ArgumentParser(description="Pain Radar OS Cloud Multi-Dimensional Collector")
    parser.add_argument("--mode", choices=["full", "incremental"], default="full")
    parser.add_argument("--lookback-days", type=int, default=14)
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        run_self_test()
        return

    run_collection(mode=args.mode, lookback_days=args.lookback_days, out_dir=args.out_dir)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
fetch_hprc_index.py — Select the 4 requested HPRC assemblies by exact
assembly_name from the official Release 2 index, and write a compact
manifest with the numeric haplotype column preserved.

The official index uses NUMERIC haplotype values (1, 2).
Do NOT key on "maternal"/"paternal" strings — the official CSV does not
use them.  Derive human-readable labels from the assembly_name pattern.

Official columns: sample_id, haplotype, phasing, assembly_method,
assembly_method_version, assembly_date, assembly_name, source,
genbank_accession, assembly_md5, assembly_fai, assembly_gzi, assembly

Network hardening (Khoi):
  - Retries with exponential backoff + full jitter on transient failures.
  - Falls back to mirror URLs if the primary raw.githubusercontent URL is
    unreachable (proxy / rate limit / DNS).
  - --verify-urls does a cheap HEAD against each selected assembly's S3
    HTTPS mirror so a moved or renamed object is caught BEFORE we burn an
    hour on downloads.
  - --cache reuses a previously downloaded index (offline, or DNAnexus
    workers with no egress).

Usage:
    python3 scripts/fetch_hprc_index.py
    python3 scripts/fetch_hprc_index.py --verify-urls
    python3 scripts/fetch_hprc_index.py --cache work/manifests/index_raw.csv
    # Writes work/manifests/hprc_selected.csv
    # Exits non-zero if any of the 4 requested assemblies are missing.
"""
import argparse
import csv
import json
import os
import random
import sys
import time
import urllib.error
import urllib.request

# Primary and fallback locations for the official Release 2 index.
HPRC_INDEX_URLS = [
    ("raw.githubusercontent (refs/heads/main)",
     "https://raw.githubusercontent.com/human-pangenomics/"
     "hprc_intermediate_assembly/refs/heads/main/data_tables/"
     "assemblies_release2_v1.0.index.csv"),
    ("raw.githubusercontent (main)",
     "https://raw.githubusercontent.com/human-pangenomics/"
     "hprc_intermediate_assembly/main/data_tables/"
     "assemblies_release2_v1.0.index.csv"),
    ("github.com raw redirect",
     "https://github.com/human-pangenomics/hprc_intermediate_assembly/"
     "raw/main/data_tables/assemblies_release2_v1.0.index.csv"),
]

# Canonical assembly names from config/samples.yaml.
REQUESTED_ASSEMBLY_NAMES = [
    "HG00673_mat_hprc_r2_v1.0.1",
    "HG00673_pat_hprc_r2_v1.0.1",
    "HG00733_mat_hprc_r2_v1.0.1",
    "HG00733_pat_hprc_r2_v1.0.1",
]

MANIFEST_DIR = "work/manifests"
MANIFEST_PATH = os.path.join(MANIFEST_DIR, "hprc_selected.csv")
RAW_CACHE_PATH = os.path.join(MANIFEST_DIR, "index_raw.csv")
PROVENANCE_PATH = os.path.join(MANIFEST_DIR, "index_provenance.json")

USER_AGENT = "linear-scaling-pangenome-sv/0.1 (+BCM2026)"

# Transient HTTP codes worth retrying. 403 is included because
# raw.githubusercontent rate-limits with 403 rather than 429.
RETRYABLE_HTTP = {403, 408, 425, 429, 500, 502, 503, 504}


def _sleep_backoff(attempt, base=1.5, cap=20.0):
    """Exponential backoff with full jitter."""
    delay = min(cap, base * (2 ** attempt))
    time.sleep(random.uniform(0, delay))


def _http_get(url, timeout=60):
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode("utf-8")


def _http_head_ok(url, timeout=30):
    """HEAD a URL. Returns (ok, detail). Used only for --verify-urls."""
    req = urllib.request.Request(url, method="HEAD",
                                 headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            size = resp.headers.get("Content-Length", "?")
            return True, "HTTP %s, %s bytes" % (resp.status, size)
    except urllib.error.HTTPError as e:
        return False, "HTTP %s" % e.code
    except Exception as e:
        return False, str(e)


def s3_to_https(s3_uri):
    """Convert s3://human-pangenomics/... to the public HTTPS mirror."""
    if not s3_uri.startswith("s3://"):
        return s3_uri
    without_scheme = s3_uri[len("s3://"):]
    bucket, _, key = without_scheme.partition("/")
    return "https://%s.s3.amazonaws.com/%s" % (bucket, key)


def fetch_index(urls=HPRC_INDEX_URLS, retries=4, cache_path=None):
    """Download and parse the official HPRC Release 2 index CSV.

    Tries each candidate URL in order; each URL gets `retries` attempts
    with exponential backoff. Returns (rows, source_description).
    """
    if cache_path:
        if not os.path.exists(cache_path):
            print("FATAL: --cache given but not found: %s" % cache_path,
                  file=sys.stderr)
            sys.exit(1)
        print("Using cached index: %s" % cache_path)
        with open(cache_path, encoding="utf-8") as f:
            content = f.read()
        rows = list(csv.DictReader(content.splitlines()))
        print("  Loaded %d records from cache" % len(rows))
        return rows, "cache:%s" % cache_path

    last_error = None
    for label, url in urls:
        print("Fetching HPRC Release 2 index from %s:\n  %s" % (label, url))
        for attempt in range(retries):
            try:
                content = _http_get(url)
                rows = list(csv.DictReader(content.splitlines()))
                if not rows:
                    raise ValueError("index parsed to zero rows")

                # Cache raw bytes so a later run can go offline.
                os.makedirs(MANIFEST_DIR, exist_ok=True)
                with open(RAW_CACHE_PATH, "w", encoding="utf-8") as f:
                    f.write(content)

                print("  Downloaded %d records from index" % len(rows))
                print("  Raw index cached: %s" % RAW_CACHE_PATH)
                return rows, url
            except urllib.error.HTTPError as e:
                last_error = "HTTP %s from %s" % (e.code, url)
                if e.code not in RETRYABLE_HTTP:
                    print("  Non-retryable %s" % last_error, file=sys.stderr)
                    break
                print("  Attempt %d/%d failed (%s); backing off..."
                      % (attempt + 1, retries, last_error), file=sys.stderr)
            except Exception as e:
                last_error = "%s: %s" % (type(e).__name__, e)
                print("  Attempt %d/%d failed (%s); backing off..."
                      % (attempt + 1, retries, last_error), file=sys.stderr)

            if attempt < retries - 1:
                _sleep_backoff(attempt)

    print("FATAL: Could not fetch the HPRC index from any known URL.",
          file=sys.stderr)
    print("  Last error: %s" % last_error, file=sys.stderr)
    print("  If you have a previous run, retry with:\n"
          "    python3 scripts/fetch_hprc_index.py --cache %s" % RAW_CACHE_PATH,
          file=sys.stderr)
    sys.exit(1)


def _haplotype_label(assembly_name):
    """Derive human-readable label from the canonical assembly name."""
    if "_mat_" in assembly_name:
        return "maternal"
    if "_pat_" in assembly_name:
        return "paternal"
    return "unknown"


def select_assemblies(rows, requested):
    """Exact-match the requested assembly_names; fail hard on any miss."""
    lookup = {}
    for r in rows:
        aname = (r.get("assembly_name") or "").strip()
        if not aname:
            continue
        if aname in lookup:
            print("FATAL: Duplicate assembly_name in index: %s" % aname,
                  file=sys.stderr)
            sys.exit(1)
        lookup[aname] = r

    selected, missing = [], []
    for aname in requested:
        r = lookup.get(aname)
        if r is None:
            missing.append(aname)
        else:
            selected.append(r)

    if missing:
        print("FATAL: The following requested assemblies were NOT found "
              "in the official HPRC Release 2 index:", file=sys.stderr)
        for a in missing:
            print("  %s" % a, file=sys.stderr)
        # Nudge toward the right name if this looks like a version bump
        # (v1.0.1 -> v1.0.2) rather than a real removal.
        for a in missing:
            stem = a.rsplit("_v", 1)[0]
            near = sorted(k for k in lookup if k.startswith(stem))
            if near:
                print("  Similar names present for '%s': %s"
                      % (stem, ", ".join(near[:5])), file=sys.stderr)
        sys.exit(1)

    if len(selected) != len(requested):
        print("FATAL: Expected %d records, got %d. Possible duplicate."
              % (len(requested), len(selected)), file=sys.stderr)
        sys.exit(1)

    return selected


def verify_urls(selected):
    """HEAD each assembly's HTTPS mirror. Returns True if all reachable."""
    print("\nVerifying assembly URLs (HEAD requests)...")
    all_ok = True
    for r in selected:
        https = s3_to_https(r.get("assembly", ""))
        ok, detail = _http_head_ok(https)
        status = "OK  " if ok else "FAIL"
        print("  [%s] %s: %s" % (status, r.get("assembly_name", "?"), detail))
        if not ok:
            print("         %s" % https)
            all_ok = False
    if not all_ok:
        print("\nWARNING: at least one assembly URL was unreachable.\n"
              "  The S3 bucket is public; try the AWS CLI instead:\n"
              "    aws s3 --no-sign-request ls "
              "s3://human-pangenomics/working/HPRC/HG00673/assemblies/release2/",
              file=sys.stderr)
    return all_ok


def write_manifest(selected, out_path):
    parent = os.path.dirname(out_path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    fieldnames = [
        "sample_id", "haplotype", "haplotype_label", "assembly_name",
        "assembly_md5", "assembly_fai", "assembly_gzi", "assembly",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in selected:
            aname = r.get("assembly_name", "")
            w.writerow({
                "sample_id": r.get("sample_id", ""),
                "haplotype": r.get("haplotype", ""),
                "haplotype_label": _haplotype_label(aname),
                "assembly_name": aname,
                "assembly_md5": r.get("assembly_md5", ""),
                "assembly_fai": r.get("assembly_fai", ""),
                "assembly_gzi": r.get("assembly_gzi", ""),
                "assembly": r.get("assembly", ""),
            })


def write_provenance(source, n_rows, selected, urls_verified):
    """Record where the index came from — needed for the methods section."""
    os.makedirs(MANIFEST_DIR, exist_ok=True)
    payload = {
        "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "index_source": source,
        "index_record_count": n_rows,
        "requested": REQUESTED_ASSEMBLY_NAMES,
        "selected": [r.get("assembly_name", "") for r in selected],
        "urls_verified": urls_verified,
    }
    with open(PROVENANCE_PATH, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def main():
    ap = argparse.ArgumentParser(
        description="Select HPRC Release 2 assemblies from the official index.")
    ap.add_argument("--verify-urls", action="store_true",
                    help="HEAD each selected assembly URL before downloading.")
    ap.add_argument("--cache", metavar="CSV", default=None,
                    help="Use a cached raw index instead of the network "
                         "(e.g. %s)." % RAW_CACHE_PATH)
    ap.add_argument("--retries", type=int, default=4,
                    help="Attempts per URL before falling back (default 4).")
    ap.add_argument("--out", default=MANIFEST_PATH,
                    help="Manifest output path (default %s)." % MANIFEST_PATH)
    args = ap.parse_args()

    rows, source = fetch_index(retries=args.retries, cache_path=args.cache)
    selected = select_assemblies(rows, REQUESTED_ASSEMBLY_NAMES)

    urls_ok = None
    if args.verify_urls:
        urls_ok = verify_urls(selected)

    write_manifest(selected, args.out)
    write_provenance(source, len(rows), selected, urls_ok)

    print("\nSelected %d/%d assemblies — all found."
          % (len(selected), len(REQUESTED_ASSEMBLY_NAMES)))
    print("Manifest:   %s" % args.out)
    print("Provenance: %s" % PROVENANCE_PATH)
    for r in selected:
        aname = r.get("assembly_name", "")
        print("  %s haplotype=%s (%s): %s"
              % (r.get("sample_id", ""), r.get("haplotype", ""),
                 _haplotype_label(aname), aname))
        print("    S3: %s" % r.get("assembly", "N/A"))

    if urls_ok is False:
        # Manifest is still written; surface the problem via exit code so
        # `make fetch-index` fails loudly in CI.
        sys.exit(2)

    print("\nNext step:")
    print("  python3 scripts/download_hprc.py")


if __name__ == "__main__":
    main()

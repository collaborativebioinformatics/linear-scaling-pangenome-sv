#!/usr/bin/env python3
"""
download_hprc.py — Download HPRC assemblies, preserve native .fa.gz,
validate integrity with gzip -t, and keep compressed for DNAnexus storage.
"""
import argparse, csv, gzip, os, shutil, subprocess, sys

def load_manifest(path="work/manifests/hprc_selected.csv"):
    if not os.path.exists(path):
        print(f"Manifest not found: {path}", file=sys.stderr); sys.exit(1)
    with open(path) as f: rows = list(csv.DictReader(f))
    if not rows: print(f"Manifest empty: {path}", file=sys.stderr); sys.exit(1)
    return rows

def check_gzip(path):
    try:
        subprocess.run(["gzip","-t",path], check=True, capture_output=True); return True
    except: return False

def decompress(gz, fa):
    with gzip.open(gz,"rb") as fin, open(fa,"wb") as fout: shutil.copyfileobj(fin,fout)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--execute", action="store_true"); p.add_argument("--decompress", action="store_true")
    p.add_argument("--dest", default="work/downloads"); a = p.parse_args()
    rows = load_manifest(); os.makedirs(a.dest, exist_ok=True)
    print(f"Found {len(rows)} assemblies\n")
    for r in rows:
        sm, hap, name = r.get("sample_id","?"), r.get("haplotype","?"), r.get("assembly_name","?")
        uri = r.get("assembly","").strip(); gz = os.path.join(a.dest, f"{name}.fa.gz")
        if not uri: print(f"  SKIP {name}"); continue
        if a.execute:
            if os.path.exists(gz):
                print(f"  EXISTS {name}.fa.gz ({os.path.getsize(gz)/1e6:.0f} MB)")
            else:
                print(f"  DOWNLOAD {name} ({sm} {hap})..."); sys.stdout.flush()
                try: subprocess.run(["aws","s3","cp","--no-sign-request",uri,gz], check=True)
                except subprocess.CalledProcessError: continue
                except FileNotFoundError: print("FAILED: aws CLI not found"); sys.exit(1)
            if os.path.exists(gz):
                if check_gzip(gz): print(f"    OK ({os.path.getsize(gz)/1e6:.0f} MB, gzip valid)")
                else: print(f"    gzip FAILED, removing"); os.remove(gz)
        elif a.decompress:
            fa = os.path.join(a.dest, f"{name}.fa")
            if not os.path.exists(gz): print(f"  SKIP decompress {name}"); continue
            if os.path.exists(fa): print(f"  EXISTS {name}.fa"); continue
            decompress(gz, fa); print(f"  DECOMPRESS -> {name}.fa ({os.path.getsize(fa)/1e6:.0f} MB)")
        else:
            print(f"  QUEUED {name}: aws s3 cp --no-sign-request {uri} {gz}")
    if not a.execute and not a.decompress:
        print(f"\nDownload: {sys.argv[0]} --execute\nDecompress: {sys.argv[0]} --decompress")
        print("Keep .fa.gz compressed for DNAnexus storage. Decompress only for local use.")

if __name__ == "__main__": main()
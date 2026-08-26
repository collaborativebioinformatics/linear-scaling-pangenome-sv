"""
chunk_manifest.py — Utilities for chunk manifest management.
TODO (Quang): Add chunk status tracking, progress reporting.
"""
import csv
import os


def load_manifest(path: str = "work/chunks/chunk_manifest.tsv") -> list:
    if not os.path.exists(path):
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def update_chunk_status(chunk_id: str, status: str,
                        path: str = "work/chunks/chunk_manifest.tsv") -> bool:
    chunks = load_manifest(path)
    for c in chunks:
        if c["chunk_id"] == chunk_id:
            c["status"] = status
            break
    if chunks:
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=chunks[0].keys(), delimiter="\t")
            w.writeheader()
            w.writerows(chunks)
    return True


def main():
    chunks = load_manifest()
    print(f"Chunks in manifest: {len(chunks)}")
    for c in chunks:
        print(f"  {c['chunk_id']}: {c['status']}")


if __name__ == "__main__":
    main()
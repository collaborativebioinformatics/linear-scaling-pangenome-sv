"""
extract_regions.py — Extract orthologous regions from assemblies.

Uses mapping from map_chromosome.py to extract correct intervals.
Never naively assumes GRCh38 coordinates apply to de novo assemblies.

TODO (Khoi): Implement region extraction with samtools faidx when assemblies available.
"""


def main():
    print("=== Region Extraction ===")
    print("REAL-DATA STEP: Requires mapped sequences.")
    print("Run: python3 pipeline/prepare/map_chromosome.py")


if __name__ == "__main__":
    main()
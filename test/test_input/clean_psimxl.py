"""Clean pSimXL test data to keep only 20 spectra per type.

Usage:
    python test/clean_psimxl.py

This script:
    1. Reads Simulation.csv
    2. Selects 20 entries from each type (0=regular, 1=mono, 2=loop, 3=xlink)
    3. Writes a new CSV with only those 80 entries
    4. Cleans the MGF to only contain those spectra
"""
import os

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_input', 'pSimXL_test')
MAX_PER_TYPE = 20


def main():
    csv_path = os.path.join(TEST_DIR, 'Simulation.csv')
    mgf_path = os.path.join(TEST_DIR, 'Simulation.mgf')

    # Group entries by type
    types = {0: [], 1: [], 2: [], 3: []}
    with open(csv_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(',')
            if len(parts) < 8:
                continue
            t = int(parts[5])
            if t in types:
                types[t].append(line)

    # Select up to MAX_PER_TYPE from each type
    selected = []
    type_names = {0: 'regular', 1: 'mono', 2: 'loop', 3: 'xlink'}
    for t in [0, 1, 2, 3]:
        entries = types[t][:MAX_PER_TYPE]
        selected.extend(entries)
        print(f"Type {t} ({type_names[t]}): selected {len(entries)}/{len(types[t])}")

    # Write new CSV
    bak_csv = csv_path + '.bak'
    os.rename(csv_path, bak_csv)
    with open(csv_path, 'w') as f:
        for line in selected:
            f.write(line + '\n')
    print(f"\nWrote {len(selected)} entries to {csv_path}")

    # Collect titles to keep
    keep_titles = set()
    for line in selected:
        title = line.split(',')[0]
        keep_titles.add(title)

    # Clean MGF
    clean_mgf(mgf_path, keep_titles)

    # Remove backup
    os.remove(bak_csv)
    print("CSV backup removed.")


def clean_mgf(mgf_path, keep_titles):
    """Rewrite MGF to only contain spectra whose title is in keep_titles."""
    keep_lower = {t.lower() for t in keep_titles}

    bak_path = mgf_path + '.bak'
    os.rename(mgf_path, bak_path)

    kept_count = 0
    total_count = 0
    current_block = []
    in_block = False
    current_title = None

    with open(bak_path, 'r', encoding='utf-8') as fin, \
         open(mgf_path, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith('BEGIN IONS'):
                in_block = True
                current_block = [line]
                current_title = None
                total_count += 1
            elif in_block:
                current_block.append(line)
                if line.startswith('TITLE='):
                    current_title = line[6:].strip()
                elif line.startswith('END IONS'):
                    if current_title and current_title.lower() in keep_lower:
                        fout.writelines(current_block)
                        kept_count += 1
                    in_block = False
                    current_block = []
                    current_title = None

    print(f"MGF: kept {kept_count}/{total_count} spectra")

    # Verify all expected spectra were found
    found_lower = set()
    with open(mgf_path, 'r') as f:
        for line in f:
            if line.startswith('TITLE='):
                found_lower.add(line[6:].strip().lower())

    missing = keep_lower - found_lower
    if missing:
        print(f"WARNING: {len(missing)} spectra NOT found in MGF:")
        for m in sorted(missing)[:10]:
            print(f"  - {m}")
        print(f"Backup kept at {bak_path}")
    else:
        os.remove(bak_path)
        print("MGF backup removed (all spectra found).")


if __name__ == '__main__':
    main()

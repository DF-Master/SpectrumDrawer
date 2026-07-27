"""Clean MGF files to only contain spectra referenced in .plabel files.

Usage:
    python test/clean_mgf.py

For each test directory under test/test_input/, this script:
    1. Collects all spectrum names from .plabel files
    2. Rewrites corresponding .mgf files to only keep matching spectra
    3. Creates .bak backups (auto-removed if all spectra are found)

Add new test directories to TEST_DIRS as needed.
"""
import os
import re
import sys

TEST_INPUT = os.path.join(os.path.dirname(__file__), 'test_input')

# Add test directories here to include them in MGF cleanup
TEST_DIRS = [
    'pLink_BDG_test',
    'pLink_BS3_test',
    'pLink_SDA_test',
    # 'pSimXL_test',  # Already small, usually no cleanup needed
]


def collect_label_spectra(label_dir):
    """Collect all spectrum names from .plabel files in a directory."""
    spectra = set()
    for fname in os.listdir(label_dir):
        if fname.endswith('.plabel'):
            fpath = os.path.join(label_dir, fname)
            with open(fpath, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('name='):
                        name = line[5:].strip()
                        spectra.add(name)
    return spectra


def clean_mgf(mgf_path, keep_titles):
    """Rewrite MGF to only contain spectra whose title (case-insensitive) is in keep_titles."""
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

    print(f'  {os.path.basename(mgf_path)}: kept {kept_count}/{total_count} spectra')

    # Verify all expected spectra were found
    found_lower = set()
    with open(mgf_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('TITLE='):
                found_lower.add(line[6:].strip().lower())

    missing = keep_lower - found_lower
    if missing:
        print(f'  WARNING: {len(missing)} spectra from label files NOT found in MGF:')
        for m in sorted(missing)[:10]:
            print(f'    - {m}')
        if len(missing) > 10:
            print(f'    ... and {len(missing) - 10} more')
        print(f'  Backup kept at {bak_path}')
    else:
        os.remove(bak_path)
        print(f'  Backup removed (all spectra found).')


def main():
    for test_dir in TEST_DIRS:
        label_dir = os.path.join(TEST_INPUT, test_dir)
        if not os.path.isdir(label_dir):
            print(f'[{test_dir}] Directory not found, skipping.')
            continue

        mgf_files = [f for f in os.listdir(label_dir) if f.endswith('.mgf')]
        if not mgf_files:
            print(f'[{test_dir}] No MGF file found, skipping.')
            continue

        print(f'\n[{test_dir}] Collecting spectra from .plabel files...')
        keep_spectra = collect_label_spectra(label_dir)
        print(f'  Found {len(keep_spectra)} unique spectrum names in label files.')

        for mgf_file in mgf_files:
            mgf_path = os.path.join(label_dir, mgf_file)
            print(f'  Cleaning {mgf_file}...')
            clean_mgf(mgf_path, keep_spectra)

    print('\nDone!')


if __name__ == '__main__':
    main()

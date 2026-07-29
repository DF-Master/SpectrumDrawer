"""Clean pFind test data to keep 40 diverse spectra.

Usage:
    python test/clean_pfind.py

This script:
    1. Reads pFind-Filtered.spectra
    2. Selects 40 diverse entries (different lengths, modifications)
    3. Filters entries whose File_Name matches titles in the target MGF
    4. Rewrites .spectra and .mgf to only keep those entries
"""
import os
import csv
from collections import defaultdict

TEST_DIR = os.path.join(os.path.dirname(__file__), 'test_input', 'pFind_SDA_test')
MAX_SPECTRA = 40


def select_diverse_spectra(entries, mgf_titles, max_count=40):
    """Select diverse spectra based on length and modifications.
    
    Strategy:
    1. Group by sequence length
    2. For each length, prioritize different modification types
    3. Ensure entries exist in MGF
    """
    # Filter entries that exist in MGF
    valid_entries = [e for e in entries if e['File_Name'] in mgf_titles]
    print(f'Valid entries (in MGF): {len(valid_entries)}')
    
    # Group by length
    by_length = defaultdict(list)
    for e in valid_entries:
        seq_len = len(e['Sequence'])
        by_length[seq_len].append(e)
    
    # Get unique modification types
    mod_types = set()
    for e in valid_entries:
        mods = e.get('Modification', '')
        if mods:
            for mod in mods.split(';'):
                if mod.strip() and ',' in mod:
                    mod_name = mod.split(',')[1]
                    mod_types.add(mod_name)
    
    print(f'Unique modification types: {len(mod_types)}')
    print(f'Length range: {min(by_length.keys())}-{max(by_length.keys())}')
    
    # Select diverse entries
    selected = []
    selected_titles = set()
    
    # Round 1: Select one from each length (if available)
    for length in sorted(by_length.keys()):
        if len(selected) >= max_count:
            break
        candidates = by_length[length]
        # Prefer entries with diverse modifications
        for e in candidates:
            if e['File_Name'] not in selected_titles:
                selected.append(e)
                selected_titles.add(e['File_Name'])
                break
    
    # Round 2: Fill remaining slots with diverse modifications
    if len(selected) < max_count:
        # Sort by modification diversity
        def mod_diversity(e):
            mods = e.get('Modification', '')
            if not mods:
                return 0
            return len([m for m in mods.split(';') if m.strip()])
        
        remaining = [e for e in valid_entries if e['File_Name'] not in selected_titles]
        remaining.sort(key=mod_diversity, reverse=True)
        
        for e in remaining:
            if len(selected) >= max_count:
                break
            if e['File_Name'] not in selected_titles:
                selected.append(e)
                selected_titles.add(e['File_Name'])
    
    print(f'Selected: {len(selected)} entries')
    return selected


def main():
    spectra_path = os.path.join(TEST_DIR, 'pFind-Filtered.spectra')
    mgf_path = os.path.join(TEST_DIR, '20260055_1_HCDFT.mgf')

    # Collect MGF titles
    mgf_titles = set()
    with open(mgf_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('TITLE='):
                mgf_titles.add(line[6:].strip())
    print(f'MGF: {len(mgf_titles)} spectra')

    # Read spectra file
    entries = []
    with open(spectra_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            entries.append(row)
    print(f'Spectra: {len(entries)} entries')

    # Select diverse entries
    selected = select_diverse_spectra(entries, mgf_titles, MAX_SPECTRA)
    
    # Write new spectra file
    bak_spectra = spectra_path + '.bak'
    os.rename(spectra_path, bak_spectra)
    
    with open(spectra_path, 'w', encoding='utf-8', newline='') as f:
        fieldnames = entries[0].keys()
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for entry in selected:
            writer.writerow(entry)
    
    print(f'Wrote {len(selected)} entries to spectra file')

    # Clean MGF
    keep_titles = {e['File_Name'] for e in selected}
    clean_mgf(mgf_path, keep_titles)

    # Remove backup
    os.remove(bak_spectra)
    print('Spectra backup removed.')


def clean_mgf(mgf_path, keep_titles):
    """Rewrite MGF to only contain spectra whose title is in keep_titles."""
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
                    if current_title and current_title in keep_titles:
                        fout.writelines(current_block)
                        kept_count += 1
                    in_block = False
                    current_block = []
                    current_title = None

    print(f'MGF: kept {kept_count}/{total_count} spectra')

    # Verify
    found = set()
    with open(mgf_path, 'r') as f:
        for line in f:
            if line.startswith('TITLE='):
                found.add(line[6:].strip())

    missing = keep_titles - found
    if missing:
        print(f'WARNING: {len(missing)} spectra NOT found in MGF')
        print(f'Backup kept at {bak_path}')
    else:
        os.remove(bak_path)
        print('MGF backup removed (all spectra found).')


if __name__ == '__main__':
    main()

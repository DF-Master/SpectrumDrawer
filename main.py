"""SpectrumDrawer — MS/MS spectrum visualization tool.

Usage:
    python main.py --mgf spectra.mgf --ident results.csv --parser psimxl
    python main.py --mgf spectra.mgf --ident results.csv --parser psimxl -c my_config.yaml
    # For pLink .plabel files, --linker and --types are auto-detected:
    python main.py --mgf spectra.mgf --ident results.plabel --parser plink
"""

import argparse
import os
import re
import sys

# Ensure the parent directory is on sys.path so the package is importable
_pkg_parent = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _pkg_parent not in sys.path:
    sys.path.insert(0, _pkg_parent)

from SpectrumDrawer.core import SpectrumDrawer


def _infer_types_from_plabel(ident_path: str):
    """Infer spectrum types from a .plabel filename.

    Examples
    --------
    'xxx.regular.plabel'          -> [0]
    'xxx.mono-linked.BS3.plabel'  -> [1]
    'xxx.loop-linked.BS3.plabel'  -> [2]
    'xxx.cross-linked.BS3.plabel' -> [3]
    """
    basename = os.path.basename(ident_path).lower()
    if '.regular.' in basename or basename.endswith('regular.plabel'):
        return [0]
    if '.mono-linked.' in basename:
        return [1]
    if '.loop-linked.' in basename:
        return [2]
    if '.cross-linked.' in basename:
        return [3]
    return None


def _infer_linker_from_plabel(ident_path: str):
    """Read the [xlink] section from a .plabel file to get the crosslinker.

    Returns the xlink name, or None if xlink=NULL or not found.
    """
    with open(ident_path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if line == '[xlink]':
                for next_line in f:
                    next_line = next_line.strip()
                    m = re.match(r'^xlink=(\S+)', next_line)
                    if m:
                        name = m.group(1).strip()
                        return None if name.upper() == 'NULL' else name
                break
    return None


def main():
    parser = argparse.ArgumentParser(
        description='SpectrumDrawer — Draw annotated MS/MS spectra',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py --mgf spectra.mgf --ident results.csv --parser psimxl
  python main.py --mgf spectra.mgf --ident results.csv --parser psimxl -o ./output/
  python main.py --mgf spectra.mgf --ident results.csv --parser psimxl -c my_config.yaml --tol 10.0
  python main.py --mgf spectra.mgf --ident results.csv --parser psimxl --types 0 1
        """,
    )

    # Required arguments
    parser.add_argument('--mgf', '--spectrum', dest='spectrum_path',
                        required=True,
                        help='Path to MGF spectrum file')
    parser.add_argument('--ident', '--identification', dest='ident_path',
                        required=True,
                        help='Path to identification results file (e.g. pSimXL CSV)')
    parser.add_argument('--parser', dest='parser_name', default='psimxl',
                        help='Identification parser name (default: psimxl)')

    # Optional arguments
    parser.add_argument('-o', '--out-dir', dest='out_dir', default='./output',
                        help='Output directory for PNG files (default: ./output)')
    parser.add_argument('-c', '--config', dest='config_path', default=None,
                        help='Path to custom YAML config file')
    parser.add_argument('--linker', dest='linker_name', default=None,
                        help='Crosslinker name for mono-link mass lookup')
    parser.add_argument('--types', dest='spec_types', nargs='+', type=int,
                        default=None,
                        help='Spectrum types to draw (0=regular, 1=mono). '
                             'Default: all.')
    parser.add_argument('--tol', dest='tol_ppm', type=float, default=None,
                        help='Override mass tolerance (ppm) from config')
    parser.add_argument('--max-charge', dest='max_charge', type=int, default=None,
                        help='Max charge state for b/y fragment ions (default: 2)')

    args = parser.parse_args()

    # Validate input files
    if not os.path.isfile(args.spectrum_path):
        print(f'Error: spectrum file not found: {args.spectrum_path}')
        sys.exit(1)
    if not os.path.isfile(args.ident_path):
        print(f'Error: identification file not found: {args.ident_path}')
        sys.exit(1)

    # Auto-detect linker and types from .plabel files
    ident_lower = args.ident_path.lower()
    is_plabel = ident_lower.endswith('.plabel')

    if is_plabel:
        if args.linker_name is None:
            detected = _infer_linker_from_plabel(args.ident_path)
            if detected:
                print(f'Auto-detected crosslinker from .plabel: {detected}')
                args.linker_name = detected
        if args.spec_types is None:
            detected = _infer_types_from_plabel(args.ident_path)
            if detected:
                type_names = {0: 'regular', 1: 'mono-link', 2: 'loop-link', 3: 'cross-link'}
                labels = [type_names.get(t, str(t)) for t in detected]
                print(f'Auto-detected type(s) from .plabel: {", ".join(labels)}')
                args.spec_types = detected

    drawer = SpectrumDrawer(config_path=args.config_path)

    # Apply CLI overrides
    cli_overrides = {}
    if args.tol_ppm is not None:
        cli_overrides['processing.tol_ppm'] = args.tol_ppm
    if args.max_charge is not None:
        cli_overrides['processing.max_charge'] = args.max_charge
    if cli_overrides:
        drawer.config.apply_cli_overrides(**cli_overrides)

    drawer.run(
        spectrum_path=args.spectrum_path,
        ident_path=args.ident_path,
        parser=args.parser_name,
        out_dir=args.out_dir,
        linker_name=args.linker_name,
        spec_types=args.spec_types,
    )


if __name__ == '__main__':
    main()

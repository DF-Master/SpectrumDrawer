"""CrosslinkerProps 管理脚本（CLI 调度程序）。

用户通过命令行给出一个或多个交联剂（名称 + SMILES），程序:
  - 计算 tPSA（Å²）与 cLogP 等理化性质，并给出水溶性 / 透膜性的经验预测
  - 与常见交联剂比较，绘制 tPSA–cLogP 散点图
  - 导出包含常见交联剂数据的汇总 CSV（默认输出到本目录下的 test/）

用法示例::

    # 预测 4 个新型交联剂（Py3.13 环境）
    py -3.13 crosslinker_props_cli.py ^
        --xlinker "BDG=CC1(N=N1)CNC(CCCC(NCC2(N=N2)C)=O)=O" ^
        --xlinker "EBDA=CC1(N=N1)CCC(NCCNC(CCC2(N=N2)C)=O)=O" ^
        --xlinker "DCD=CC1(N=N1)CCNC(CCCCCCC(NCCC2(N=N2)C)=O)=O" ^
        --xlinker "DCDS=CC1(N=N1)CCNC(CCCC(NCCC2(N=N2)C)=O)=O"

    # 自定义输出目录 / 配置文件，跳过绘图
    py -3.13 crosslinker_props_cli.py --xlinker "BDG=...=O" -o my_out -c my_config.yaml --no-plot

    # 重新生成参考数据表 xlinkers.csv（当内置列表更新后）
    py -3.13 crosslinker_props_cli.py --init-data
"""

import argparse
import os
import sys

import pandas as pd

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

from crosslinker_props import (
    ConfigLoader, build_reference_csv, load_reference, compute_descriptors,
    add_predictions, plot_comparison, save_csv, _OUTPUT_COLUMNS,
)


def _parse_xlinker_spec(spec):
    """解析 'NAME=SMILES' -> (name, smiles)。"""
    if '=' not in spec:
        raise argparse.ArgumentTypeError(
            '格式应为 NAME=SMILES，得到: %r' % spec)
    name, smiles = spec.split('=', 1)
    name = name.strip()
    smiles = smiles.strip()
    if not name:
        raise argparse.ArgumentTypeError('交联剂名称为空: %r' % spec)
    if not smiles:
        raise argparse.ArgumentTypeError('%s 的 SMILES 为空' % name)
    return name, smiles


def _build_arg_parser():
    parser = argparse.ArgumentParser(
        description='CrosslinkerProps：预测新交联剂的水溶性与透膜性，'
                    '并与常见交联剂比较（独立工具）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument('--xlinker', action='append', default=[],
                        metavar='NAME=SMILES',
                        help='交联剂名称与 SMILES（可多次使用，一次可给多个）')
    parser.add_argument('-o', '--out-dir', type=str, default=None,
                        help='输出目录（默认取配置 output.dir，即本目录下的 test/）')
    parser.add_argument('-c', '--config', type=str, default=None,
                        help='自定义 YAML 配置路径')
    parser.add_argument('--no-plot', action='store_true',
                        help='不绘制比较图，仅输出 CSV 与命令行结果')
    parser.add_argument('--legend', action='store_true',
                        help='在比较图内绘制图例（默认关闭）')
    parser.add_argument('--init-data', action='store_true',
                        help='由内置列表重新生成参考数据表 xlinkers.csv'
                             '（无需 --xlinker）')
    return parser


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()

    cfg = ConfigLoader(args.config)

    csv_name = cfg.get('data', 'csv', default='xlinkers.csv')
    csv_path = os.path.join(_SCRIPT_DIR, csv_name)

    # ── --init-data：重建参考数据表 ─────────────────────────────
    if args.init_data:
        df = build_reference_csv(csv_path)
        print('参考数据表已重建: %s（共 %d 条）' % (csv_path, len(df)))
        if not args.xlinker:
            return 0

    # ── 解析用户输入 ─────────────────────────────────────────────
    specs = [_parse_xlinker_spec(s) for s in args.xlinker]
    if not specs:
        parser.error('请通过 --xlinker "NAME=SMILES" 提供至少一个交联剂'
                     '（或用 --init-data 重建数据表）')

    # ── 计算用户交联剂描述符 ─────────────────────────────────────
    user_rows = []
    for name, smiles in specs:
        props = compute_descriptors(smiles)
        if props is None:
            print('  [ERROR] SMILES 无法解析: %s = %s' % (name, smiles))
            continue
        user_rows.append({'name': name, 'smiles': smiles, **props,
                          'is_common': 0, 'origin': 'user'})
    if not user_rows:
        parser.error('所有输入的 SMILES 均无法解析')

    # ── 汇总：参考数据 + 用户输入 ────────────────────────────────
    ref = load_reference(csv_path)
    if 'origin' not in ref.columns:
        ref['origin'] = 'reference'
    combined = add_predictions(
        pd.concat([ref, pd.DataFrame(user_rows)], ignore_index=True), cfg)
    user_df = combined[combined['origin'] == 'user']

    # ── 命令行输出 ───────────────────────────────────────────────
    show_cols = ['name', 'tpsa', 'clogp', 'mol_weight', 'hbd', 'hba',
                 'solubility_pred', 'permeability_pred']
    print()
    print('===== 本次输入的交联剂预测结果 =====')
    print(user_df[show_cols].rename(columns={
        'name': '名称', 'tpsa': 'tPSA(Å²)', 'clogp': 'cLogP',
        'mol_weight': 'MW', 'hbd': 'HBD', 'hba': 'HBA',
        'solubility_pred': '水溶性(预测)', 'permeability_pred': '透膜性(预测)',
    }).to_string(index=False))
    print()
    print('说明: tPSA 单位 Å²；cLogP 为 Wildman-Crippen logP（RDKit）。')
    print('      水溶性由 cLogP 推断，透膜性由 tPSA 推断（经验阈值，见配置）。')

    # ── 输出 ─────────────────────────────────────────────────────
    out_dir = args.out_dir or os.path.join(
        _SCRIPT_DIR, cfg.get('output', 'dir', default='test'))
    os.makedirs(out_dir, exist_ok=True)
    out_csv = os.path.join(out_dir, cfg.get('output', 'csv_name',
                                            default='xlinkers_output.csv'))
    save_csv(combined, out_csv)
    print('汇总 CSV（含 %d 条常见/已收录 + %d 条本次输入）: %s'
          % (len(ref), len(user_df), os.path.abspath(out_csv)))

    if not args.no_plot:
        out_png = os.path.join(out_dir, cfg.get('output', 'plot_name',
                                                default='xlinkers_comparison.png'))
        user_names = list(user_df['name'])
        show_legend = args.legend or bool(cfg.get('plot', 'legend', default=False))
        plot_comparison(combined, user_names, out_png, cfg, show_legend=show_legend)
        print('比较图: %s' % os.path.abspath(out_png))
    return 0


if __name__ == '__main__':
    sys.exit(main())

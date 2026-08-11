"""
special_ions_pipeline.py — special_ions 流水线统一入口
======================================================
集中启动脚本 A（scan_report_ions.py，MGF 报告离子扫描）、
脚本 B（pparse_quant.py，pParse 定量）、脚本 C（group_by_plink.py，pLink 分组）
与对比脚本（compare_mono_regular.py，mono/regular × 实验/对照 分层对比）。

用法示例:
    # 只跑脚本 A（全谱图扫描）
    python special_ions_pipeline.py --step a --mgf-dir raw --ions-file xxx.ini \
        --all-spectra --jobs 8 --out-a out_a

    # 只跑脚本 B（pParse 定量，可按 title 清单过滤扫描）
    python special_ions_pipeline.py --step b --pparse-dir raw --ions-file xxx.ini \
        --prefix 20260512_BDG_plus --titles-b out_c/BDG/mono_link.fdr100.titles.txt \
        --out-b out_b

    # 只跑脚本 C（pLink 分组 + FDR -> title 清单）
    python special_ions_pipeline.py --step c --plink-dir pLink --out-c out_c

    # 分组重跑 A / B：用脚本 C 输出的各交联剂×组 fdr100 title 集逐组扫描
    python special_ions_pipeline.py --step a --mgf-dir raw --ions-file xxx.ini \
        --titles-file out_c/BDG/mono_link.fdr100.titles.txt --out-a out_c_runs_fdr100/BDG_mono
    python special_ions_pipeline.py --step b --pparse-dir raw --ions-file xxx.ini \
        --titles-b out_c/BDG/mono_link.fdr100.titles.txt --out-b out_b_runs_fdr100/BDG_mono
    # ...（4 交联剂 × {mono,regular} 共 8 次）

    # 分层对比：mode a 用脚本 A 结果（检出率/相对强度），mode b 用脚本 B 结果（+Ratio 分位）
    python special_ions_pipeline.py --step compare --compare-mode a \
        --runs-dir out_c_runs_fdr100 --out-c out_c --thresholds 0.10,0.50,1.00 \
        --out-compare out_c_runs_fdr100
    python special_ions_pipeline.py --step compare --compare-mode b \
        --runs-dir out_b_runs_fdr100 --out-c out_c --thresholds 0.10,0.50,1.00 \
        --out-compare out_b_runs_fdr100

    # 依次跑 A + B（all 模式，两脚本共用 --ions-file）
    python special_ions_pipeline.py --step all --mgf-dir raw --pparse-dir raw \
        --ions-file xxx.ini --all-spectra --out-a out_a --out-b out_b

依赖: Python 3.9+（仅标准库）
"""

import argparse
import os
import subprocess
import sys

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# scripts/special_ions -> scripts -> SpectrumDrawer -> sFFP_tools
_PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(_SCRIPT_DIR)))

_SCRIPT_A = 'scan_report_ions.py'
_SCRIPT_B = 'pparse_quant.py'
_SCRIPT_C = 'group_by_plink.py'
_SCRIPT_COMPARE = 'compare_mono_regular.py'

_DEFAULT_IONS = os.path.join(_PROJECT_ROOT, 'SpectrumDrawer', 'database',
                             'special_ions-jiangyida.ini')


def _run(script, args_list, step):
    """以子进程运行目标脚本，返回退出码。"""
    script_path = os.path.join(_SCRIPT_DIR, script)
    cmd = [sys.executable, '-B', script_path] + args_list
    print(f'\n[step {step}] {script}')
    print('  ' + ' '.join(cmd))
    return subprocess.run(cmd).returncode


def _build_a(args):
    """构造脚本 A 的命令行参数。"""
    argv = ['--special-ions-file', args.ions_file]
    if args.mgf_dir:
        argv += ['--mgf-dir', args.mgf_dir]
    if args.all_spectra:
        argv += ['--all']
    if args.titles_file:
        argv += ['--titles-file', args.titles_file]
    for flag, attr in (('--intensity-threshold', 'threshold'), ('--tol', 'tol'),
                       ('--top-n', 'top_n'), ('--jobs', 'jobs')):
        v = getattr(args, attr, None)
        if v is not None:
            argv += [flag, str(v)]
    if args.tol_mode:
        argv += ['--tol-mode', args.tol_mode]
    if args.ions:
        argv += ['--ions', args.ions]
    if args.files:
        argv += ['--files', args.files]
    argv += ['--out', args.out_a]
    return argv


def _build_b(args):
    """构造脚本 B 的命令行参数。"""
    argv = ['--ions-file', args.ions_file]
    if args.pparse_dir:
        argv += ['--pparse-dir', args.pparse_dir]
    if args.engine:
        argv += ['--engine', args.engine]
    if args.prefix:
        argv += ['--prefix', args.prefix]
    if args.titles_b:
        argv += ['--titles-file', args.titles_b]
    for flag, attr in (('--precursor-tol', 'precursor_tol'),
                       ('--ms1-window', 'ms1_window'),
                       ('--fragment-tol', 'fragment_tol')):
        v = getattr(args, attr, None)
        if v is not None:
            argv += [flag, str(v)]
    argv += ['--out', args.out_b]
    return argv


def _build_c(args):
    """构造脚本 C 的命令行参数。"""
    argv = ['--plink-dir', args.plink_dir]
    if args.agents:
        argv += ['--agents', args.agents]
    if args.thresholds:
        argv += ['--thresholds', args.thresholds]
    if args.score_column:
        argv += ['--score-column', args.score_column]
    if args.groups:
        argv += ['--groups', args.groups]
    if args.write_psm:
        argv += ['--write-psm']
    argv += ['--out', args.out_c]
    return argv


def _build_compare(args):
    """构造对比脚本的命令行参数。"""
    argv = ['--mode', args.compare_mode or 'a', '--runs-dir', args.runs_dir]
    if args.out_c:
        argv += ['--titles-dir', args.out_c]
    if args.thresholds:
        argv += ['--thresholds', args.thresholds]
    if args.agents:
        argv += ['--agents', args.agents]
    if args.groups:
        argv += ['--groups', args.groups]
    if args.match_mode:
        argv += ['--match-mode', args.match_mode]
    if args.exp_suffix:
        argv += ['--exp-suffix', args.exp_suffix]
    if args.ctrl_files:
        argv += ['--ctrl-files', args.ctrl_files]
    if args.agent_ions:
        argv += ['--agent-ions', args.agent_ions]
    argv += ['--ions-file', args.ions_file]
    if args.out_compare:
        argv += ['--out', args.out_compare]
    return argv


def main():
    ap = argparse.ArgumentParser(
        description='special_ions 流水线：统一启动脚本 A/B/C/对比')
    ap.add_argument('--step', required=True,
                    choices=['a', 'b', 'c', 'all', 'compare'],
                    help='a=脚本A(MGF扫描) b=脚本B(pParse定量) c=脚本C(pLink分组) '
                         'all=A+B compare=mono/regular 分层对比')
    # 公共
    ap.add_argument('--ions-file', default=_DEFAULT_IONS,
                    help=f'特殊离子 ini（A/B 共用；默认 {_DEFAULT_IONS}）')
    # 脚本 A 参数
    ap.add_argument('--mgf-dir', default=None, help='MGF 目录（step a/all 需要）')
    ap.add_argument('--all-spectra', action='store_true', help='A：处理目录内全部谱图')
    ap.add_argument('--titles-file', default=None, help='A：title 清单文件')
    ap.add_argument('--threshold', type=float, default=None, help='A：强度阈值')
    ap.add_argument('--tol-mode', default=None, choices=['ppm', 'da'], help='A：容差模式')
    ap.add_argument('--tol', type=float, default=None, help='A：da 容差值')
    ap.add_argument('--top-n', type=int, default=None, help='A：相关峰数量')
    ap.add_argument('--ions', default=None, help='A：只处理指定离子（逗号分隔）')
    ap.add_argument('--files', default=None, help='A：只处理指定文件前缀（逗号分隔）')
    ap.add_argument('--jobs', type=int, default=None, help='A：并行进程数')
    ap.add_argument('--out-a', default='out_a', help='A：输出目录')
    # 脚本 B 参数
    ap.add_argument('--pparse-dir', default=None, help='pParse 产物目录（step b/all 需要）')
    ap.add_argument('--engine', default=None, help='B：定量引擎')
    ap.add_argument('--prefix', default=None, help='B：只处理指定文件前缀（逗号分隔）')
    ap.add_argument('--titles-b', default=None,
                    help='B：title 清单文件（只处理清单内扫描，对应 B --titles-file）')
    ap.add_argument('--precursor-tol', type=float, default=None, help='B：前体容差 ppm')
    ap.add_argument('--ms1-window', type=int, default=None, help='B：MS1 搜索窗口')
    ap.add_argument('--fragment-tol', type=float, default=None, help='B：报告离子容差 ppm')
    ap.add_argument('--out-b', default='out_b', help='B：输出目录')
    # 脚本 C 参数
    ap.add_argument('--plink-dir', default='pLink', help='C：pLink 报告根目录')
    ap.add_argument('--agents', default=None, help='C：交联剂列表（逗号分隔）')
    ap.add_argument('--thresholds', default=None, help='C：FDR 阈值列表（逗号分隔）')
    ap.add_argument('--score-column', default=None, help='C：排序依据列（默认 Re-score）')
    ap.add_argument('--groups', default=None,
                    help='C/对比：分组（逗号分隔；C 需 Peptide_Type 名如 mono_link,regular；'
                         '对比可带简称如 "mono:mono_link,regular:regular"）')
    ap.add_argument('--write-psm', action='store_true', help='C：额外写出 PSM 明细')
    ap.add_argument('--out-c', default='out_c', help='C：输出目录')
    # 对比脚本参数（step compare）
    ap.add_argument('--compare-mode', default=None, choices=['a', 'b'],
                    help='对比：a=脚本A per_spectrum（检出率/相对强度） '
                         'b=脚本B per_scan（检出率+Ratio 分位）')
    ap.add_argument('--runs-dir', default=None,
                    help='对比：逐谱图结果根目录（含 {agent}_{分组简称}/ 子目录）')
    ap.add_argument('--out-compare', default=None, help='对比：输出目录')
    ap.add_argument('--agent-ions', default=None,
                    help='对比：各 agent 专属离子（CLI 内联'
                         '"agent:离子名=mz,离子名=mz;agent:..."，'
                         '如 "BDG:BDG154=154.086,BDG211=211.144"；缺省用 ini 前缀归属）')
    ap.add_argument('--match-mode', default=None, choices=['token', 'exact', 'contains'],
                    help='对比：文件前缀分层匹配方式（默认 token）')
    ap.add_argument('--exp-suffix', default=None,
                    help='对比：实验组文件后缀（必填，逗号分隔，与 agent 拼成 {agent}_后缀，如 plus,minus）')
    ap.add_argument('--ctrl-files', default=None,
                    help='对比：对照组文件前缀（必填，逗号分隔，如 blank_plus,blank_minus）')
    args = ap.parse_args()

    if not os.path.exists(args.ions_file):
        print(f'错误: --ions-file 不存在: {args.ions_file}')
        sys.exit(1)

    codes = []
    if args.step in ('a', 'all'):
        if not args.mgf_dir:
            print('错误: step a/all 需要 --mgf-dir')
            sys.exit(1)
        if not (args.all_spectra or args.titles_file):
            print('错误: 脚本 A 需要 --all-spectra 或 --titles-file 之一')
            sys.exit(1)
        codes.append(_run(_SCRIPT_A, _build_a(args), 'a'))
    if args.step in ('b', 'all'):
        if not args.pparse_dir:
            print('错误: step b/all 需要 --pparse-dir')
            sys.exit(1)
        codes.append(_run(_SCRIPT_B, _build_b(args), 'b'))
    if args.step == 'c':
        codes.append(_run(_SCRIPT_C, _build_c(args), 'c'))
    if args.step == 'compare':
        if not args.runs_dir:
            print('错误: step compare 需要 --runs-dir（脚本 A/B 逐谱图结果根目录）')
            sys.exit(1)
        codes.append(_run(_SCRIPT_COMPARE, _build_compare(args), 'compare'))

    if any(c != 0 for c in codes):
        print('\n流水线: 存在失败步骤')
        return 1
    print('\n流水线: 全部步骤完成')
    return 0


if __name__ == '__main__':
    sys.exit(main())

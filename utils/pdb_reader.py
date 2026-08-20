"""PDB 结构读取：解析各链氨基酸序列与 CA 原子坐标，计算交联位点间的 CA-CA 距离。

用于在 CSV 报告中为 cross-link 结果附加交联位点（Cα-Cα）的空间距离（单位 Å）。
肽段在链序列中进行严格（精确）匹配；一条肽段可能命中多个位置（如同源多聚体），
存在多个匹配时返回所有组合中的最小距离，找不到匹配或缺少坐标时返回 None。
"""

import math
import os

# 标准氨基酸三字母 → 单字母
_AA3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z',
}


class ProteinStructure:
    """蛋白质结构：链序列（SEQRES）+ 每链残基 CA 坐标（ATOM）。

    Attributes
    ----------
    chains : dict
        {chain_id: {'seq': str, 'ca': {pos: (x, y, z)}}}
        pos 为链序列中的 1-based 残基位置。
    """

    def __init__(self, pdb_path: str):
        self.pdb_path = pdb_path
        self.chains: dict = {}
        self._parse(pdb_path)

    # ── 解析 ──────────────────────────────────────────────────────

    def _parse(self, pdb_path: str):
        seqres = self._read_seqres(pdb_path)
        ca_atoms = self._read_ca_atoms(pdb_path)  # {chain: {resseq: (name, xyz)}}
        for chain, residues in ca_atoms.items():
            seq = seqres.get(chain, '')
            if not seq:
                # 无 SEQRES（AlphaFold/pymol 导出等）：按 ATOM 残基顺序重建链序列
                seq = ''.join(_AA3TO1.get(name.upper(), 'X')
                              for name, _ in residues.values())
            if not seq:
                continue
            ca = {}
            resseqs = sorted(residues.keys())
            if len(resseqs) == len(seq):
                for i, rs in enumerate(resseqs):
                    ca[i + 1] = residues[rs][1]
            else:
                for pos in range(1, len(seq) + 1):
                    if pos in residues:
                        ca[pos] = residues[pos][1]
            self.chains[chain] = {'seq': seq, 'ca': ca}

    @staticmethod
    def _read_seqres(pdb_path: str) -> dict:
        """读取 SEQRES 记录 → {chain_id: 单字母序列}。"""
        seqres = {}
        with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line.startswith('SEQRES') or len(line) < 19:
                    continue
                chain = line[11]
                residues = line[19:70].split()
                seq = ''.join(_AA3TO1.get(r.upper(), 'X') for r in residues)
                seqres[chain] = seqres.get(chain, '') + seq
        return seqres

    @staticmethod
    def _read_ca_atoms(pdb_path: str) -> dict:
        """读取 ATOM 记录中每链 CA 原子 → {chain: {resseq: (resname, (x, y, z))}}。

        多模型（NMR）时取第一个模型的坐标；按 resseq 排序即链残基顺序。
        """
        atoms = {}
        with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line.startswith('ATOM') or len(line) < 54:
                    continue
                if line[12:16].strip() != 'CA':
                    continue
                chain = line[21]
                resname = line[17:20].strip()
                try:
                    resseq = int(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                atoms.setdefault(chain, {}).setdefault(
                    resseq, (resname, (x, y, z)))
        return atoms

    # ── 距离计算 ──────────────────────────────────────────────────

    def get_ca_distance(self, alpha_seq: str, alpha_site: int,
                        beta_seq: str, beta_site: int):
        """返回 α/β 交联位点 CA-CA 距离（Å）。

        两条肽段在各链序列中严格匹配；所有 (α 命中, β 命中) 组合中取最小距离。
        任一肽段无匹配、位点无效或缺少 CA 坐标时返回 None。
        """
        alpha_hits = self._find_sites(alpha_seq, alpha_site)
        beta_hits = self._find_sites(beta_seq, beta_site)
        if not alpha_hits or not beta_hits:
            return None
        best = None
        for chain_a, pos_a in alpha_hits:
            coord_a = self.chains[chain_a]['ca'].get(pos_a)
            if coord_a is None:
                continue
            for chain_b, pos_b in beta_hits:
                coord_b = self.chains[chain_b]['ca'].get(pos_b)
                if coord_b is None:
                    continue
                # 欧氏距离（避免 math.dist 的 Python 3.8+ 版本要求）
                d = math.sqrt((coord_a[0] - coord_b[0]) ** 2 +
                              (coord_a[1] - coord_b[1]) ** 2 +
                              (coord_a[2] - coord_b[2]) ** 2)
                if best is None or d < best:
                    best = d
        return best

    def _find_sites(self, seq: str, site: int):
        """返回 [(chain_id, 链内 1-based 位置)] —— 肽段 seq 的第 site 个残基。"""
        if not seq or site <= 0 or site > len(seq):
            return []
        hits = []
        for chain, info in self.chains.items():
            chain_seq = info['seq']
            start = chain_seq.find(seq)
            while start != -1:
                hits.append((chain, start + site))
                start = chain_seq.find(seq, start + 1)
        return hits

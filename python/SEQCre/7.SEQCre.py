#!/usr/bin/env python
# coding: utf-8
"""
SEQCre.py

模块化重构后的 SeqCre 处理脚本，暴露 run_seqcre 接口
用法示例：
    from SEQCre import run_seqcre
    run_seqcre("输入目录", "输出目录", "subtype.xlsx")
"""
import os
import re
import sqlite3
import pandas as pd
import requests
from glob import glob
from time import sleep
from tqdm import tqdm
from urllib.parse import quote

# ====================
# 常量配置
# ====================
UNIPROT_API = "https://rest.uniprot.org/uniprotkb/search"
CACHE_DB    = "bio_cache.db"

# ====================
# 缓存初始化与查询函数
# ====================
def init_cache():
    """初始化本地缓存数据库"""
    conn = sqlite3.connect(CACHE_DB)
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS proteins
                 (name TEXT PRIMARY KEY, sequence TEXT)''')
    c.execute('''CREATE TABLE IF NOT EXISTS drugs
                 (name TEXT PRIMARY KEY, smiles TEXT)''')
    conn.commit()
    conn.close()


def get_protein_sequence(protein_name):
    """从 UniProt 获取蛋白质序列，含本地 SQLite 缓存"""
    conn = sqlite3.connect(CACHE_DB)
    try:
        # 1) 本地缓存查询
        cached = conn.execute(
            "SELECT sequence FROM proteins WHERE name=?", 
            (protein_name,)
        ).fetchone()
        if cached:
            return cached[0]
        # 2) 远程请求
        params = {
            "query": f"(gene:{protein_name}) AND (reviewed:true)",
            "format": "json",
            "fields": "sequence",
            "size": 1
        }
        resp = requests.get(UNIPROT_API, params=params)
        resp.raise_for_status()
        data = resp.json()
        if data.get("results"):
            seq = data["results"][0]["sequence"]["value"]
            conn.execute(
                "INSERT OR REPLACE INTO proteins VALUES (?, ?)",
                (protein_name, seq)
            )
            conn.commit()
            return seq
    except Exception as e:
        print(f"[Error] Protein query for {protein_name}: {e}")
    finally:
        conn.close()
    return None


def get_drug_smiles(drug_name):
    """从 PubChem 获取药物 SMILES，含本地 SQLite 缓存"""
    # 1) 只删除末尾的“_数字”
    clean = re.sub(r'_[0-9]+$', '', str(drug_name))
    # 2) 去掉末尾括号及其内容
    clean = re.sub(r'\s*\(.*\)\s*$', '', clean)
    # 3) 将中横线改为空格
    clean = clean.replace('-', ' ').strip()

    # —— 手动纠正规则 —— #
    # 4) BDP 前缀：去掉所有前导 0，并去掉空格，如 "BDP 00009066" / "BDP00009066" → "BDP9066"
    m = re.match(r'^(BDP)\s*0+(\d+)$', clean, flags=re.IGNORECASE)
    if m:
        clean = f"{m.group(1).upper()}{m.group(2)}"

    # 5) Picolinici acid 拼写改正
    if re.fullmatch(r'Picolinici acid', clean, flags=re.IGNORECASE):
        clean = "Picolinic acid"

    # 6) 跳过纯大写字母+数字（疑似蛋白名）
    if re.fullmatch(r'[A-Z0-9]+', clean):
        return None
    # —— 清洗完成 —— #

    # 7) URL 编码
    encoded = quote(clean)

    conn = sqlite3.connect(CACHE_DB)
    try:
        # 缓存查询
        cached = conn.execute(
            "SELECT smiles FROM drugs WHERE name=?", 
            (clean,)
        ).fetchone()
        if cached:
            return cached[0]

        # 网络请求
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{encoded}/property/CanonicalSMILES/txt"
        )
        resp = requests.get(url, verify=False)
        resp.raise_for_status()
        smiles = resp.text.strip()

        # 缓存写入
        conn.execute(
            "INSERT OR REPLACE INTO drugs VALUES (?, ?)",
            (clean, smiles)
        )
        conn.commit()
        return smiles

    except Exception as e:
        print(f"[Error] Drug query for {clean!r}: {e}")
        return None

    finally:
        conn.close()

# ====================
# 对外主流程函数
# ====================
def run_seqcre(input_base: str,
               output_base: str,
               subtype_file: str) -> None:
    """
    对外接口：扫描 input_base 下所有 DRN_info_*.txt，
    按照 subtype.xlsx 中的分型列表，
    获取蛋白序列与药物 SMILES，并写出 CSV。

    参数:
      input_base   - 原始 DRN_info 文件根目录
      output_base  - 结果输出根目录
      subtype_file - 包含 Sample, Subtype 的 Excel 文件路径
    返回:
      None
    """
    # 1) 初始化缓存
    init_cache()

    # 2) 读取分型
    subtype_df = pd.read_excel(subtype_file)
    required = {'Sample', 'Subtype'}
    if not required.issubset(subtype_df.columns):
        raise ValueError(f"subtype.xlsx 必须包含列: {required}")
    phenotypes = subtype_df['Subtype'].unique().tolist()

    # 3) 收集所有 DRN_info 文件
    pattern = os.path.join(input_base, "**", "DRN_info_*.txt")
    files = glob(pattern, recursive=True)
    if not files:
        print(f"未在 {input_base} 下找到任何 DRN_info 文件")
        return

    # 4) 逐文件处理
    for fp in tqdm(files, desc="Processing DRN_info files"):
        try:
            parts = fp.split(os.sep)
            if len(parts) < 5:
                continue
            phenotype     = parts[-5]
            module_dir    = parts[-4]
            network_method= parts[-3]
            module        = parts[-2]
            if phenotype not in phenotypes:
                continue

            df = pd.read_csv(fp, sep="\t")
            if not {'node','type'}.issubset(df.columns):
                continue

            # 蛋白部分
            prot = df[df['type']=='protein'].copy()
            prot['sequence'] = prot['node'].apply(
                lambda x: get_protein_sequence(x) or 'Not Found'
            )

            # 药物部分
            drug = df[df['type']=='drug'].copy()
            drug['node'] = drug['node'].apply(lambda x: re.sub(r'_\d+$','',str(x)))
            drug['SMILES'] = drug['node'].apply(
                lambda x: get_drug_smiles(x) or 'Not Found'
            )

            # 输出目录
            out_dir = os.path.join(
                output_base,
                phenotype,
                module_dir,
                network_method,
                module
            )
            os.makedirs(out_dir, exist_ok=True)

            # 文件名
            seq_file    = f"seq_file_{module_dir}_{network_method}_{phenotype}_{module}.csv"
            smiles_file = f"smiles_file_{module_dir}_{network_method}_{phenotype}_{module}.csv"

            if not prot.empty:
                prot[['node','sequence']].to_csv(
                    os.path.join(out_dir, seq_file),
                    index=False
                )
            if not drug.empty:
                drug[['node','SMILES']].to_csv(
                    os.path.join(out_dir, smiles_file),
                    index=False
                )

            sleep(1)
        except Exception as e:
            print(f"[Error] Processing {fp}: {e}")

# 指定公开接口
__all__ = ['run_seqcre', 'get_protein_sequence', 'get_drug_smiles']




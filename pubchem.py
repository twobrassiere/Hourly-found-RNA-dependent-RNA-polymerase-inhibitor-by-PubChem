import os
import re
import time
import json
import requests
import pandas as pd
import pubchempy as pcp
import gspread
from datetime import datetime, timezone
from google.oauth2.service_account import Credentials

# 病毒株清單
ORGANISM_LIST = [
    'Hepatitis B virus', 'Varicella zoster virus', 'Human cytomegalovirus',
    'Herpes simplex virus', 'Epstein–Barr virus', 'Variola virus',
    'Molluscum contagiosum virus', 'Mpox virus', 'Human papillomavirus',
    'SARS-CoV-2', 'MERS-CoV', 'SARS-CoV-1', 'Human coronavirus 229E',
    'Human coronavirus NL63', 'Hepatitis C virus', 'Dengue virus',
    'Yellow fever virus', 'Zika virus', 'West Nile virus',
    'Japanese encephalitis virus', 'Hepatitis A virus', 'Coxsackie virus',
    'Norovirus', 'Chikungunya virus', 'Rubella virus', 'Influenza virus',
    'Respiratory syncytial virus', 'Mumps virus', 'Measles virus',
    'Rabies virus', 'Ebola virus', 'Lassa virus',
    'Crimean–Congo hemorrhagic fever virus', 'Rotavirus', 'HIV'
]

# RdRp 蛋白質關鍵字
TARGET_KEYWORDS = [
    "RNA-directed RNA polymerase", "RNA-dependent RNA polymerase", "RNA dependent RNA polymerase",
    "RdRp", "Replicase", "RNA polymerase", "NSP12", "NS5B", "NS5", "L protein", "Large protein",
    "PB1", "VP63", "Replicase polyprotein 1ab", "Replicase polyprotein 1a", "Replicase polyprotein",
    "Genome polyprotein", "Polyprotein P1234", "Polyprotein", "nsp12", "Non-structural protein 12",
    "NS5B RNA-dependent RNA polymerase", "Non-structural protein 5", "Polymerase acidic protein",
    "Polymerase basic protein 1", "Polymerase basic protein 2", "Polymerase basic protein", "PB2",
    "PA protein", "RSV polymerase", "RNA polymerase L", "Large structural protein", "Protein L",
    "3D", "3Dpol", "3D polymerase", "Protein 3D", "Protein 3CD"
]

def clean_sheet_title(name: str) -> str:
    """清理工作表名稱，符合 Google Sheets 規範（不超過 31 字元）"""
    return re.sub(r'[\\/*?\[\]:]', '_', name)[:31].strip()

def get_next_organism(organism_list: list) -> str:
    """依據當前 UTC 時間的小時數自動輪播病毒株"""
    env_override = os.environ.get("TARGET_ORGANISM")
    if env_override:
        return env_override
    utc_now = datetime.now(timezone.utc)
    hour_epoch = int(utc_now.timestamp() // 3600)
    selected_index = hour_epoch % len(organism_list)
    return organism_list[selected_index]

def search_pubchem_assays(organism: str) -> list:
    """透過 NCBI Entrez / PubChem API 搜尋符合病毒株與 RdRp 關鍵字的 BioAssay AID 列表"""
    print(f"🔍 搜尋 PubChem BioAssays: {organism} + RdRp 關鍵字...")
    aids = set()
    
    # 核心關鍵字
    core_keywords = ["RdRp", "RNA-dependent RNA polymerase", "NS5B", "nsp12", "replicase", "L protein"]
    
    for kw in core_keywords:
        term = f'"{organism}"[Organism] AND "{kw}"[Assay Name]'
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        params = {
            "db": "pcassay",
            "term": term,
            "retmode": "json",
            "retmax": 50
        }
        try:
            res = requests.get(url, params=params, timeout=10)
            if res.status_code == 200:
                id_list = res.json().get("esearchresult", {}).get("idlist", [])
                aids.update(id_list)
            time.sleep(0.3)  # 遵守 API 限速規範
        except Exception as e:
            print(f"⚠️ 查詢 API 失敗 ({kw}): {e}")
            
    return list(aids)

def fetch_compounds_safe(cids: list) -> list:
    """安全地使用 PubChemPy 批量獲取 Compound 資訊（加入分批與限速控制）"""
    if not cids:
        return []
    
    compounds = []
    batch_size = 10  # 每批次最多查詢 10 個 CID，防止 API 超時或 429
    
    for i in range(0, len(cids), batch_size):
        batch = cids[i:i + batch_size]
        try:
            # 呼叫 PubChemPy 取得化合物物件
            batch_compounds = pcp.get_compounds(batch)
            compounds.extend(batch_compounds)
            time.sleep(0.4)  # 控制請求頻率低於 5 req/sec
        except Exception as e:
            print(f"⚠️ PubChemPy 提取批次 {batch} 失敗: {e}")
            
    return compounds

def get_pubchem_rdrp_data(organism: str) -> pd.DataFrame:
    """獲取特定病毒株之 RdRp 活性化合物與結構數據"""
    aids = search_pubchem_assays(organism)
    if not aids:
        print(f"📊 針對 [{organism}] 未搜尋到對應的 PubChem BioAssay。")
        return pd.DataFrame()

    print(f"🧪 找到 {len(aids)} 個 BioAssays，開始提取活性化合物數據...")
    records = []
    
    for aid in aids[:5]:  # 限制讀取前 5 個 Assay
        try:
            # 取得 Assay 內標示為 Active 的化合物 CIDs
            assay_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/{aid}/cids/JSON?assay_outcome=active"
            res = requests.get(assay_url, timeout=10)
            if res.status_code != 200:
                continue
            
            cids = res.json().get("IdentifierList", {}).get("CID", [])
            if not cids:
                continue

            print(f"  - Assay AID {aid}: 發現 {len(cids)} 個活性化合物 (CID)")

            # 使用安全批次模式呼叫 PubChemPy
            target_cids = cids[:20]  # 單一 Assay 採樣前 20 個 CID
            compounds = fetch_compounds_safe(target_cids)
            
            for comp in compounds:
                records.append({
                    "Organism": organism,
                    "PubChem AID": str(aid),
                    "PubChem CID": str(comp.cid),
                    "Compound Name": comp.iupac_name or getattr(comp, 'synonyms', ['N/A'])[0] if getattr(comp, 'synonyms', None) else "N/A",
                    "Molecular Formula": comp.molecular_formula or "N/A",
                    "Molecular Weight": comp.molecular_weight or "N/A",
                    "Canonical SMILES": comp.canonical_smiles or "N/A",
                    "InChIKey": comp.inchikey or "N/A"
                })
        except Exception as e:
            print(f"⚠️ 擷取 AID {aid} 失敗: {e}")

    df = pd.DataFrame(records)
    if not df.empty:
        # 依據 PubChem CID 與 AID 去重複
        df.drop_duplicates(subset=["PubChem AID", "PubChem CID"], inplace=True)
    return df

def update_google_sheet(df: pd.DataFrame, organism_name: str):
    """將結果同步至 Google Sheets 並進行去重與欄寬自動調整"""
    creds_json = os.environ.get("GCP_SERVICE_ACCOUNT_KEY")
    spreadsheet_id = os.environ.get("SPREADSHEET_ID")

    if not creds_json or not spreadsheet_id:
        print("⚠️ 缺乏 GCP Credentials 或 SPREADSHEET_ID，跳過 Sheet 同步。")
        return

    scopes = ["https://www.googleapis.com/auth/spreadsheets", "https://www.googleapis.com/auth/drive"]
    try:
        creds_dict = json.loads(creds_json)
        creds = Credentials.from_service_account_info(creds_dict, scopes=scopes)
    except Exception as e:
        print(f"❌ GCP Credentials 解析失敗: {e}")
        return

    gc = gspread.authorize(creds)
    sh = gc.open_by_key(spreadsheet_id)

    sheet_title = clean_sheet_title(organism_name)
    existing_worksheets = {ws.title: ws for ws in sh.worksheets()}

    default_headers = [
        "Organism", "PubChem AID", "PubChem CID", "Compound Name",
        "Molecular Formula", "Molecular Weight", "Canonical SMILES", "InChIKey"
    ]

    if sheet_title not in existing_worksheets:
        worksheet = sh.add_worksheet(title=sheet_title, rows="1000", cols="20")
        worksheet.append_row(default_headers)
        existing_keys = set()
    else:
        worksheet = existing_worksheets[sheet_title]
        all_values = worksheet.get_all_values()
        existing_keys = set()
        if len(all_values) > 1:
            header = all_values[0]
            aid_idx = header.index("PubChem AID") if "PubChem AID" in header else 1
            cid_idx = header.index("PubChem CID") if "PubChem CID" in header else 2
            for row in all_values[1:]:
                if len(row) > max(aid_idx, cid_idx):
                    existing_keys.add(f"{row[aid_idx]}_{row[cid_idx]}")

    if df.empty:
        print(f"ℹ️ [{sheet_title}]: 無新數據需要寫入。")
        return

    # 過濾已存在於 Sheet 中的記錄
    new_rows = []
    for _, row in df.iterrows():
        key = f"{row['PubChem AID']}_{row['PubChem CID']}"
        if key not in existing_keys:
            new_rows.append(row.fillna('N/A').tolist())

    if new_rows:
        worksheet.append_rows(new_rows)
        print(f"✨ 成功新增 {len(new_rows)} 筆 PubChem 數據至 [{sheet_title}]！")
        
        # 自動調整欄寬
        try:
            sh.batch_update({
                "requests": [{
                    "autoResizeDimensions": {
                        "dimensions": {
                            "sheetId": worksheet.id,
                            "dimension": "COLUMNS",
                            "startIndex": 0
                        }
                    }
                }]
            })
        except Exception as e:
            print(f"⚠️ 自動調整欄寬失敗: {e}")
    else:
        print(f"ℹ️ [{sheet_title}]: 數據皆已存在，無需更新。")

if __name__ == "__main__":
    current_organism = get_next_organism(ORGANISM_LIST)
    print(f"⏰ 本次輪播病毒株: '{current_organism}'")
    df_result = get_pubchem_rdrp_data(current_organism)
    if not df_result.empty:
        update_google_sheet(df_result, current_organism)
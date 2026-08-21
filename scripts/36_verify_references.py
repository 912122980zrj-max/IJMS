#!/usr/bin/env python
"""36_verify_references.py

逐条核验 PLOS 手稿的 28 条参考文献真实性（只读网络 + 手稿文件）：
  - 提取 References 章节中 "N. ..." 条目与 DOI；
  - 对每条 DOI 查询 Crossref，比对标题（词集合 Jaccard）、年份、期刊、
    卷/期/页码（或文章号）、第一作者姓氏；
  - 对无 DOI 的软件条目（R Core Team）做 URL 可达性检查；
  - 输出 PASS/FAIL 清单，不修改任何文件。
"""

from __future__ import annotations

import os
import re
import sys
import time
import unicodedata

import requests

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass


MD = r"E:/sheng xin/ObstructiveNephropathy_MRG/submission/plos/manuscript/manuscript_plos.md"

JOURNAL_MAP = {
    "annu rev physiol": "annual review of physiology",
    "kidney int": "kidney international",
    "nat rev mol cell biol": "nature reviews molecular cell biology",
    "jci insight": "jci insight",
    "j transl med": "journal of translational medicine",
    "plos one": "plos one",
    "sci rep": "scientific reports",
    "sci data": "scientific data",
    "nat commun": "nature communications",
    "nat biotechnol": "nature biotechnology",
    "bmc bioinformatics": "bmc bioinformatics",
    "nucleic acids res": "nucleic acids research",
    "j stat softw": "journal of statistical software",
    "mach learn": "machine learning",
    "cell syst": "cell systems",
    "proc natl acad sci u s a": (
        "proceedings of the national academy of sciences"),
    "cell rep": "cell reports",
    "nat protoc": "nature protocols",
}


def norm(s: str) -> str:
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode()
    s = re.sub(r"[^a-z0-9 ]", " ", s.lower())
    return " ".join(s.split())


def words(s: str):
    return set(w for w in norm(s).split() if len(w) > 2
               and w not in {"the", "and", "for", "with", "from", "into",
                             "using", "data", "analysis"})


def jaccard(a, b) -> float:
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)


def parse_refs(md_text: str):
    section = md_text.split("## References")[1].split("## Table 1")[0]
    refs = []
    for m in re.finditer(r"(?m)^(\d+)\.\s+(.*)$", section):
        refs.append((int(m.group(1)), m.group(2)))
    return refs


def check_crossref(body: str) -> tuple[bool, list[str]]:
    fails = []
    doi = re.search(r"doi:\s*(\S+)", body, re.I)
    if not doi:
        return False, ["no DOI found in entry"]
    doi = doi.group(1).rstrip(".")
    url = "https://api.crossref.org/works/" + doi
    for attempt in range(3):
        try:
            r = requests.get(url, timeout=40,
                             headers={"User-Agent": "ref-audit/1.0 (mailto:audit@example.org)"})
            if r.status_code == 200:
                msg = r.json()["message"]
                break
        except requests.RequestException:
            time.sleep(2 * (attempt + 1))
    else:
        return False, [f"Crossref request failed (status {getattr(r, 'status_code', '?')})"]

    # 年份（优先印刷版年份，其次在线年份）
    try:
        year = ((msg.get("published-print") or msg.get("issued") or {})
                .get("date-parts", [[None]])[0][0])
    except Exception:
        year = None
    m_year = re.search(r"\.\s*(\d{4})[;.]", body)
    if year and m_year and str(year) != m_year.group(1):
        fails.append(f"year manuscript={m_year.group(1)} vs crossref={year}")

    # 标题与期刊：年份前的文本按 ". " 分段，最长段=标题，末段=期刊
    pre = body[:m_year.start()] if m_year else body
    segs = [s.strip() for s in pre.split(". ") if s.strip()]
    has_volume = re.search(r"[.;]\s*\d+\s*[(:]", body) is not None
    m_title = (segs[-2] if has_volume and len(segs) >= 2
               else (max(segs, key=len) if segs else ""))
    journal = segs[-1] if segs else ""
    cx_title = (msg.get("title") or [""])[0]
    j = jaccard(words(m_title), words(cx_title))
    is_book = not has_volume
    subtitle_ok = is_book and norm(cx_title) in norm(m_title)
    if j < 0.55 and not subtitle_ok:
        fails.append(f"title similarity {j:.2f} (ms='{m_title[:50]}' cx='{cx_title[:50]}')")

    # 期刊
    cx_j = (msg.get("container-title") or [""])[0]
    if journal and has_volume:
        abbr = norm(journal)
        mapped = None
        for k, v in JOURNAL_MAP.items():
            if abbr.startswith(k) or k.startswith(abbr):
                mapped = v
                break
        if mapped:
            if norm(cx_j) != norm(mapped):
                fails.append(f"journal '{journal}' vs crossref '{cx_j}'")
        else:
            if norm(cx_j) != abbr and not norm(cx_j).startswith(abbr):
                fails.append(f"journal '{journal}' vs crossref '{cx_j}'")

    # 卷/期/页/文章号
    m_vol = re.search(r"[.;]\s*(\d+)\s*\((\d+)\)\s*:\s*([\d\-–—eE]+)", body)
    if m_vol:
        v, iss, pg = m_vol.groups()
        if str(msg.get("volume") or "") != v:
            fails.append(f"volume manuscript={v} vs crossref={msg.get('volume')}")
        if str(msg.get("issue") or "") != iss:
            fails.append(f"issue manuscript={iss} vs crossref={msg.get('issue')}")
        cx_pg = str(msg.get("page") or msg.get("article-number") or "")
        cx_pg_n = cx_pg.replace("–", "-").replace("—", "-")
        pg_n = pg.replace("–", "-").replace("—", "-")
        if cx_pg and cx_pg_n != pg_n and not cx_pg_n.startswith(pg_n + "-"):
            fails.append(f"pages manuscript={pg} vs crossref={cx_pg}")
    else:
        # 仅卷号（如 BMC Bioinformatics 14:7）
        m_vol2 = re.search(r"[.;]\s*(\d+)\s*:\s*([\d\-–—eE]+)", body)
        if m_vol2:
            v2, art = m_vol2.groups()
            if str(msg.get("volume") or "") != v2:
                fails.append(f"volume manuscript={v2} vs crossref={msg.get('volume')}")
            cx_art = str(msg.get("article-number") or msg.get("page") or "")
            cx_art_n = cx_art.replace("–", "-").replace("—", "-")
            art_n = art.replace("–", "-").replace("—", "-")
            if cx_art and cx_art_n != art_n and not cx_art_n.startswith(art_n + "-"):
                fails.append(f"article manuscript={art} vs crossref={cx_art}")

    # 第一作者姓氏
    m_au = re.match(r"([A-Z][A-Za-z\-' ]+?)(?:,| [A-Z])", body)
    cx_first = ((msg.get("author") or [{}])[0].get("family") or "")
    if m_au and cx_first:
        ms_surname = m_au.group(1).strip().split()[-1]
        if norm(ms_surname) != norm(re.sub(r"[,.]", "", cx_first)):
            fails.append(f"first author '{ms_surname}' vs crossref '{cx_first}'")

    return (not fails), fails


def main() -> int:
    md = open(MD, encoding="utf-8").read()
    refs = parse_refs(md)
    print(f"parsed {len(refs)} references")
    n_fail = 0
    for num, body in refs:
        if "doi:" in body.lower():
            try:
                ok, fails = check_crossref(body)
            except Exception as e:
                ok, fails = False, [f"checker error: {type(e).__name__}: {e}"]
        else:
            # 无 DOI：URL 可达性
            url = re.search(r"https?://\S+", body)
            ok, fails = False, ["no DOI"]
            if url:
                try:
                    r = requests.get(url.group(0), timeout=30,
                                     headers={"User-Agent": "ref-audit/1.0"})
                    ok = r.status_code == 200
                    fails = [] if ok else [f"URL status {r.status_code}"]
                except requests.RequestException as e:
                    fails = [f"URL unreachable: {type(e).__name__}"]
        status = "PASS" if ok else "FAIL"
        if not ok:
            n_fail += 1
        print(f"[{status}] ref {num}: {body[:68]}...")
        for f in fails:
            print(f"        - {f}")
    print("=" * 60)
    print(f"TOTAL FAILURES: {n_fail}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

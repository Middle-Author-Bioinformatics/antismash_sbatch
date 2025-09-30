#!/usr/bin/env python3
import re
from pathlib import Path
from html.parser import HTMLParser
import re
import csv
from pathlib import Path
import argparse
from typing import List, Dict, Optional, Tuple
from Bio import SeqIO  # Requires Biopython
import csv
import re
from typing import Optional, List


from typing import List, Dict
from Bio import SeqIO
import re

def _fmt_tfbs_maincol(tfbs_rows):
    """
    Build compact TFBS summary for the main CSV.
    Each site becomes: <type>(<strength>)(<orf-call>)
    Fallbacks ensure we render something even if a field name differs.
    """
    labels = []
    for r in tfbs_rows or []:
        tf_type   = (r.get("tfbs_type") or r.get("tf_type") or r.get("name") or "TFBS").strip()
        strength  = (r.get("strength") or r.get("site_strength") or r.get("score_class") or "NA").strip()
        # try to pick an ORF/locus tag that you already compute when writing the per-site TFBS CSV
        orf_call  = (r.get("orf") or r.get("orf_call") or r.get("locus_tag") or r.get("nearest_orf") or "?").strip()
        labels.append(f"{tf_type}({strength})({orf_call})")
    return "; ".join(labels)

class _MiniHTMLTable(HTMLParser):
    def __init__(self):
        super().__init__()
        self._in_table = False
        self._in_td = False
        self._cur_cell = []
        self._cur_row = []
        self.rows = []

    def handle_starttag(self, tag, attrs):
        t = tag.lower()
        if t == "table":
            if not self._in_table:
                self._in_table = True
        elif self._in_table and t in ("td", "th"):
            self._in_td = True
            self._cur_cell = []
        elif self._in_table and t in ("br",):
            if self._in_td:
                self._cur_cell.append("\n")

    def handle_endtag(self, tag):
        t = tag.lower()
        if t == "table" and self._in_table:
            self._in_table = False
        elif self._in_table and t in ("td", "th"):
            txt = "".join(self._cur_cell).strip()
            self._cur_row.append(txt)
            self._in_td = False
        elif self._in_table and t == "tr":
            if self._cur_row:
                self.rows.append(self._cur_row)
            self._cur_row = []

    def handle_data(self, data):
        if self._in_table and self._in_td and data:
            self._cur_cell.append(data)


def _iter_mibig_hit_html_files(kcb_root: Path):
    if not kcb_root or not kcb_root.exists():
        return
    for region_dir in sorted(kcb_root.glob("region*")):
        if not region_dir.is_dir():
            continue
        for hp in sorted(region_dir.glob("*_mibig_hits.html")):
            yield region_dir.name, hp

_mibig_re = re.compile(r"^(?P<orf>.+?)_mibig_hits\.html$", re.IGNORECASE)

def _orf_from_html_filename(html_path: Path) -> str:
    m = _mibig_re.match(html_path.name)
    return m.group("orf") if m else html_path.stem

def _parse_mibig_hits_table(html_text: str):
    p = _MiniHTMLTable()
    p.feed(html_text)
    return p.rows

def _write_mibig_hits_csv(out_prefix: Path, kcb_root: Path, verbose: bool=False) -> Path:
    import csv
    out_csv = out_prefix.with_suffix("")
    out_csv = out_csv.parent / (out_csv.name + ".mibig_hits.csv")
    n_rows = 0

    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["ORF", "MIBiG Protein", "Description", "MIBiG Cluster", "MiBiG Product"])

        for region_name, hp in _iter_mibig_hit_html_files(kcb_root):
            print(region_name, hp)
            try:
                text = hp.read_text(encoding="utf-8", errors="replace")
            except Exception as e:
                if verbose:
                    print(f"[MiBIG] WARN: cannot read {hp}: {e}")
                continue

            rows = _parse_mibig_hits_table(text)
            if not rows or len(rows) < 2:
                if verbose:
                    print(f"[MiBIG] INFO: no data rows in {hp}")
                continue

            # find header row if present; otherwise assume first row is header
            header_idx = 0
            for i, r in enumerate(rows[:5]):
                hdr_line = "|".join([c.lower() for c in r])
                if ("mibig protein" in hdr_line) and ("description" in hdr_line) and ("mibig cluster" in hdr_line):
                    header_idx = i
                    break

            data_rows = rows[header_idx+1:]
            orf = _orf_from_html_filename(hp)

            for r in data_rows:
                if len(r) < 4:
                    if verbose:
                        print(f"[MiBIG] WARN: malformed row in {hp}: {r!r}")
                    continue
                w.writerow([orf, r[0].strip(), r[1].strip(), r[2].strip(), r[3].strip()])
                n_rows += 1

    if verbose:
        print(f"[MiBIG] Wrote {n_rows} rows -> {out_csv}")
    return out_csv

def _resolve_knownclusterblast_dir(args, verbose: bool=False) -> Optional[Path]:
    """Always look for a folder named 'knownclusterblast' near the run dir.
    Priority:
      1) <run_dir>/knownclusterblast   (if you pass the run dir positionally)
      2) dirname(-o)/knownclusterblast (if -o points inside the run dir)
      3) ./knownclusterblast           (current working directory)
    """
    # candidates: list[Path] = []
    candidates: List[Path] = []

    # Try common positional names your v2 might be using (adjust freely):
    for attr in ("run_dir", "input", "indir", "antismash_dir", "run", "path"):
        if hasattr(args, attr) and getattr(args, attr):
            rd = Path(getattr(args, attr))
            candidates.append(rd / "knownclusterblast")

    # Output directory sibling (useful if you run "-o PSMB4MC5.csv PSMB4MC5")
    if hasattr(args, "output") and args.output:
        candidates.append(Path(args.output).parent / "knownclusterblast")

    # Fallback: CWD
    candidates.append(Path.cwd() / "knownclusterblast")

    for p in candidates:
        if p.exists() and p.is_dir():
            if verbose:
                print(f"[MiBIG] using knownclusterblast at: {p}")
            return p

    if verbose:
        tried = "\n  ".join(str(x) for x in candidates)
        print(f"[MiBIG] INFO: could not find 'knownclusterblast'. Tried:\n  {tried}")
    return None

def _nearest_locus_tag(rec, tfbs_feat) -> str:
    """Find the overlapping or nearest CDS and return its locus_tag/gene if present."""
    try:
        tfbs_start = int(tfbs_feat.location.start)
        tfbs_end = int(tfbs_feat.location.end)
        tfbs_mid = (tfbs_start + tfbs_end) // 2
    except Exception:
        return ""

    best = ("", 10**12)  # (locus_tag, distance)
    for feat in rec.features:
        if feat.type != "CDS":
            continue
        quals = feat.qualifiers or {}
        locus = " ".join(quals.get("locus_tag", [])) or " ".join(quals.get("gene", []))
        if not locus:
            continue
        try:
            cds_start = int(feat.location.start)
            cds_end = int(feat.location.end)
        except Exception:
            continue

        if not (tfbs_end < cds_start or tfbs_start > cds_end):
            # overlap
            return locus

        # distance to nearest edge
        if tfbs_mid < cds_start:
            dist = cds_start - tfbs_mid
        elif tfbs_mid > cds_end:
            dist = tfbs_mid - cds_end
        else:
            dist = 0

        if dist < best[1]:
            best = (locus, dist)

    return best[0]


def _extract_tfbs_rows(gbk_path: "Path", verbose: bool=False) -> List[Dict[str, str]]:
    """
    Build per-site TFBS rows: contig,locus_tag,tfbs_type,strength.
    - tfbs_type: parsed from notes like 'TFBS match to X'
    - strength: first try /score or /binding_site_strength; else 'confidence' if present; else blank
    - locus_tag: from feature if present; else nearest/overlapping CDS
    """
    out: List[Dict[str, str]] = []
    try:
        for rec in SeqIO.parse(str(gbk_path), "genbank"):
            contig = rec.id
            for feat in rec.features:
                ftype = (feat.type or "").lower()
                quals = feat.qualifiers or {}

                # Notes can be a list; flatten for searches
                notes = " ".join(quals.get("note", []))

                # Detect TFBS features
                is_tfbs = False
                if "tfbs match to" in notes.lower():
                    is_tfbs = True
                elif ftype in {"tf_binding_site"}:
                    is_tfbs = True
                elif ftype == "regulatory":
                    rc = " ".join(quals.get("regulatory_class", [])).lower()
                    if "tf" in rc and "binding" in rc:
                        is_tfbs = True

                if not is_tfbs:
                    continue

                # Parse TFBS type (name after 'TFBS match to')
                tfbs_type = ""
                m = re.search(r"TFBS match to\s+(.+?)(?:[,;]|\)|$)", notes, flags=re.IGNORECASE | re.DOTALL)
                if m:
                    tfbs_type = m.group(1).strip()

                # Parse strength/score/confidence
                strength = (
                    " ".join(quals.get("score", [])) or
                    " ".join(quals.get("binding_site_strength", []))
                ).strip()

                if not strength:
                    # fall back to confidence field in notes
                    m2 = re.search(r"confidence:\s*([a-zA-Z]+)", notes, flags=re.IGNORECASE)
                    if m2:
                        strength = m2.group(1).strip().lower()

                # Determine locus_tag
                locus_tag = (
                    " ".join(quals.get("locus_tag", [])) or
                    " ".join(quals.get("gene", []))
                ).strip()
                if not locus_tag:
                    locus_tag = _nearest_locus_tag(rec, feat) or ""

                out.append({
                    "contig": contig,
                    "locus_tag": locus_tag,
                    "tfbs_type": tfbs_type,
                    "strength": strength,
                })
    except Exception as e:
        if verbose:
            print(f"[TFBS] Failed to parse {gbk_path.name}: {e}")
    return out


def _read_text(path: Path) -> str:
    return path.read_text(errors="ignore")

def _kcb_find_blocks(text: str):
    pat = re.compile(r'^\s*(?P<marker>>>\s*)?(?P<rank>\d+)\.\s*(?P<bgc>BGC[0-9.]+)\s*$', re.MULTILINE)
    blocks = []
    headers = list(pat.finditer(text))
    for i, m in enumerate(headers):
        start = m.end()
        end = headers[i+1].start() if i+1 < len(headers) else len(text)
        blocks.append({
            "topflag": bool(m.group('marker')),
            "rank": int(m.group('rank')),
            "bgc": m.group('bgc'),
            "text": text[start:end]
        })
    return blocks

def _kcb_parse_details(block_text: str) -> Dict[str, str]:
    def grab(pattern):
        mm = re.search(pattern, block_text, re.IGNORECASE | re.MULTILINE)
        return mm.group(1).strip() if mm else ""
    return {
        "BGC Source": grab(r'^\s*Source:\s*(.+)$'),
        "BGC Type": grab(r'^\s*Type:\s*(.+)$'),
        "BGC Proteins with BLAST hits": grab(r'^\s*Number of proteins with BLAST hits to this cluster:\s*([\d.]+)'),
        "BGC Cumulative BLAST score": grab(r'^\s*Cumulative BLAST score:\s*([\d.]+)'),
    }

def _parse_knownclusterblast(txt_path: Optional[Path], verbose: bool=False) -> Dict[str, str]:
    """Parse KnownClusterBlast text and return *all* significant hits
    concatenated with semicolons in the BGC columns (Accession, Source, Type,
    Proteins with BLAST hits, Cumulative BLAST score). Preserves original behavior
    if no blocks are found.
    """
    empty = {
        "BGC Accession": "",
        "BGC Source": "",
        "BGC Type": "",
        "BGC Proteins with BLAST hits": "",
        "BGC Cumulative BLAST score": "",
    }
    if not txt_path or not txt_path.exists():
        if verbose:
            print(f"[KCB] No file found for mapping.")
        return empty

    text = _read_text(txt_path)
    blocks = _kcb_find_blocks(text)
    if not blocks:
        if verbose:
            print(f"[KCB] No blocks found in {txt_path}")
        return empty

    # Sort by 'rank' (1..N); '>>' markers keep relative order.
    blocks = sorted(blocks, key=lambda b: b.get("rank", 10**9))

    accs: List[str] = []
    srcs: List[str] = []
    types: List[str] = []
    prots: List[str] = []
    scores: List[str] = []

    for b in blocks:
        accs.append(b.get("bgc", ""))
        d = _kcb_parse_details(b.get("text", ""))
        srcs.append(d.get("BGC Source", ""))
        types.append(d.get("BGC Type", ""))
        prots.append(d.get("BGC Proteins with BLAST hits", ""))
        scores.append(d.get("BGC Cumulative BLAST score", ""))

    def join_semicol(vals: List[str]) -> str:
        cleaned = [v.strip() for v in vals if v is not None and str(v).strip() != ""]
        return "; ".join(cleaned)

    out = {
        "BGC Accession": join_semicol(accs),
        "BGC Source": join_semicol(srcs),
        "BGC Type": join_semicol(types),
        "BGC Proteins with BLAST hits": join_semicol(prots),
        "BGC Cumulative BLAST score": join_semicol(scores),
    }
    if verbose:
        print(f"[KCB] Aggregated {len(blocks)} hits from {txt_path.name}")
    return out

def _find_knowncluster_dirs(root: Path) -> List[Path]:
    # Case-insensitive search for 'knownclusterblast' near the GBKs
    candidates = []
    for base in [root, root.parent, root.parent.parent if root.parent else None]:
        if not base: continue
        for p in base.iterdir():
            if p.is_dir() and 'knownclusterblast' in p.name.lower():
                candidates.append(p)
    return list(dict.fromkeys(candidates))  # de-dup

def _find_knowncluster_file(gbk_file: Path, known_dirs: List[Path]) -> Optional[Path]:
    stem = gbk_file.stem  # e.g., k105_371.region001
    m = re.search(r'^(?P<prefix>.+)\.region(?P<reg>\d+)$', stem)
    prefix = m.group('prefix') if m else stem
    reg = int(m.group('reg')) if m else None

    names = [f"{stem}.txt"]
    if reg is not None:
        names += [f"{prefix}_c{reg}.txt", f"{prefix}_c{int(reg)}.txt"]
    globs = [f"{stem}*.txt"]
    if reg is not None:
        globs += [f"{prefix}_c{reg}*.txt"]
    globs += [f"{prefix}*.txt"]

    for kdir in known_dirs:
        for n in names:
            cand = kdir / n
            if cand.exists():
                return cand
        for pat in globs:
            found = list(kdir.glob(pat))
            if found:
                return found[0]

    for n in names:
        cand = gbk_file.parent / n
        if cand.exists():
            return cand
    if reg is not None:
        cand = gbk_file.parent / f"{prefix}_c{reg}.txt"
        if cand.exists():
            return cand
    for pat in globs:
        for cand in gbk_file.parent.glob(pat):
            return cand
    return None

def parse_antismash_gbk(gbk_file: Path, known_dirs: List[Path], verbose: bool=False) -> Tuple[dict, List[Dict[str, str]]]:
    text = _read_text(gbk_file)
    lines = text.splitlines()

    region = ""
    m = re.search(r"\.region(\d+)\b", gbk_file.name)
    if m:
        try:
            region = int(m.group(1))
        except ValueError:
            region = m.group(1)

    total_len = ""
    m = re.search(r"^LOCUS\s+\S+\s+(\d+)\s+bp", "\n".join(lines[:25]), re.MULTILINE)
    if m:
        total_len = int(m.group(1))

    header_block = "\n".join(lines[:300])

    def grab(pattern: str, block=header_block, flags=0, default=""):
        mm = re.search(pattern, block, flags)
        return mm.group(1) if mm else default

    product = grab(r'/product="([^"]+)"')
    category = grab(r'/category="([^"]+)"')
    contig_edge = grab(r'/contig_edge="([^"]+)"')

    # NEW: structured per-site TFBS rows (used for both the per-site CSV and the compact main-CSV column)
    tfbs_rows = _extract_tfbs_rows(gbk_file, verbose=verbose)
    tfbs_types = _fmt_tfbs_maincol(tfbs_rows)   # -> "Type(strength)(orf); Type2(strength)(orf2)"
    tfbs_hits = len(tfbs_rows)

    # KnownClusterBlast mapping
    kcb_file = _find_knowncluster_file(gbk_file, known_dirs)
    if verbose:
        print(f"[MAP] {gbk_file.name} -> {kcb_file if kcb_file else 'NO MATCH'}")
    kcb = _parse_knownclusterblast(kcb_file, verbose=verbose) if kcb_file else {
        "BGC Accession": "",
        "BGC Source": "",
        "BGC Type": "",
        "BGC Proteins with BLAST hits": "",
        "BGC Cumulative BLAST score": "",
    }

    row = {
        "File": gbk_file.name,
        "Region": region if region != "" else 1,
        "Product": product,
        "Category": category,
        "Contig Edge": contig_edge,
        "Total Length (bp)": total_len,
        "BGC Accession": kcb["BGC Accession"],
        "BGC Source": kcb["BGC Source"],
        "BGC Type": kcb["BGC Type"],
        "BGC Proteins with BLAST hits": kcb["BGC Proteins with BLAST hits"],
        "BGC Cumulative BLAST score": kcb["BGC Cumulative BLAST score"],
        # TFBS LAST:
        "TFBS Hits": tfbs_hits,
        "TFBS Types": tfbs_types,
    }

    return row, tfbs_rows


def _collect_known_dirs(dir_path: Path, gbks: List[Path]) -> List[Path]:
    known_dirs = _find_knowncluster_dirs(dir_path)
    if not known_dirs:
        for g in gbks[:1]:
            known_dirs = _find_knowncluster_dirs(g.parent)
            if known_dirs:
                break
    return known_dirs

def parse_directory(
    dir_path: Path,
    out_csv: Path,
    known_dir_override: Optional[Path]=None,
    verbose: bool=False,
    tfbs_csv: Optional[Path]=None,
) -> None:
    rows = []
    all_tfbs_rows: List[Dict[str, str]] = []

    gbks = [g for g in (sorted(dir_path.glob("*.gbk")) or sorted(dir_path.rglob("*.gbk")))
            if g.name != f"{dir_path.name}.gbk"]

    if known_dir_override and known_dir_override.exists():
        known_dirs = [known_dir_override]
    else:
        known_dirs = _collect_known_dirs(dir_path, gbks)

    if verbose:
        print("[INFO] knownclusterblast dirs:", [str(p) for p in known_dirs] or ["<none>"])

    for gbk in gbks:
        try:
            parsed = parse_antismash_gbk(gbk, known_dirs, verbose=verbose)

            # Support: tuple return (row, tfbs_rows)
            if isinstance(parsed, tuple) and len(parsed) == 2:
                row, tfbs_rows = parsed
                rows.append(row)
                if tfbs_rows:
                    all_tfbs_rows.extend(tfbs_rows)
            else:
                # If the function hasn't been updated, keep working and extract TFBS here
                row = parsed
                rows.append(row)
                all_tfbs_rows.extend(_extract_tfbs_rows(gbk, verbose=verbose))

        except Exception as e:
            if verbose:
                print(f"[ERROR] {gbk.name}: {e}")
            rows.append({
                "File": gbk.name,
                "Region": "",
                "Product": "",
                "Category": "",
                "Contig Edge": "",
                "Total Length (bp)": "",
                "BGC Accession": "",
                "BGC Source": "",
                "BGC Type": "",
                "BGC Proteins with BLAST hits": "",
                "BGC Cumulative BLAST score": "",
                "TFBS Hits": "",
                "TFBS Types": f"Error: {e}",
            })

    # Write main CSV (unchanged)
    fieldnames = [
        "File", "Region", "Product", "Category", "Contig Edge",
        "Total Length (bp)",
        "BGC Accession", "BGC Source", "BGC Type",
        "BGC Proteins with BLAST hits", "BGC Cumulative BLAST score",
        "TFBS Hits", "TFBS Types"
    ]
    with out_csv.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    # Write TFBS CSV (new)
    tfbs_out = tfbs_csv or out_csv.with_name(out_csv.stem + ".tfbs.csv")
    tfbs_fields = ["contig", "locus_tag", "tfbs_type", "strength"]
    with tfbs_out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=tfbs_fields)
        writer.writeheader()
        for r in all_tfbs_rows:
            writer.writerow({
                "contig": r.get("contig", ""),
                "locus_tag": r.get("locus_tag", ""),
                "tfbs_type": r.get("tfbs_type", ""),
                "strength": r.get("strength", ""),
            })

    if verbose:
        print(f"[OK] Wrote main CSV: {out_csv}")
        print(f"[OK] Wrote TFBS CSV: {tfbs_out}")


def main():
    ap = argparse.ArgumentParser(description="Summarize antiSMASH .gbk files (KnownClusterBlast + TFBS).")
    ap.add_argument("gbk_directory", help="Directory containing .gbk files")
    ap.add_argument("-o", "--output", default="summary_with_knownclusterblast.csv", help="Output CSV path")
    ap.add_argument("-k", "--knowncluster-dir", default=None, help="Path to knownclusterblast directory (optional)")
    ap.add_argument("-v", "--verbose", action="store_true", help="Verbose mapping/debug output")
    args = ap.parse_args()

    out_path = Path(args.output)  # your main CSV path

    parse_directory(
        Path(args.gbk_directory),
        Path(args.output),
        known_dir_override=Path(args.knowncluster_dir) if args.knowncluster_dir else None,
        verbose=args.verbose
    )

    out_path = Path(args.output)  # e.g., HSM67VEI.csv
    run_id = out_path.with_suffix("").name  # -> "HSM67VEI"

    # Hardcode your layout: ./HSM67VEI/knownclusterblast
    kcb_root = Path.cwd() / run_id / "knownclusterblast"
    if not kcb_root.is_dir():
        # If you cd into the run dir, fallback to ./knownclusterblast
        kcb_root = Path("knownclusterblast")

    # Debug prints (optional)
    if getattr(args, "verbose", False):
        print(f"[MiBIG] kcb_root candidate: {kcb_root} exists={kcb_root.is_dir()}")
        count = len(list(kcb_root.glob("region*/**/*_mibig_hits.html"))) if kcb_root.is_dir() else 0
        print(f"[MiBIG] found {count} mibig html(s)")

    # Always attempt to write the CSV if the folder exists
    if kcb_root.is_dir():
        _write_mibig_hits_csv(out_prefix=out_path, kcb_root=kcb_root, verbose=getattr(args, "verbose", False))
    else:
        # Fail loudly so this can't silently return None again
        raise SystemExit(f"[MiBIG] ERROR: expected knownclusterblast at {kcb_root} but it does not exist")

if __name__ == "__main__":
    main()




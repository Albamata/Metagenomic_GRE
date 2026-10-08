#!/usr/bin/env python3
"""
01_aggregate_orfs_by_ko.py — Agrupa ORFs por KEGG ID para la Fracción Celular (FC) y la
Fracción Vesicular (FV), calcula TPM con denominador compartido y genera
tres archivos de salida:

  • <salida_tpm>.txt       → Raw + TPM por estación, una fila FC y una FV por KEGG
                             Incluye gene_name tomado del TXT KEGG.
  • <salida_tax>.tsv       → KEGG_ID × Origen × Estación × Taxón único (formato largo para R)
  • <salida_kegg_ref>.tsv  → KEGG_ID × gene_name × función × ruta(s) del TXT KEGG

Uso:
python3 scripts/01_aggregate_orfs_by_ko.py --fc annotation/13.FC_Plan_nacional.orftable.prok.tsv --fv VESICULAS/13.Z1_V14.orftable VESICULAS/13.Z2_V10.orftable VESICULAS/13.Z3_V6.orftable VESICULAS/13.Z4_V2.orftable --kegg_txt SLACK/docs/keggfun2.txt --out_tpm GUADALQUIVIR/resultado/kegg_tpm.tsv --out_tax GUADALQUIVIR/resultado/kegg_tax.tsv --out_kegg_ref GUADALQUIVIR/resultado/kegg_ref.tsv


Notas:
  • FC ya está pre-filtrado a procariontes.
  • FV se filtra a Molecule=CDS y Tax que empiece por k_Bacteria o k_Archaea.
  • El denominador TPM es compartido: suma RPK de todos los CDS válidos de FC + FV.
  • Cada archivo FV tiene una única estación; el nombre se infiere del nombre
    del archivo (e.g. 13.Z2_V10.orftable → Z2_V10).
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

# ---------------------------------------------------------------------------
# Constantes
# ---------------------------------------------------------------------------

# Estaciones FC (en el orftable grande todas juntas)
FC_STATIONS: list[str] = [
    "Z1_15", "Z1_14", "Z1_13",
    "Z2_11", "Z2_10", "Z2_9",
    "Z3_7",  "Z3_6",  "Z3_5",  "Z3_2",
]

# Prefijos taxonómicos aceptados para FV
PROK_PREFIXES: tuple[str, ...] = ("k_Bacteria", "k_Archaea")

# Columnas del orftable (comunes a FC y FV)
COL_ORF      = "ORF ID"
COL_MOLECULE = "Molecule"
COL_LENGTH   = "Length NT"
COL_TAX      = "Tax"
COL_KEGG_ID   = "KEGG ID"
COL_KEGGFUN   = "KEGGFUN"
COL_KEGGPATH  = "KEGGPATH"
RAW_PREFIX   = "Raw read count "

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _to_float(val: str) -> float | None:
    try:
        return float(val.strip())
    except (ValueError, AttributeError):
        return None


def _clean(val: str) -> str:
    return val.strip().rstrip("*").strip()


def _join_unique(values: list[str]) -> str:
    seen: dict[str, None] = {}
    for v in values:
        v = v.strip()
        if v:
            seen[v] = None
    return "; ".join(seen)


def _infer_station(path: Path) -> str:
    """
    Infiere el nombre de la estación desde el nombre del archivo FV.
    '13.Z2_V10.orftable' → 'Z2_V10'
    """
    m = re.search(r"(Z\d+_V\d+)", path.name)
    if m:
        return m.group(1)
    # Fallback: usar el stem completo
    return path.stem


# Categorías principales que aparecen en las rutas KEGG.
# Se usan solo como plan B si el TXT no está perfectamente tabulado.
KEGG_MAIN_CATEGORIES: tuple[str, ...] = (
    "Metabolism",
    "Genetic Information Processing",
    "Environmental Information Processing",
    "Cellular Processes",
    "Organismal Systems",
    "Human Diseases",
    "Drug Development",
    "Brite Hierarchies",
)


def _split_pathway_levels(pathway: str) -> dict[str, str]:
    """Divide una ruta KEGG completa en niveles jerárquicos."""
    levels = [x.strip() for x in pathway.split(";") if x.strip()]
    return {
        "Pathway_L1": levels[0] if len(levels) > 0 else "",
        "Pathway_L2": levels[1] if len(levels) > 1 else "",
        "Pathway_L3": levels[2] if len(levels) > 2 else "",
        "Pathway_L4": levels[3] if len(levels) > 3 else "",
    }


def parse_kegg_txt(path: Path) -> dict[str, dict]:
    """
    Lee un TXT de KEGG tipo:
      K00001    adh    alcohol dehydrogenase [EC:1.1.1.1]    Metabolism; ... | Metabolism; ...

    Devuelve:
      dict[KEGG_ID] = {
          "gene_name": str,
          "function_txt": str,
          "pathways": list[str],
      }

    Importante:
      - gene_name, function_txt y pathways salen del TXT, no de la orftable.
      - Si hay varias rutas para un KO, se conservan todas.
    """
    kegg_ref: dict[str, dict] = {}
    category_regex = r"\s+(?=(" + "|".join(map(re.escape, KEGG_MAIN_CATEGORIES)) + r");)"

    with path.open(encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or not line.startswith("K"):
                continue

            # Caso ideal: columnas separadas por tabuladores o 2+ espacios.
            # K00002 | AKR1A1, adh | alcohol dehydrogenase... | Metabolism; ...
            fields = re.split(r"\t+|\s{2,}", line, maxsplit=3)

            kegg_id = ""
            gene_name = ""
            function_txt = ""
            path_text = ""

            if len(fields) >= 4:
                kegg_id, gene_name, function_txt, path_text = [x.strip() for x in fields[:4]]
            elif len(fields) == 3:
                kegg_id, gene_name, function_txt = [x.strip() for x in fields]
            elif len(fields) in (1, 2):
                # Plan B para líneas que no vengan bien tabuladas:
                # separar K del resto y localizar el inicio de la ruta.
                if len(fields) == 2:
                    kegg_id, rest = fields[0].strip(), fields[1].strip()
                else:
                    m_ko = re.match(r"^(K\d{5})\s+(.+)$", fields[0].strip())
                    if not m_ko:
                        continue
                    kegg_id, rest = m_ko.group(1), m_ko.group(2).strip()

                m = re.search(category_regex, rest)
                if m:
                    annot = rest[:m.start()].strip()
                    path_text = rest[m.start():].strip()
                else:
                    annot = rest

                # Plan B conservador:
                # - Si hay alias separados por coma, conserva todo el bloque de genes.
                #   Ej.: "AKR1A1, adh alcohol dehydrogenase" → gene_name = "AKR1A1, adh"
                # - Si no hay coma, usa el primer token como gene_name.
                m_gene_alias = re.match(r"^(.+?,\s*\S+)\s+(.+)$", annot)
                if m_gene_alias:
                    gene_name = m_gene_alias.group(1).strip().rstrip(",")
                    function_txt = m_gene_alias.group(2).strip()
                else:
                    annot_parts = annot.split(maxsplit=1)
                    if len(annot_parts) == 1:
                        gene_name = annot_parts[0].strip().rstrip(",")
                    elif len(annot_parts) == 2:
                        gene_name = annot_parts[0].strip().rstrip(",")
                        function_txt = annot_parts[1].strip()
            else:
                continue

            kegg_id = _clean(kegg_id)
            gene_name = gene_name.strip().rstrip("*").strip()
            function_txt = function_txt.strip().rstrip("*").strip()

            pathways = [p.strip() for p in path_text.split("|") if p.strip()]
            pathways = list(dict.fromkeys(pathways))  # quitar duplicados manteniendo orden

            if kegg_id:
                kegg_ref[kegg_id] = {
                    "gene_name": gene_name,
                    "function_txt": function_txt,
                    "pathways": pathways,
                }

    return kegg_ref


def _make_bucket() -> dict:
    return {
        "orfs":       [],
        "lengths":    [],
        "raw":        defaultdict(float),
        "keggfun":     [],
        "keggpath":    [],
        # taxa_by_station: dict[station, set[str]]
        # Permite reconstruir kegg × Estación × Taxón único
        "taxa_by_st": defaultdict(set),
    }

# ---------------------------------------------------------------------------
# Lectura FC
# ---------------------------------------------------------------------------

def read_fc(path: Path) -> tuple[dict[str, dict], float]:
    """
    Lee el orftable de Fracción Celular (ya pre-filtrado a procariontes).
    Devuelve:
        kegg_data  : dict[kegg_id] → bucket
        rpk_total : contribución al denominador TPM compartido
    """
    kegg_data: dict[str, dict] = defaultdict(_make_bucket)
    rpk_total = 0.0

    with path.open(encoding="utf-8", errors="ignore") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        header = reader.fieldnames or []

        _check_cols(header, [COL_ORF, COL_MOLECULE, COL_LENGTH, COL_KEGG_ID,
                              COL_KEGGFUN, COL_KEGGPATH, COL_TAX], path)

        raw_cols = {st: f"{RAW_PREFIX}{st}" for st in FC_STATIONS}
        missing_raw = [v for v in raw_cols.values() if v not in header]
        if missing_raw:
            sys.exit(f"[ERROR] Columnas Raw no encontradas en FC: {missing_raw}")

        n_cds = 0
        for row in reader:
            if row.get(COL_MOLECULE, "").strip() != "CDS":
                continue
            n_cds += 1

            length = _to_float(row[COL_LENGTH])
            if length and length > 0:
                length_kb = length / 1000.0
                for st in FC_STATIONS:
                    raw_val = _to_float(row[raw_cols[st]]) or 0.0
                    rpk_total += raw_val / length_kb

            kegg_id = _clean(row[COL_KEGG_ID])
            if not kegg_id:
                continue

            b = kegg_data[kegg_id]
            b["orfs"].append(row[COL_ORF].strip())
            if length:
                b["lengths"].append(length)
            for st in FC_STATIONS:
                b["raw"][st] += _to_float(row[raw_cols[st]]) or 0.0
            if row[COL_KEGGFUN].strip():
                b["keggfun"].append(row[COL_KEGGFUN])
            if row[COL_KEGGPATH].strip():
                b["keggpath"].append(row[COL_KEGGPATH])
            tax = row[COL_TAX].strip()
            if tax:
                # Solo añadir el taxón a las estaciones donde este ORF
                # tiene raw reads > 0, así la variación por estación es real
                for st in FC_STATIONS:
                    if (_to_float(row[raw_cols[st]]) or 0.0) > 0:
                        b["taxa_by_st"][st].add(tax)

    print(f"  [FC] CDS leídos: {n_cds} | KEGG IDs únicos: {len(kegg_data)}")
    return dict(kegg_data), rpk_total

# ---------------------------------------------------------------------------
# Lectura FV (un archivo por estación)
# ---------------------------------------------------------------------------

def read_fv(paths: list[Path]) -> tuple[dict[str, dict], float, list[str]]:
    """
    Lee los orftables de Fracción Vesicular (uno por estación).
    Filtra a Molecule=CDS y Tax que empiece por k_Bacteria / k_Archaea.
    Devuelve:
        kegg_data   : dict[kegg_id] → bucket  (raw keyed por nombre de estación)
        rpk_total  : contribución al denominador TPM compartido
        stations   : lista ordenada de nombres de estación inferidos
    """
    kegg_data: dict[str, dict] = defaultdict(_make_bucket)
    rpk_total = 0.0
    stations: list[str] = []

    for path in paths:
        station = _infer_station(path)
        stations.append(station)
        print(f"  [FV] Leyendo {path.name} → estación '{station}'")

        with path.open(encoding="utf-8", errors="ignore") as fh:
            # Saltar líneas de comentario que empiezan por '#'
            lines = (ln for ln in fh if not ln.startswith("#"))
            reader = csv.DictReader(lines, delimiter="\t")
            header = reader.fieldnames or []

            _check_cols(header, [COL_ORF, COL_MOLECULE, COL_LENGTH, COL_KEGG_ID,
                                  COL_KEGGFUN, COL_KEGGPATH, COL_TAX], path)

            raw_col = f"{RAW_PREFIX}{station}"
            if raw_col not in header:
                sys.exit(f"[ERROR] Columna '{raw_col}' no encontrada en {path.name}\n"
                         f"  Header: {header}")

            n_cds = n_prok = 0
            for row in reader:
                if row.get(COL_MOLECULE, "").strip() != "CDS":
                    continue
                n_cds += 1

                tax = row[COL_TAX].strip()
                if not any(tax.startswith(p) for p in PROK_PREFIXES):
                    continue
                n_prok += 1

                length = _to_float(row[COL_LENGTH])
                if length and length > 0:
                    length_kb = length / 1000.0
                    raw_val = _to_float(row[raw_col]) or 0.0
                    rpk_total += raw_val / length_kb

                kegg_id = _clean(row[COL_KEGG_ID])
                if not kegg_id:
                    continue

                b = kegg_data[kegg_id]
                b["orfs"].append(row[COL_ORF].strip())
                if length:
                    b["lengths"].append(length)
                b["raw"][station] += _to_float(row[raw_col]) or 0.0
                if row[COL_KEGGFUN].strip():
                    b["keggfun"].append(row[COL_KEGGFUN])
                if row[COL_KEGGPATH].strip():
                    b["keggpath"].append(row[COL_KEGGPATH])
                if tax:
                    # En FV cada archivo tiene su propia estación
                    b["taxa_by_st"][station].add(tax)

            print(f"         CDS totales: {n_cds} | CDS prok retenidos: {n_prok}")

    print(f"  [FV] KEGG IDs únicos: {len(kegg_data)}")
    return dict(kegg_data), rpk_total, stations

# ---------------------------------------------------------------------------
# Helpers de validación
# ---------------------------------------------------------------------------

def _check_cols(header: list[str], required: list[str], path: Path) -> None:
    missing = [c for c in required if c not in header]
    if missing:
        sys.exit(f"[ERROR] Columnas no encontradas en {path.name}: {missing}\n"
                 f"  Header disponible: {header}")

# ---------------------------------------------------------------------------
# Construcción de tablas de salida
# ---------------------------------------------------------------------------

def build_tpm_rows(
    fc_data: dict[str, dict],
    fv_data: dict[str, dict],
    rpk_denom: float,
    fv_stations: list[str],
    kegg_ref: dict[str, dict],
) -> list[dict]:
    """Genera filas para la tabla TPM (una FC + una FV por KEGG_ID)."""

    all_keggs = sorted(set(fc_data) | set(fv_data))
    rows = []

    for kegg_id in all_keggs:
        for origen, data, stations in (
            ("FC", fc_data, FC_STATIONS),
            ("FV", fv_data, fv_stations),
        ):
            if kegg_id not in data:
                continue

            d = data[kegg_id]
            lengths = d["lengths"]
            raw     = d["raw"]
            mean_len_kb = (sum(lengths) / len(lengths) / 1000.0) if lengths else None

            ref = kegg_ref.get(kegg_id, {})

            row: dict = {
                "KEGG_ID": kegg_id,
                "gene_name": ref.get("gene_name", ""),
                "Origen": origen,
                "ORFs":   "; ".join(d["orfs"]),
                "n_ORFs": len(d["orfs"]),
            }

            for st in stations:
                row[f"Raw_{st}"] = raw.get(st, 0.0)

            for st in stations:
                raw_val = raw.get(st, 0.0)
                if mean_len_kb and mean_len_kb > 0 and rpk_denom > 0:
                    rpk = raw_val / mean_len_kb
                    row[f"TPM_{st}"] = rpk / rpk_denom * 1e6
                else:
                    row[f"TPM_{st}"] = float("nan")

            row["KEGGFUN"]  = _join_unique(d["keggfun"])
            row["KEGGPATH"] = _join_unique(d["keggpath"])
            rows.append(row)

    return rows


def build_kegg_ref_rows(
    fc_data: dict[str, dict],
    fv_data: dict[str, dict],
    kegg_ref: dict[str, dict],
) -> list[dict]:
    """
    Genera una tabla de referencia desde el TXT KEGG.

    Formato largo: una fila por KEGG_ID × ruta.
    Así no se pierde información cuando un mismo KEGG aparece en varias rutas.
    """
    rows = []
    all_keggs = sorted(set(fc_data) | set(fv_data))

    for kegg_id in all_keggs:
        ref = kegg_ref.get(kegg_id, {})
        gene_name = ref.get("gene_name", "")
        function_txt = ref.get("function_txt", "")
        pathways = ref.get("pathways", [])

        if not pathways:
            row = {
                "KEGG_ID": kegg_id,
                "gene_name": gene_name,
                "function_txt": function_txt,
                "Pathway_full": "",
                "Pathway_L1": "",
                "Pathway_L2": "",
                "Pathway_L3": "",
                "Pathway_L4": "",
            }
            rows.append(row)
            continue

        for pathway in pathways:
            row = {
                "KEGG_ID": kegg_id,
                "gene_name": gene_name,
                "function_txt": function_txt,
                "Pathway_full": pathway,
            }
            row.update(_split_pathway_levels(pathway))
            rows.append(row)

    return rows


def build_tpm_fieldnames(fv_stations: list[str]) -> list[str]:
    """
    Orden fijo para la tabla principal.
    Primero salen las columnas FC que pediste; si hay FV, se añaden detrás
    para no perder las estaciones vesiculares.
    """
    stations = FC_STATIONS + [st for st in fv_stations if st not in FC_STATIONS]
    return (
        ["KEGG_ID", "gene_name", "Origen", "ORFs", "n_ORFs"]
        + [f"Raw_{st}" for st in stations]
        + [f"TPM_{st}" for st in stations]
        + ["KEGGFUN", "KEGGPATH"]
    )


# Prefijos linneanos reconocidos (n_ se ignora)
_TAX_RANKS: dict[str, str] = {
    "k": "Kingdom",
    "p": "Phylum",
    "c": "Class",
    "o": "Order",
    "f": "Family",
    "g": "Genus",
    "s": "Species",
}
_NO_DEF = "No definido"


def parse_taxonomy(tax_str: str) -> dict[str, str]:
    """
    Parsea una cadena taxonómica de SqueezeMeta en columnas linneanas.
    - Prefijos reconocidos: k, p, c, o, f, g, s
    - Prefijo 'n_' (clade informal) → ignorado
    - Niveles no presentes → 'No definido'

    Ejemplo:
      'k_Archaea;n_TACK group;p_Nitrososphaerota;c_Nitrososphaeria'
      → {Kingdom: Archaea, Phylum: Nitrososphaerota, Class: Nitrososphaeria,
         Order: No definido, Family: No definido, Genus: No definido, Species: No definido}
    """
    result = {rank: _NO_DEF for rank in _TAX_RANKS.values()}
    for token in tax_str.split(";"):
        token = token.strip()
        if not token or "_" not in token:
            continue
        prefix, _, value = token.partition("_")
        prefix = prefix.strip().lower()
        if prefix in _TAX_RANKS:
            result[_TAX_RANKS[prefix]] = value.strip()
        # prefix == "n" → ignorado
    return result


def build_tax_rows(
    fc_data: dict[str, dict],
    fv_data: dict[str, dict],
) -> list[dict]:
    """
    Genera filas para la tabla taxonómica en formato largo:
        KEGG_ID | Origen | Estacion | Kingdom | Phylum | Class | Order | Family | Genus | Species
    Una fila por combinación única KEGG × Origen × Estación × Taxón.
    Óptimo para dplyr::group_by en R (riqueza, abundancia, etc.).
    """
    rows = []
    all_keggs = sorted(set(fc_data) | set(fv_data))

    for kegg_id in all_keggs:
        for origen, data in (("FC", fc_data), ("FV", fv_data)):
            if kegg_id not in data:
                continue
            taxa_by_st: dict[str, set] = data[kegg_id]["taxa_by_st"]
            for station in sorted(taxa_by_st):
                for tax in sorted(taxa_by_st[station]):
                    row = {
                        "KEGG_ID":   kegg_id,
                        "Origen":   origen,
                        "Estacion": station,
                    }
                    row.update(parse_taxonomy(tax))
                    rows.append(row)

    return rows

# ---------------------------------------------------------------------------
# Escritura TSV
# ---------------------------------------------------------------------------

def write_tsv(
    rows: list[dict],
    out_path: Path,
    label: str,
    fieldnames: list[str] | None = None,
) -> None:
    if not rows:
        print(f"[WARN] No hay filas para {label}.")
        return

    if fieldnames is None:
        fieldnames = list(rows[0].keys())

    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)

    print(f"[OK] {label}: {out_path}  ({len(rows)} filas)")

# ---------------------------------------------------------------------------
# Resumen TPM
# ---------------------------------------------------------------------------

def print_tpm_summary(rows: list[dict]) -> None:
    # Recoger todas las columnas TPM presentes
    tpm_cols = [k for k in rows[0] if k.startswith("TPM_")] if rows else []
    sums: dict[str, float] = defaultdict(float)
    counts: dict[str, int] = defaultdict(int)

    for row in rows:
        origen = row["Origen"]
        for col in tpm_cols:
            v = row.get(col, float("nan"))
            if v == v:   # NaN-safe
                sums[f"{origen}_{col}"] += v
                counts[f"{origen}_{col}"] += 1

    print("\nSuma TPM por fracción y estación:")
    for k, s in sorted(sums.items()):
        print(f"  {k}: {s:.4f}")

# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Agrupa ORFs por KEGG ID (FC + FV) y calcula TPM compartido."
    )
    p.add_argument("--fc",  required=True,  metavar="FILE",
                   help="Orftable de Fracción Celular (.prok.tsv)")
    p.add_argument("--fv",  required=True,  nargs="+", metavar="FILE",
                   help="Orftables de Fracción Vesicular (uno por estación)")
    p.add_argument("--kegg_txt", required=True, metavar="FILE",
                   help="TXT de referencia KEGG con KEGG ID, gene name, función y rutas")
    p.add_argument("--out_tpm", required=True, metavar="FILE",
                   help="Archivo de salida con Raw + TPM")
    p.add_argument("--out_tax", required=True, metavar="FILE",
                   help="Archivo de salida con KEGG_ID × Taxonomía")
    p.add_argument("--out_kegg_ref", required=True, metavar="FILE",
                   help="Archivo de salida con KEGG_ID × gene_name × función × rutas del TXT")
    return p.parse_args()


def main() -> None:
    args = parse_args()

    fc_path  = Path(args.fc)
    fv_paths = [Path(p) for p in args.fv]
    kegg_txt_path = Path(args.kegg_txt)

    for p in [fc_path, kegg_txt_path] + fv_paths:
        if not p.exists():
            sys.exit(f"[ERROR] No existe el archivo: {p}")

    # ---- Referencia KEGG desde TXT ----
    print("[INFO] Leyendo referencia KEGG desde TXT...")
    kegg_ref = parse_kegg_txt(kegg_txt_path)
    print(f"  [KEGG TXT] KEGG IDs leídos: {len(kegg_ref)}")

    # ---- Lectura ----
    print("[INFO] Leyendo FC...")
    fc_data, rpk_fc = read_fc(fc_path)

    print("[INFO] Leyendo FV...")
    fv_data, rpk_fv, fv_stations = read_fv(fv_paths)

    rpk_denom = rpk_fc + rpk_fv
    print(f"\n[INFO] Denominador TPM compartido (suma RPK total): {rpk_denom:.2f}")
    print(f"         FC: {rpk_fc:.2f}  |  FV: {rpk_fv:.2f}")

    # ---- Construcción de tablas ----
    tpm_rows = build_tpm_rows(fc_data, fv_data, rpk_denom, fv_stations, kegg_ref)
    tax_rows = build_tax_rows(fc_data, fv_data)
    kegg_ref_rows = build_kegg_ref_rows(fc_data, fv_data, kegg_ref)

    # ---- Escritura ----
    write_tsv(tpm_rows, Path(args.out_tpm), "Tabla TPM", build_tpm_fieldnames(fv_stations))
    write_tsv(tax_rows, Path(args.out_tax), "Tabla Taxonomía")
    write_tsv(kegg_ref_rows, Path(args.out_kegg_ref), "Tabla referencia KEGG TXT")

    print_tpm_summary(tpm_rows)


if __name__ == "__main__":
    main()
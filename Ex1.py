#!/usr/bin/env python3
"""
Trabajo Práctico 1 - Ejercicio 1: Procesamiento de Secuencias
Introducción a la Bioinformática - UTN

Lee una o más secuencias en formato GenBank (.gbk), traduce los 6 marcos
de lectura posibles (+1, +2, +3, -1, -2, -3), identifica los Marcos Abiertos
de Lectura (ORFs) y genera un archivo FASTA con las secuencias traducidas.
Informa en cuál de los 6 marcos está codificada la secuencia funcional.
"""

import argparse
import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="TP1 - Ejercicio 1: Traduccion a 6 marcos y deteccion de ORFs"
    )
    parser.add_argument(
        "-i", "--input",
        default="data/inputs/NM_000513.gbk",
        help="Ruta al archivo GenBank de entrada (por defecto: data/inputs/NM_000513.gbk)"
    )
    parser.add_argument(
        "-o", "--output",
        default="data/outputs/Ex1_orfs.fas",
        help="Ruta al archivo FASTA de salida con los ORFs (por defecto: data/outputs/Ex1_orfs.fas)"
    )
    parser.add_argument(
        "-m", "--min-len",
        type=int,
        default=30,
        help="Longitud minima en aminoacidos para exportar un ORF (por defecto: 30)"
    )
    parser.add_argument(
        "--best-only",
        action="store_true",
        help="Exportar unicamente el ORF funcional principal en lugar de todos los marcos"
    )
    return parser.parse_args()


def extract_annotated_cds(record):
    """Extrae informacion de la caracteristica CDS anotada en el GenBank si existe."""
    for feature in record.features:
        if feature.type == "CDS":
            protein_id = feature.qualifiers.get("protein_id", ["Desconocido"])[0]
            product = feature.qualifiers.get("product", ["Proteina no especificada"])[0]
            translation = feature.qualifiers.get("translation", [""])[0]
            location = str(feature.location)
            return {
                "exists": True,
                "protein_id": protein_id,
                "product": product,
                "translation": translation,
                "length": len(translation),
                "location": location
            }
    return {"exists": False}


def find_orfs_in_frame(nucleotide_seq, frame_name, strand, shift, min_len=30):
    """
    Traduce una secuencia ajustada al marco de lectura y localiza ORFs iniciados en Met (ATG).
    Calcula coordenadas biologicas (1-based) relativas a la hebra original.
    """
    total_nt = len(nucleotide_seq)
    sub_seq = nucleotide_seq[shift:]
    valid_len = len(sub_seq) - (len(sub_seq) % 3)
    trimmed_seq = sub_seq[:valid_len]
    translation = trimmed_seq.translate()
    trans_str = str(translation)

    orfs = []
    # Buscar cada codón de inicio 'M' seguido eventualmente de un stop '*'
    i = 0
    while i < len(trans_str):
        if trans_str[i] == "M":
            stop_pos = trans_str.find("*", i)
            if stop_pos != -1:
                peptide = trans_str[i:stop_pos]
                if len(peptide) >= min_len:
                    # Calculo de coordenadas sobre la hebra analizada
                    start_nt = shift + (i * 3) + 1
                    end_nt = shift + ((stop_pos + 1) * 3)
                    orfs.append({
                        "frame": frame_name,
                        "strand": strand,
                        "aa_start": i,
                        "aa_end": stop_pos,
                        "nt_start": start_nt,
                        "nt_end": end_nt,
                        "length": len(peptide),
                        "seq": peptide
                    })
            else:
                # ORF sin codon de parada hasta el final de la secuencia
                peptide = trans_str[i:]
                if len(peptide) >= min_len:
                    start_nt = shift + (i * 3) + 1
                    end_nt = shift + valid_len
                    orfs.append({
                        "frame": frame_name,
                        "strand": strand,
                        "aa_start": i,
                        "aa_end": len(trans_str),
                        "nt_start": start_nt,
                        "nt_end": end_nt,
                        "length": len(peptide),
                        "seq": peptide
                    })
        i += 1

    return {
        "frame": frame_name,
        "strand": strand,
        "shift": shift,
        "total_aa": len(trans_str),
        "stop_count": trans_str.count("*"),
        "translation": trans_str,
        "orfs": orfs
    }


def analyze_record(record, min_len=30):
    """Analiza los 6 marcos de lectura para un registro GenBank."""
    seq_dna = record.seq
    rc_dna = seq_dna.reverse_complement()
    frames_data = []

    # Hebra directa: Marcos +1, +2, +3
    for shift in range(3):
        frame_name = f"+{shift + 1}"
        data = find_orfs_in_frame(seq_dna, frame_name, "+", shift, min_len)
        frames_data.append(data)

    # Hebra reversa complementaria: Marcos -1, -2, -3
    for shift in range(3):
        frame_name = f"-{shift + 1}"
        data = find_orfs_in_frame(rc_dna, frame_name, "-", shift, min_len)
        frames_data.append(data)

    return frames_data


def print_report(record, frames_data, annotated_cds, best_orf):
    """Imprime el informe detallado en consola sobre los 6 marcos."""
    print("=" * 80)
    print("INFORME DE PROCESAMIENTO DE SECUENCIAS - TRABAJO PRACTICO 1")
    print("=" * 80)
    print(f"Registro:     {record.id}")
    print(f"Descripcion:  {record.description}")
    print(f"Longitud:     {len(record.seq)} nucleotidos")
    print(f"Organismo:    {record.annotations.get('organism', 'No especificado')}")

    if annotated_cds["exists"]:
        print("\nAnotacion CDS de referencia en GenBank:")
        print(f"  Producto:   {annotated_cds['product']}")
        print(f"  Proteina:   {annotated_cds['protein_id']}")
        print(f"  Ubicacion:  {annotated_cds['location']}")
        print(f"  Longitud:   {annotated_cds['length']} aminoacidos")

    print("\n" + "-" * 80)
    print("ANALISIS DE LOS 6 MARCOS DE LECTURA (READING FRAMES)")
    print("-" * 80)
    print(f"{'Marco':<8}{'Hebra':<8}{'Despl.':<8}{'Total aa':<12}{'Stops (*)':<12}{'ORFs >= min':<14}{'Mayor ORF (aa)'}")
    print("-" * 80)

    for fd in frames_data:
        max_len = max([o["length"] for o in fd["orfs"]], default=0)
        print(
            f"{fd['frame']:<8}"
            f"{fd['strand']:<8}"
            f"{fd['shift']:<8}"
            f"{fd['total_aa']:<12}"
            f"{fd['stop_count']:<12}"
            f"{len(fd['orfs']):<14}"
            f"{max_len}"
        )

    print("-" * 80)
    print("\nCONCLUSION:")
    if best_orf:
        print(f"-> El marco de lectura correcto es: MARCO {best_orf['frame']}")
        print(f"-> Cadena:                     {best_orf['strand']}")
        print(f"-> Longitud del ORF principal: {best_orf['length']} aminoacidos")
        print(f"-> Coordenadas aproximadas nt: {best_orf['nt_start']} .. {best_orf['nt_end']}")

        if annotated_cds["exists"]:
            matches = best_orf["seq"] == annotated_cds["translation"]
            print(f"-> Coincidencia con CDS NCBI:  {'EXACTA (100%)' if matches else 'PARCIAL'}")
    print("=" * 80)


def export_fasta(output_path, record, frames_data, best_orf, best_only=False):
    """Exporta los ORFs seleccionados a un archivo FASTA."""
    fasta_records = []

    if best_only and best_orf:
        header_id = f"{record.id}_ORF_principal"
        desc = (
            f"frame={best_orf['frame']} strand={best_orf['strand']} "
            f"length={best_orf['length']}aa nt={best_orf['nt_start']}..{best_orf['nt_end']}"
        )
        rec = SeqRecord(Seq(best_orf["seq"]), id=header_id, description=desc)
        fasta_records.append(rec)
    else:
        # Exportar el ORF mas largo de cada uno de los 6 marcos para permitir evaluacion comparativa
        for fd in frames_data:
            if fd["orfs"]:
                longest = max(fd["orfs"], key=lambda o: o["length"])
                is_best = (longest == best_orf)
                tag = "CODIFICANTE_REAL" if is_best else "candidato"
                header_id = f"{record.id}_frame_{fd['frame']}"
                desc = (
                    f"estado={tag} frame={fd['frame']} strand={fd['strand']} "
                    f"length={longest['length']}aa nt={longest['nt_start']}..{longest['nt_end']}"
                )
                rec = SeqRecord(Seq(longest["seq"]), id=header_id, description=desc)
                fasta_records.append(rec)

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    count = SeqIO.write(fasta_records, output_path, "fasta")
    return count


def main():
    args = parse_arguments()
    input_path = Path(args.input)
    output_path = Path(args.output)

    if not input_path.exists():
        print(f"Error: no se encontro el archivo de entrada '{input_path}'", file=sys.stderr)
        sys.exit(1)

    records = list(SeqIO.parse(input_path, "genbank"))
    if not records:
        print(f"Error: el archivo '{input_path}' no contiene registros GenBank validos", file=sys.stderr)
        sys.exit(1)

    for record in records:
        annotated_cds = extract_annotated_cds(record)
        frames_data = analyze_record(record, min_len=args.min_len)

        # Seleccionar el mejor ORF de todos los marcos
        all_orfs = [orf for fd in frames_data for orf in fd["orfs"]]
        best_orf = max(all_orfs, key=lambda o: o["length"]) if all_orfs else None

        print_report(record, frames_data, annotated_cds, best_orf)
        count = export_fasta(output_path, record, frames_data, best_orf, best_only=args.best_only)
        print(f"\nSe escribieron {count} secuencias en formato FASTA en:")
        print(f"  {output_path.resolve()}\n")


if __name__ == "__main__":
    main()

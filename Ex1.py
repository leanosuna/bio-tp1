from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

INPUT_GBK = "data/inputs/NM_000513.gbk"
OUTPUT_ORFS = "data/outputs/Ex1_orfs.fas"
OUTPUT_BEST = "data/outputs/OPN1MW_best_orf.fas"


def obtener_mayor_orf(secuencia_nt):
    largo_multiplo = len(secuencia_nt) - (len(secuencia_nt) % 3)
    traduccion = str(secuencia_nt[:largo_multiplo].translate())

    mejor_orf = ""
    for i, aa in enumerate(traduccion):
        if aa == "M":
            fin = traduccion.find("*", i)
            candidato = traduccion[i:fin] if fin != -1 else traduccion[i:]
            if len(candidato) > len(mejor_orf):
                mejor_orf = candidato
    return mejor_orf


for record in SeqIO.parse(INPUT_GBK, "genbank"):
    print(f"Secuencia: {record.id} - {record.description} ({len(record.seq)} pb)\n")

    hebras = [
        ("+", record.seq),
        ("-", record.seq.reverse_complement()),
    ]

    registros_fasta = []
    mejor_marco = None
    mejor_secuencia = ""

    for signo, seq_hebra in hebras:
        for desplazamiento in range(3):
            marco = f"{signo}{desplazamiento + 1}"
            orf = obtener_mayor_orf(seq_hebra[desplazamiento:])
            print(f"Marco {marco}: mayor ORF de {len(orf)} aminoacidos")

            registros_fasta.append(
                SeqRecord(
                    Seq(orf),
                    id=f"{record.id}_frame_{marco}",
                    description=f"frame={marco} length={len(orf)}aa",
                )
            )

            if len(orf) > len(mejor_secuencia):
                mejor_secuencia = orf
                mejor_marco = marco

    SeqIO.write(registros_fasta, OUTPUT_ORFS, "fasta")
    SeqIO.write([r for r in registros_fasta if f"frame_{mejor_marco}" in r.id], OUTPUT_BEST, "fasta")

    print(f"\nMarco de lectura codificante: {mejor_marco} ({len(mejor_secuencia)} aa)")
    print(f"Archivos guardados en {OUTPUT_ORFS} y {OUTPUT_BEST}")

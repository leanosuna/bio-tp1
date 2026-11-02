from io import StringIO
import urllib.request
from Bio import SeqIO
from pyfamsa import Aligner, Sequence

QUERY_FASTA = "data/outputs/OPN1MW_best_orf.fas"
HOMOLOGS_FASTA = "data/outputs/Ex3_secuencias.fas"
MSA_FASTA = "data/outputs/Ex3_msa.fas"
MSA_REPORT = "data/outputs/Ex3_msa.txt"

ACCESIONES_ESPECIES = [
    ("O35478", "Sciurus_carolinensis"),
    ("O18912", "Equus_caballus"),
    ("O18913", "Felis_catus"),
    ("O35599", "Mus_musculus"),
    ("P22329", "Gallus_gallus"),
]

query_record = SeqIO.read(QUERY_FASTA, "fasta")
query_record.id = "Homo_sapiens_OPN1MW"
query_record.description = "Query NM_000513.2_frame_+2"
registros = [query_record]

print("Descargando secuencias homologas desde UniProt/Swiss-Prot...")
for acc, especie in ACCESIONES_ESPECIES:
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    fasta_txt = urllib.request.urlopen(url).read().decode()
    rec = SeqIO.read(StringIO(fasta_txt), "fasta")
    rec.id = f"{especie}_{acc}"
    registros.append(rec)
    print(f"  Descargada: {rec.id} ({len(rec.seq)} aa)")

SeqIO.write(registros, HOMOLOGS_FASTA, "fasta")

secuencias_famsa = [
    Sequence(r.id.encode(), str(r.seq).encode())
    for r in registros
]

aligner = Aligner()
msa = aligner.align(secuencias_famsa)

alineadas = [
    (seq.id.decode(), seq.sequence.decode())
    for seq in msa
]

with open(MSA_FASTA, "w") as f:
    for seq_id, seq_str in alineadas:
        f.write(f">{seq_id}\n{seq_str}\n")

largo_aln = len(alineadas[0][1])
simbolos_conservacion = []
identicas = 0

for col in range(largo_aln):
    columna = [seq_str[col] for _, seq_str in alineadas]
    if "-" in columna:
        simbolos_conservacion.append(" ")
    elif len(set(columna)) == 1:
        simbolos_conservacion.append("*")
        identicas += 1
    elif len(set(columna)) == 2:
        simbolos_conservacion.append(":")
    else:
        simbolos_conservacion.append(".")

linea_cons = "".join(simbolos_conservacion)
pct_conservacion = (identicas / largo_aln) * 100

lineas_txt = [
    "ALINEAMIENTO MULTIPLE DE SECUENCIAS (MSA - FAMSA)",
    f"Secuencias alineadas: {len(alineadas)}",
    f"Longitud del alineamiento: {largo_aln} columnas",
    f"Posiciones 100% identicas (*): {identicas} ({pct_conservacion:.1f}%)\n",
]

bloque = 60
for inicio in range(0, largo_aln, bloque):
    for seq_id, seq_str in alineadas:
        fragmento = seq_str[inicio : inicio + bloque]
        lineas_txt.append(f"{seq_id:<28} {fragmento}")
    lineas_txt.append(f"{'':<28} {linea_cons[inicio : inicio + bloque]}\n")

reporte_final = "\n".join(lineas_txt)
with open(MSA_REPORT, "w") as f:
    f.write(reporte_final)

print(f"\nLongitud total del MSA: {largo_aln} posiciones")
print(f"Residuos estrictamente conservados (*): {identicas}/{largo_aln} ({pct_conservacion:.1f}%)")
print(f"Archivos guardados en {HOMOLOGS_FASTA}, {MSA_FASTA} y {MSA_REPORT}")

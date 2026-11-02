import os
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

INPUT_FASTA = "data/outputs/OPN1MW_best_orf.fas"
OUTPUT_REPORT = "data/outputs/blast.out"
OUTPUT_XML = "data/outputs/blast.xml"
BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"

if not os.path.exists(OUTPUT_XML):
    with open(INPUT_FASTA, "r") as f:
        secuencia = f.read()

    datos = urllib.parse.urlencode({
        "email": "losuna@frba.utn.edu.ar",
        "program": "blastp",
        "stype": "protein",
        "database": "uniprotkb_swissprot",
        "sequence": secuencia,
    }).encode()

    print("Enviando consulta BLASTP contra Swiss-Prot...")
    job_id = urllib.request.urlopen(f"{BASE_URL}/run", data=datos).read().decode().strip()
    print(f"Job ID: {job_id}")

    while True:
        time.sleep(2)
        estado = urllib.request.urlopen(f"{BASE_URL}/status/{job_id}").read().decode().strip()
        print(f"Estado: {estado}")
        if estado == "FINISHED":
            break
        if estado in ("ERROR", "FAILURE", "NOT_FOUND"):
            raise RuntimeError(f"Fallo en la ejecucion de BLAST: {estado}")

    reporte_txt = urllib.request.urlopen(f"{BASE_URL}/result/{job_id}/out").read().decode()
    with open(OUTPUT_REPORT, "w") as f:
        f.write(reporte_txt)

    reporte_xml = urllib.request.urlopen(f"{BASE_URL}/result/{job_id}/xml").read().decode()
    with open(OUTPUT_XML, "w") as f:
        f.write(reporte_xml)

ns = {"ebi": "http://www.ebi.ac.uk/schema"}
root = ET.parse(OUTPUT_XML).getroot()
hits = root.findall(".//ebi:hit", ns)

print(f"Total de hits en Swiss-Prot: {len(hits)}\n")
for hit in hits[:15]:
    ac = hit.attrib["ac"]
    hit_id = hit.attrib["id"]
    desc = hit.attrib["description"]
    aln = hit.find(".//ebi:alignment", ns)
    score = aln.find("ebi:score", ns).text
    bits = aln.find("ebi:bits", ns).text
    evalue = aln.find("ebi:expectation", ns).text
    identity = aln.find("ebi:identity", ns).text
    positives = aln.find("ebi:positives", ns).text

    print(
        f"{ac:<8} {hit_id:<12} | Score: {score:>4} | Bits: {bits:>6} | "
        f"E-val: {evalue:>6} | Ident: {identity:>5}% | Pos: {positives:>5}% | {desc[:42]}"
    )

print(f"\nReportes disponibles en {OUTPUT_REPORT} y {OUTPUT_XML}")

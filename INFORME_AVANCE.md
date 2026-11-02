# Informe de avance: TP 1 Parte 1 Bioinformática

| Estudiante | Legajo |
| :-- |  :--- |
| D'antoni, Priscila  | 167.275-7 |
| Da Ruos, Ailén | 169.265-3|
| Mazzeo, Melina | 146.473-5| 
| Osuna, Leandro | 167.973-9|

 
---

## 1. Selección de patología y gen en OMIM

Para el desarrollo del trabajo práctico se seleccionó el **daltonismo parcial deutan**, cuyo código es **OMIM 303800**. Esta condición hereditaria ligada al cromosoma X afecta la percepción del eje cromático rojo-verde por ausencia o alteración funcional de los conos verdes, responsables de longitudes de onda medias.

Se identificó el gen responsable **OPN1MW**, ubicado en el brazo largo del cromosoma X en la región Xq28. Desde la base de datos NCBI Nucleotide se descargó el transcripto maduro de referencia **NM_000513.2** en formato GenBank. Se eligió el ARN mensajero maduro porque ya atravesó el proceso de corte y empalme y carece de intrones, lo que permite traducir de forma continua la región codificante de la proteína. La secuencia tiene una longitud total de 1998 pares de bases.

---

## 2. Ejercicio 1: Procesamiento de secuencias y detección del ORF



Se implementó el script `Ex1.py` utilizando Biopython para leer el archivo, generar las seis traducciones completas y localizar en cada marco el mayor marco abierto de lectura, definido como la secuencia continua más larga de aminoácidos que inicia en un codón de metionina `M` y termina antes de un codón de parada `*`. El script exporta el candidato más largo de cada marco y guarda el marco funcional principal.

### Resultados obtenidos

| Marco de lectura | Hebra | Desplazamiento inicial | Longitud del mayor ORF |
| :--- | :--- | :--- | :--- |
| **+1** | Directa | 0 nucleótidos | 102 aminoácidos |
| **+2** | Directa | 1 nucleótido | **364 aminoácidos** |
| **+3** | Directa | 2 nucleótidos | 87 aminoácidos |
| **-1** | Reversa complementaria | 0 nucleótidos | 146 aminoácidos |
| **-2** | Reversa complementaria | 1 nucleótido | 104 aminoácidos |
| **-3** | Reversa complementaria | 2 nucleótidos | 41 aminoácidos |

### Interpretación y conclusión del Ejercicio 1

El marco de lectura codificante real es el **Marco +2**, que se extiende entre los nucleótidos 83 y 1177 del transcripto y genera un polipéptido continuo de **364 aminoácidos**. Los primeros 82 nucleótidos corresponden a la región 5 prima no traducida y los nucleótidos posteriores al 1177 conforman la región 3 prima no traducida. Los otros cinco marcos contienen múltiples codones de parada internos que interrumpen prematuramente cualquier traducción, dejando péptidos cortos sin función biológica. La secuencia de 364 aminoácidos obtenida en el Marco +2 coincide en un 100% con la proteína de referencia `NP_000504.1` documentada en NCBI.

---

## 3. Ejercicio 2: Búsqueda por BLASTP e interpretación estadística


Se desarrolló el script `Ex2.py` para ejecutar una búsqueda blastp remota enviando el ORF funcional de 364 aminoácidos contra la base de datos curada **UniProtKB/Swiss-Prot**. Se utilizó Swiss-Prot porque contiene secuencias proteicas revisadas manualmente y anotadas con precisión experimental. La búsqueda empleó la matriz de sustitución `BLOSUM62` con penalización de 11 por apertura de gap y 1 por extensión.

### Resultados principales en Swiss-Prot

La consulta arrojó 500 hits significativos en vertebrados. A continuación se resumen los principales resultados representativos:

| Accesión | Identificador | Organismo y proteína | Score | Bit-score | E-value | Identidad | Positivos |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **P04001** | `OPSG_HUMAN` | *Homo sapiens*, opsina onda media 1 | 1944 | 753.4 | 0.0 | 100.0% | 100.0% |
| **P04000** | `OPSR_HUMAN` | *Homo sapiens*, opsina onda larga 1 | 1887 | 731.5 | 0.0 | 96.4% | 98.1% |
| **P34989** | `OPSL_CALJA` | *Callithrix jacchus*, tití común | 1806 | 700.3 | 0.0 | 94.3% | 98.3% |
| **O35478** | `OPSG_SCICA` | *Sciurus carolinensis*, ardilla gris | 1800 | 698.0 | 0.0 | 90.4% | 95.9% |
| **O18912** | `OPSR_HORSE` | *Equus caballus*, caballo | 1799 | 697.6 | 0.0 | 90.1% | 95.3% |
| **O18913** | `OPSR_FELCA` | *Felis catus*, gato doméstico | 1799 | 697.6 | 0.0 | 90.1% | 95.3% |
| **O18914** | `OPSR_CANLF` | *Canis lupus familiaris*, perro | 1787 | 693.0 | 0.0 | 89.6% | 94.8% |
| **Q9BGI7** | `OPSR_BOVIN` | *Bos taurus*, bovino | 1785 | 692.2 | 0.0 | 89.0% | 94.5% |
| **O18910** | `OPSG_RABIT` | *Oryctolagus cuniculus*, conejo | 1756 | 681.0 | 0.0 | 87.9% | 94.0% |
| **O35599** | `OPSG_MOUSE` | *Mus musculus*, ratón | 1735 | 672.9 | 0.0 | 88.2% | 94.1% |
| **P22329** | `OPSR_CHICK` | *Gallus gallus*, gallo | 1635 | 634.0 | 0.0 | 83.2% | 90.9% |
| **Q9W6A7** | `OPSR1_DANRE` | *Danio rerio*, pez cebra | 1541 | 597.0 | 0.0 | 77.7% | 87.6% |



## Conclusión de los resultados de BLAST

El resultado de `BLASTP` confirma que la secuencia del Marco +2 es la opsina humana de ondas medias `P04001`, con un 100% de identidad, un Score de 1944, un Bit-score de 753.4 y un E-value de `0.0`. Como el E-value mide la probabilidad de encontrar un alineamiento igual por simple azar y acá vale cero, queda demostrado que las secuencias encontradas comparten un origen evolutivo real. El dato más llamativo es el hit `P04000`, que corresponde a la opsina del cono rojo humano `OPN1LW` y tiene un 96.4% de identidad con nuestra secuencia. Que ambos genes sean casi idénticos y estén pegados en el cromosoma X explica por qué ocurre el daltonismo: durante la meiosis es muy fácil que se recombinen mal o se pierda el gen `OPN1MW`. En el resto de las especies animales, la identidad baja un poco pero el porcentaje de positivos sigue por encima del 87% porque la mayoría de los cambios son sustituciones conservativas entre aminoácidos parecidos que no alteran la función de la proteína.


---

## 4. Ejercicio 3: Alineamiento múltiple de secuencias

Se implementó el script `Ex3.py` para descargar automáticamente desde UniProt las secuencias completas en formato FASTA de cinco especies distintas presentes en los resultados del BLAST:

1. *Sciurus carolinensis*, ardilla gris oriental, accesión `O35478`, 364 aminoácidos.
2. *Equus caballus*, caballo, accesión `O18912`, 364 aminoácidos.
3. *Felis catus*, gato doméstico, accesión `O18913`, 364 aminoácidos.
4. *Mus musculus*, ratón, accesión `O35599`, 359 aminoácidos.
5. *Gallus gallus*, gallo, accesión `P22329`, 362 aminoácidos.

Las cinco secuencias se combinaron con la secuencia humana `OPN1MW` en `data/outputs/Ex3_secuencias.fas` y se alinearon mediante el motor de alineamiento múltiple FAMSA, generando el archivo alineado `data/outputs/Ex3_msa.fas` y el reporte con símbolos de conservación `data/outputs/Ex3_msa.txt`.

### Alineamiento múltiple obtenido

```text
ALINEAMIENTO MULTIPLE DE SECUENCIAS (MSA - FAMSA)
Secuencias alineadas: 6
Longitud del alineamiento: 366 columnas
Posiciones 100% identicas (*): 272 (74.3%)

Homo_sapiens_OPN1MW          MAQQWSLQRLAGRHPQDSYEDSTQSSIFTYTNSNSTRGPFEGPNYHIAPRWVYHLTSVWM
Sciurus_carolinensis_O35478  MAQRWDPQRLAGGQPQDSHEDSTQSSIFTYTNSNATRGPFEGPNYHIAPRWVYHITSTWM
Equus_caballus_O18912        MAQRWGPQKLAGGQPQAGFEDSTQASIFTYTNNNATRDPFEGPNYHIAPRWVYHVTSAWM
Felis_catus_O18913           MTQRWGPQRLAGGQPHAGLEDSTRASIFTYTNSNATRGPFEGPNYHIAPRWVYHVTSAWM
Mus_musculus_O35599          -----MAQRLTGEQTLDHYEDSTHASIFTYTNSNSTKGPFEGPNYHIAPRWVYHLTSTWM
Gallus_gallus_P22329         -MAAWEAAFAARRRHEE--EDTTRDSVFTYTNSNNTRGPFEGPNYHIAPRWVYNLTSVWM
                                   ..:.:::.....  **:*..*:*****:*.*::***************:.**.**

Homo_sapiens_OPN1MW          IFVVIASVFTNGLVLAATMKFKKLRHPLNWILVNLAVADLAETVIASTISVVNQVYGYFV
Sciurus_carolinensis_O35478  IIVVIASVFTNGLVLVATMKFKKLRHPLNWILVNLAIADLAETVIASTISVVNQLYGYFV
Equus_caballus_O18912        IFVVIASVFTNGLVLAATMRFKKLRHPLNWILVNLAVADLAETIIASTISVVNQIYGYFV
Felis_catus_O18913           IFVVIASVFTNGLVLAATMKFKKLRHPLNWILVNLAVADLAETIIASTISVVNQIYGYFV
Mus_musculus_O35599          ILVVVASVFTNGLVLAATMRFKKLRHPLNWILVNLAVADLAETIIASTISVVNQIYGYFV
Gallus_gallus_P22329         IFVVAASVFTNGLVLVATWKFKKLRHPLNWILVNLAVADLGETVIASTISVINQISGYFI
                             *.**.**********:**::****************:***:**:*******:**.:***:

Homo_sapiens_OPN1MW          LGHPMCVLEGYTVSLCGITGLWSLAIISWERWMVVCKPFGNVRFDAKLAIVGIAFSWIWA
Sciurus_carolinensis_O35478  LGHPLCVVEGYTVSVCGITGLWSLAIISWERWLVVCKPFGNMRFDAKLAIVGIAFSWIWS
Equus_caballus_O18912        LGHPMCVVEGYTVSLCGITGLWSLAIISWERWMVVCKPFGNVRFDAKLAVAGIAFSWIWA
Felis_catus_O18913           LGHPMCVLEGYTVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLAIAGIAFSWIWA
Mus_musculus_O35599          LGHPLCVIEGYIVSLCGITGLWSLAIISWERWLVVCKPFGNVRFDAKLATVGIVFSWVWA
Gallus_gallus_P22329         LGHPMCVVEGYTVSACGITALWSLAIISWERWFVVCKPFGNIKFDGKLAVAGILFSWLWS
                             ****:**.***:**.****:************.********.:**:***.:**.***.*:

Homo_sapiens_OPN1MW          AVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMVTCCITPLSIIVLCYL
Sciurus_carolinensis_O35478  AVWTAPPIFGWSRYWPYGLKTSCGPDVFSGTSYPGVQSYMMVLMVTCCIIPLSIIILCYL
Equus_caballus_O18912        AVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMITCCIIPLSVIVLCYL
Felis_catus_O18913           AVWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSYPGVQSYMIVLMITCCIIPLSVIVLCYL
Mus_musculus_O35599          AIWTAPPIFGWSRYWPYGLKTSCGPDVFSGTSYPGVQSYMMVLMVTCCIFPLSIIVLCYL
Gallus_gallus_P22329         CAWTAPPIFGWSRYWPHGLKTSCGPDVFSGSSDPGVQSYMVVLMVTCCFFPLAIIILCYL
                             :.**************:*************:*:*******.***:***:.**::*:****

Homo_sapiens_OPN1MW          QVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVLAFCFCWGPYAFFACFAAANPGYPFH
Sciurus_carolinensis_O35478  QVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVFAYCLCWGPYTFFACFATANPGYAFH
Equus_caballus_O18912        QVWLAIRAVAKQQKESESTQKAEKEVTRMVMVMVFAFCLCWGPYTFFACFAAAHPGYAFH
Felis_catus_O18913           QVWLAIRAVAKQQKESESTQKAEKEVTRMVMVMIFAYCVCWGPYTFFACFAAAHPGYAFH
Mus_musculus_O35599          QVWLAIRAVAKQQKESESTQKAEKEVTRMVVVMVFAYCLCWGPYTFFACFATAHPGYAFH
Gallus_gallus_P22329         QVWLAIRAVAAQQKESESTQKAEKEVSRMVVVMIVAYCFCWGPYTFFACFAAANPGYAFH
                             **********:***************:***:**:.*:*.*****:******:*:***:**

Homo_sapiens_OPN1MW          PLMAALPAFFAKSATIYNPVIYVFMNRQFRNCILQLFGKKVDDGSELSSASKTEVSSVS-
Sciurus_carolinensis_O35478  PLVAALPAYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVDDTSELSSASKTEASSVS-
Equus_caballus_O18912        PLVAALPAYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVDDSSELSSVSKTEASSVS-
Felis_catus_O18913           PLVAALPAYFAKSATIYNPIIYVFMNRQFRNCIMQLFGKKVDDGSELSSASRTEASSVS-
Mus_musculus_O35599          PLVASLPSYFAKSATIYNPIIYVFMNRQFRNCILHLFGKKVDDSSELSSTSKTEVSSVS-
Gallus_gallus_P22329         PLAAALPAYFAKSATIYNPIIYVFMNRQFRNCILQLFGKKVDDGSEVST-SRTEVSSVSN
                             **.*:**::**********:*************::********.**:*: *:**:**** 

Homo_sapiens_OPN1MW          -SVSPA
Sciurus_carolinensis_O35478  -SVSPA
Equus_caballus_O18912        -SVSPA
Felis_catus_O18913           -SVSPA
Mus_musculus_O35599          -SVSPA
Gallus_gallus_P22329         SSVSPA
                              *****
```

## Conclusión del alineamiento múltiple

El alineamiento múltiple entre las seis especies es consistente porque muestra un 74.3% de posiciones totalmente idénticas y bloques largos sin huecos que corresponden a los dominios funcionales de la proteína. Para proteínas con la misma estructura, los pocos gaps aparecen únicamente en los extremos de la secuencia, sobre todo en los primeros aminoácidos del ratón y del gallo, mientras que la zona central donde están las hélices transmembrana y el sitio de unión al retinal se mantiene conservada sin interrupciones a lo largo de la evolución.

# Informe de avance: Trabajo práctico 1 de bioinformática

**Alumno:** Leandro Osuna  
**Correo:** losuna@frba.utn.edu.ar  
**Materia:** Introducción a la Bioinformática, UTN FRBA  

---

## 1. Selección de patología y gen en OMIM

Para el desarrollo del trabajo práctico se seleccionó el **daltonismo parcial deutan**, registrado en la base de datos Online Mendelian Inheritance in Man bajo la entrada **OMIM 303800**. Esta condición hereditaria ligada al cromosoma X afecta la percepción del eje cromático rojo-verde por ausencia o alteración funcional del fotopigmento de los conos sensibles a longitudes de onda medias, conocidos como conos verdes.

A partir de la ficha de OMIM 303800 se identificó el gen responsable **OPN1MW**, ubicado en el brazo largo del cromosoma X en la región Xq28. Desde la base de datos NCBI Nucleotide se descargó el transcripto maduro de referencia **NM_000513.2** en formato GenBank, almacenado en `data/inputs/NM_000513.gbk`. Se eligió el ARN mensajero maduro porque ya atravesó el proceso de corte y empalme y carece de intrones, lo que permite traducir de forma continua la región codificante de la proteína. La secuencia tiene una longitud total de 1998 pares de bases.

---

## 2. Ejercicio 1: Procesamiento de secuencias y detección del ORF

### Qué se hizo y por qué

Al recibir una secuencia de nucleótidos sin conocer a priori dónde comienza la traducción, existen seis marcos de lectura posibles: tres en la hebra directa desplazando el punto de partida en cero, uno o dos nucleótidos, denominados marcos +1, +2 y +3, y tres equivalentes sobre la hebra reversa complementaria, denominados marcos -1, -2 y -3.

Se implementó el script `Ex1.py` utilizando Biopython para leer el archivo `data/inputs/NM_000513.gbk`, generar las seis traducciones completas y localizar en cada marco el mayor marco abierto de lectura, definido como la secuencia continua más larga de aminoácidos que inicia en un codón de metionina `M` y termina antes de un codón de parada `*`. El script exporta el candidato más largo de cada marco en `data/outputs/Ex1_orfs.fas` y guarda el marco funcional principal en `data/outputs/OPN1MW_best_orf.fas`.

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

### Qué se hizo y por qué

Se desarrolló el script `Ex2.py` para ejecutar una búsqueda `blastp` remota enviando el ORF funcional de 364 aminoácidos contra la base de datos curada **UniProtKB/Swiss-Prot**. Se utilizó Swiss-Prot porque contiene secuencias proteicas revisadas manualmente y anotadas con precisión experimental. La búsqueda empleó la matriz de sustitución `BLOSUM62` con penalización de 11 por apertura de gap y 1 por extensión. El reporte completo se guardó en `data/outputs/blast.out` y `data/outputs/blast.xml`.

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

### Significado de los valores estadísticos

- **Score bruto.** Es el puntaje acumulado del alineamiento calculado sumando los valores que asigna la matriz de sustitución `BLOSUM62` a cada par de aminoácidos alineados y restando las penalizaciones por apertura y extensión de huecos. Depende directamente de la longitud de la secuencia y de la matriz elegida.
- **Bit-score.** Es el puntaje bruto normalizado mediante parámetros estadísticos propios de la matriz y del sistema de puntuación. Al estar expresado en unidades de información llamadas bits, independiza el resultado de la escala de la matriz y sirve para comparar directamente alineamientos hechos con distintas configuraciones.
- **E-value.** Representa el número esperado de alineamientos con un puntaje igual o mayor al obtenido que se encontrarían por simple azar en una base de datos del mismo tamaño. Cuanto más cercano a cero es el E-value, mayor es la certeza estadística de que las secuencias comparten un origen evolutivo común. Un valor de `0.0` indica un E-value inferior a $10^{-180}$, por lo que la homología es estadísticamente indiscutible.
- **Porcentaje de identidad y positivos.** La identidad indica qué fracción de las columnas alineadas presenta exactamente el mismo aminoácido en ambas secuencias. El porcentaje de positivos agrega a las coincidencias exactas aquellas sustituciones conservativas donde un aminoácido cambia por otro de propiedades fisicoquímicas equivalentes, como tamaño, carga o hidrofobicidad, lo que preserva el plegamiento tridimensional de la proteína.

### Interpretación biológica del resultado de BLAST

1. **Confirmación de la proteína de consulta.** El primer hit `P04001` corresponde a `OPSG_HUMAN` con 100.0% de identidad sobre los 364 residuos, confirmando que el Marco +2 codifica la opsina humana sensible a ondas medias.
2. **Parálogo humano y base genética del daltonismo.** El hit `P04000` corresponde a `OPSR_HUMAN`, la opsina humana sensible a ondas largas o pigmento del cono rojo, codificada por el gen vecino `OPN1LW`. Presenta un **96.4% de identidad** y un **98.1% de positivos** con `OPN1MW`, diferenciándose en apenas 13 de los 364 aminoácidos. Esta similitud extrema entre dos genes dispuestos uno al lado del otro en el cromosoma Xq28 explica por qué el daltonismo rojo-verde es tan frecuente en humanos: durante el entrecruzamiento meiótico se producen emparejamientos desiguales que eliminan el gen `OPN1MW` o generan genes híbridos entre el cono rojo y el cono verde.
3. **Conservación evolutiva en vertebrados.** Los hits siguientes corresponden a opsinas ortólogas de ondas medias y largas en primates, roedores, carnívoros, ungulados, aves, anfibios y peces teleósteos. Incluso frente a vertebrados distantes como el gallo `Gallus gallus` o el pez cebra `Danio rerio`, la similitud positiva supera el 87% y el E-value se mantiene en `0.0`.

---

## 4. Ejercicio 3: Alineamiento múltiple de secuencias

### Qué se hizo y por qué

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

### Interpretación del alineamiento múltiple

1. **Alto grado de conservación global.** De las 366 columnas del alineamiento entre cinco órdenes de mamíferos y un ave, **272 posiciones son 100% idénticas**, lo que representa un **74.3% de identidad estricta** tras más de 300 millones de años de divergencia evolutiva entre mamíferos y aves.
2. **Región N-terminal variable.** Las diferencias de longitud y la menor similitud se ubican en los primeros 25 residuos del extremo N-terminal extracelular, donde el ratón carece de los primeros 5 aminoácidos y el gallo presenta una deleción de 2 residuos. Esta cola extracelular no participa en el bolsillo de unión al cromóforo y tolera inserciones y deleciones sin afectar la función del receptor.
3. **Conservación de las siete hélices transmembrana y residuos catalíticos.** A partir del motivo `GPFEGPNYHIAPRWVY` en la posición 36, la arquitectura de siete hélices alfa transmembrana típica de los receptores acoplados a proteína G se conserva casi sin alteraciones:
   - En el tercer dominio transmembrana aparece conservado en las seis especies el bloque `WSLAIISWERW`, que contiene el motivo catalítico `ERW` responsable de activar la proteína G transducina al absorber un fotón.
   - En los bucles extracelulares se conservan las cisteínas `C126` y `C203` dentro del bloque `KTSCGPDVFS`, que forman un puente disulfuro indispensable para la estabilidad estructural de la opsina.
   - En la séptima hélice transmembrana se mantiene intacto el bloque `FAKSATIYNPVIYVFMNRQFR` en todas las especies analizadas. La lisina `K312` de ese segmento es el sitio activo exacto donde se enlaza covalentemente el cromóforo 11-cis-retinal mediante una base de Schiff protonada, mientras que el motivo `NPxxY` controla el cambio conformacional hacia el estado activo del fotorreceptor.

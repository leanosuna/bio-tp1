# Trabajo práctico 1: Introducción a la bioinformática

Repositorio para la resolución de los ejercicios del Trabajo Práctico 1 de la cátedra Introducción a la Bioinformática en la Universidad Tecnológica Nacional.

## Enfermedad y gen seleccionados

- **Patología:** Daltonismo parcial deutan, ceguera al color rojo-verde por defecto en el pigmento del cono verde.
- **Entrada OMIM:** 303800.
- **Gen:** OPN1MW, opsin 1 medium wave sensitive.
- **Transcripto de referencia:** NM_000513.2 en formato GenBank.
- **Proteína funcional:** NP_000504.1, opsina de 364 aminoácidos.

## Requisitos y configuración del entorno

El proyecto utiliza Python 3 y la librería Biopython. Debido a las políticas de entornos gestionados en Linux, se recomienda utilizar el entorno virtual local `.venv`.

1. Crear el entorno virtual si no existe:
   ```bash
   python3 -m venv .venv
   ```

2. Activar el entorno virtual:
   ```bash
   source .venv/bin/activate
   ```

3. Instalar las dependencias:
   ```bash
   pip install -r requirements.txt
   ```

## Ejercicio 1: Procesamiento de secuencias

El script `Ex1.py` toma un archivo GenBank con el ARN mensajero maduro, calcula la traducción en los seis marcos de lectura posibles, tres en la hebra directa y tres en la hebra reversa complementaria, identifica los marcos abiertos de lectura iniciados en metionina y genera un archivo FASTA con las secuencias resultantes.

### Ejecución básica

```bash
python Ex1.py
```

Por defecto lee el archivo `data/inputs/NM_000513.gbk` y escribe los marcos candidatos en `data/outputs/Ex1_orfs.fas`.

### Opciones disponibles

```bash
python Ex1.py --help
```

- `-i`, `--input`: Ruta al archivo GenBank de entrada.
- `-o`, `--output`: Ruta al archivo FASTA de salida.
- `-m`, `--min-len`: Longitud mínima en aminoácidos para considerar un ORF, valor por defecto 30.
- `--best-only`: Exporta únicamente el ORF codificante principal en vez de un candidato por marco.

### Ejemplo exportando solo el marco funcional principal

```bash
python Ex1.py --best-only -o data/outputs/OPN1MW_best_orf.fas
```

### Interpretación biológica del resultado

Al analizar el transcripto maduro de OPN1MW:

1. Los marcos +1 y +3 presentan codones de parada prematuros frecuentes, dejando péptidos máximos de 102 y 87 aminoácidos respectivamente.
2. Los tres marcos de la hebra reversa complementaria -1, -2 y -3 tampoco albergan secuencias funcionales completas.
3. El marco +2 contiene un marco abierto de lectura continuo de 364 aminoácidos entre los nucleótidos 83 y 1177, sin codones de parada internos. La traducción coincide al 100% con la proteína de referencia NP_000504.1 del NCBI.

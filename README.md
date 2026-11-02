# Trabajo práctico 1: Introducción a la bioinformática

Repositorio con la resolución completa de la Parte 1 del Trabajo Práctico de la cátedra Introducción a la Bioinformática en la Universidad Tecnológica Nacional.

El documento detallado con toda la fundamentación biológica, las tablas de resultados y las interpretaciones se encuentra en [INFORME_AVANCE.md](INFORME_AVANCE.md).

## Enfermedad y gen seleccionados

- **Patología.** Daltonismo parcial deutan, alteración de la visión cromática rojo-verde por defecto en el fotopigmento de los conos verdes.
- **Entrada OMIM.** 303800.
- **Gen.** `OPN1MW`, opsin 1 medium wave sensitive, cromosoma Xq28.
- **Transcripto de referencia.** `NM_000513.2` en formato GenBank, 1998 pares de bases.
- **Proteína funcional.** `NP_000504.1`, opsina de 364 aminoácidos.

## Configuración del entorno

### Linux y macOS

1. Crear el entorno virtual si no existe:
   ```bash
   python3 -m venv .venv
   ```

2. Activar el entorno virtual en Bash o Zsh:
   ```bash
   source .venv/bin/activate
   ```
   En terminal **fish**:
   ```fish
   source .venv/bin/activate.fish
   ```

3. Instalar las dependencias:
   ```bash
   pip install -r requirements.txt
   ```

### Windows

1. Crear el entorno virtual si no existe:
   ```cmd
   python -m venv .venv
   ```

2. Activar el entorno virtual:

   En PowerShell:
   ```powershell
   .venv\Scripts\Activate.ps1
   ```
   Si PowerShell restringe la ejecución de scripts, habilitar el permiso con el comando `Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned`.

   En el Símbolo del sistema clásico:
   ```cmd
   .venv\Scripts\activate.bat
   ```

3. Instalar las dependencias:
   ```cmd
   pip install -r requirements.txt
   ```

## Ejecución de los scripts

Todos los scripts tienen las rutas de entrada y salida definidas directamente en el código para ejecutarse sin argumentos adicionales.

### Ejercicio 1: Procesamiento de secuencias y detección del ORF

```bash
python Ex1.py
```

- **Entrada.** `data/inputs/NM_000513.gbk`.
- **Salidas.** `data/outputs/Ex1_orfs.fas` con el mayor ORF de cada uno de los seis marcos de lectura y `data/outputs/OPN1MW_best_orf.fas` con el ORF funcional del **Marco +2** de 364 aminoácidos.

### Ejercicio 2: Búsqueda de homología por BLASTP

```bash
python Ex2.py
```

- **Entrada.** `data/outputs/OPN1MW_best_orf.fas`.
- **Salidas.** `data/outputs/blast.out` en formato texto estándar de `BLASTP 2.16.0+` y `data/outputs/blast.xml` con los 500 hits encontrados en UniProtKB/Swiss-Prot.

### Ejercicio 3: Alineamiento múltiple de secuencias

```bash
python Ex3.py
```

- **Entrada.** `data/outputs/OPN1MW_best_orf.fas` más las secuencias homólogas de cinco especies descargadas desde UniProt: *Sciurus carolinensis* `O35478`, *Equus caballus* `O18912`, *Felis catus* `O18913`, *Mus musculus* `O35599` y *Gallus gallus* `P22329`.
- **Salidas.** `data/outputs/Ex3_secuencias.fas` con las seis secuencias sin alinear, `data/outputs/Ex3_msa.fas` en formato FASTA alineado y `data/outputs/Ex3_msa.txt` con el alineamiento formateado en bloques y símbolos de conservación.

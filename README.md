# Analizador de Secuencias de ADN 🧬 #

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/dardidallac/HCC/blob/main/Analisis_de_secunecias_ADN.ipynb)

Este script en lenguaje Python permite procesar archivos de secuencias biológicas en formato FASTA. Realiza un análisis que incluye la identificación de marcos abiertos de lectura (ORFs), traducción a aminoácidos, cálculo de contenido GC, clasificación de aminoácidos y detección de posibles dominios transmembrana.

### Características ###

Búsqueda de ORFs: Identificar codones de inicio (ATG) y parada (TAA, TAG, TGA) en ambas hebras (sentido y antisentido).

Traducción genética: Traduce las secuencias nucleotícida a aminoacidos utilizando la tabla de codones estándar universal.

Caracterización: Calcula el % de contenido GC.

Clasificación de los aminoacidos: Clasifica los aminoacidos de la secuencia a partir de sus características químicas (polares neutros, no polares, negativos, positivos y aromáticos).

Identificación de dominios transmembrana: Busca segmentos hidrofóbicos (no polares) de al menos 17 aminoácidos de longitud.

Visualización: Genera gráficos de la composición de aminoácidos de cada secuencia encontrada

Resulatdos: Genera reportes en .csv, archivos de secuencia en .fasta y reportes de texto .txt.

### Requisito ###

Tener instalada la librería `matplotlib` para la generación de gráficos:

```bash
pip install matplotlib
````

### Manual ###
Ejecuta el script desde la terminal indicando la ruta del archivo FASTA:

```bash
python Analisis_de_secunecias_ADN.py secuencia.fasta
````


### Resultados ###

El programa genera distintos archivos organizados por el tipo de información y la hebra analizada (**sentido** y **antisentido**). 


* `tabla_hebra_sentido.csv`: Contiene el N° de ORF, posiciones de inicio/fin, longitud del ADN, longitud de la proteína, %GC y el porcentaje de las 5 clases de aminoácidos.
* `tabla_hebra_antisentido.csv`: Misma estructura, pero para la hebra complementaria reversa.

* `sec_hebra_sentido.fasta` / `sec_hebra_antisentido.fasta`: Almacena las secuencias de nucleótidos de cada ORF y su correspondiente traducción a aminoácidos en formato FASTA.

* `transmembrana_hebra_sentido.txt` / `transmembrana_hebra_antisentido.txt`: un reporte que indica:
    * Si se encontraron fragmentos hidrofóbicos.
    * La secuencia del fragmento.
    * La posición dentro de la secuencia.

El script genera gráficos para cada ORF detectado, nombrados de la siguiente manera:
`Clases de aminoacidos - [Hebra] - Prot [N°].png`



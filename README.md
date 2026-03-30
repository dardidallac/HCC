# Analizador de Secuencias de ADN #

Este script en Python permite procesar archivos de secuencias biológicas en formato FASTA. Realiza un análisis exhaustivo que incluye la identificación de marcos abiertos de lectura (ORFs), traducción a aminoácidos, cálculo de contenido GC, clasificación de aminoácidos y detección de posibles dominios transmembrana.

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


## 📂 Estructura de Resultados ##

El programa genera un ecosistema de archivos organizados por el tipo de información y la hebra analizada (**sentido** y **antisentido**). 

### 1. Reportes Tabulares (`.csv`)
Ideales para abrir en Excel o procesar con librerías de datos como Pandas.
* `tabla_hebra_sentido.csv`: Contiene el N° de ORF, posiciones de inicio/fin, longitud del ADN, longitud de la proteína, %GC y el desglose porcentual de las 5 clases de aminoácidos.
* `tabla_hebra_antisentido.csv`: Misma estructura, pero para la hebra complementaria reversa.

### 2. Secuencias Genómicas y Proteicas (`.fasta`)
* `sec_hebra_sentido.fasta` / `sec_hebra_antisentido.fasta`: Almacena las secuencias de nucleótidos de cada ORF y su correspondiente traducción a aminoácidos en formato estándar FASTA.

### 3. Análisis de Dominios Transmembrana (`.txt`)
* `transmembrana_hebra_sentido.txt` / `transmembrana_hebra_antisentido.txt`: Reporte detallado que indica:
    * Si se encontraron fragmentos hidrofóbicos.
    * La secuencia del fragmento no polar.
    * La posición exacta (índices) dentro de la proteína.

### 4. Visualización Estadística (`.png`)
El script genera automáticamente gráficos de pastel para cada ORF detectado, nombrados siguiendo el patrón:
`Clases de aminoacidos - [Hebra] - Prot [N°].png`

| Categoría | Aminoácidos Incluidos |
| :--- | :--- |
| **Polares Neutros** | S, T, Q, N, C |
| **No Polares** | A, V, L, I, M, P, G |
| **Ácidos (Negativos)** | D, E |
| **Básicos (Positivos)** | K, R, H |
| **Aromáticos** | F, Y, W |

---

##Analizador de Secuencias de ADN##

Este script en Python permite procesar archivos de secuencias biológicas en formato FASTA. Realiza un análisis exhaustivo que incluye la identificación de marcos abiertos de lectura (ORFs), traducción a aminoácidos, cálculo de contenido GC, clasificación de aminoácidos y detección de posibles dominios transmembrana.

Características
Búsqueda de ORFs: Iidentificar codones de inicio (ATG) y parada (TAA, TAG, TGA) en ambas hebras (sentido y antisentido).

Traducción genética: Traduce las secuencias nucleotícida a aminoacidos utilizando la tabla de codones estándar universal.

Caracterización: Calcula el % de contenido GC.

Clasifica los aminoácidos en 5 categorías: polares neutros, no polares, negativos, positivos y aromáticos.

Identificación de Dominios Transmembrana: Busca segmentos hidrofóbicos (no polares) de al menos 17 aminoácidos de longitud.

Visualización: Genera gráficos de torta (pie charts) automáticos para la composición de aminoácidos de cada proteína encontrada.

Exportación de Datos: Genera reportes en .csv, archivos de secuencia en .fasta y reportes de texto .txt.

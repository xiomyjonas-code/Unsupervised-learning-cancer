# Análisis de Expresión Génica: Reducción de Dimensionalidad y Clustering ⚕️🧬 

Este proyecto implementa un flujo de trabajo de **Aprendizaje No Supervisado** (Unsupervised Learning) para analizar datos de secuenciación de ARN (RNA-Seq). El objetivo principal es evaluar la capacidad de diferentes algoritmos para distinguir entre 5 tipos distintos de cáncer basándose únicamente en la expresión de 11000 genes, sin utilizar etiquetas previas durante el entrenamiento.

## Descripción del Proyecto 📋 

Los datos genómicos suelen tener una alta dimensionalidad (miles de genes), lo que dificulta su visualización y análisis. Este script compara técnicas lineales y no lineales para reducir estas dimensiones y posteriormente aplica algoritmos de agrupamiento para validar si los patrones matemáticos coinciden con los diagnósticos biológicos reales.

### Dataset 
El conjunto de datos (`rna_cancer`) contiene muestras de pacientes diagnosticados con cinco tipos de cáncer:
* **BRCA:** Carcinoma invasivo de mama.
* **COAD:** Adenocarcinoma de colon.
* **KIRC:** Carcinoma renal de células claras.
* **LUAD:** Adenocarcinoma de pulmón.
* **PRAD:** Adenocarcinoma de próstata.

## Metodología 🛠️ 

El análisis se divide en cuatro fases principales:

### 1. Preprocesamiento
* Carga y limpieza de datos.
* Análisis de **Sparseness**: Evaluación de la calidad de los datos mediante el conteo de valores cero (genes no expresados) por muestra.

### 2. Reducción de Dimensionalidad
Se comparan cuatro algoritmos distintos para proyectar los datos de 800 dimensiones a un espacio 2D:
* **Isomap:** Mapeo isométrico que preserva la geometría geodésica (no lineal).
* **PCA (Principal Component Analysis):** Método lineal que maximiza la varianza explicada.
* **MDS (Multidimensional Scaling):** Preserva las distancias por pares (métrica Manhattan).
* **UMAP (Uniform Manifold Approximation and Projection):** Técnica de aprendizaje de variedades que preserva tanto la estructura local como la global.

### 3. Clustering (Agrupamiento)
* Se aplica el algoritmo **K-Means** sobre las coordenadas obtenidas por la proyección UMAP.
* Determinación del número óptimo de clusters ($k$) mediante el **Método del Codo (Elbow Method)**.

### 4. Validación Biológica
* Se utiliza una matriz de confusión para cruzar los clusters matemáticos generados por K-Means con las etiquetas reales de los tipos de cáncer.

## Requisitos e Instalación 📦 

El código está desarrollado en **R**. Se requiere la instalación de las siguientes librerías (el script incluye una rutina de instalación automática):

```r
# Librerías CRAN
install.packages(c("tidyverse", "uwot", "factoextra", "ggplot2", "BiocManager"))

# Librería Bioconductor
BiocManager::install("RDRToolbox")
````
## Estructura del Repositorio 📂 
`nosupervisado_analysis.R´: Código fuente completo en R.
`rna_cancer/´: carpeta con los archivos `data.csv´ (Matriz de expresión génica) y `labels.csv´ (etiquetas reales)
`plots/´: Carpeta con gráficos de los métodos de clusterización.

#### ---- 1. Preparación del Entorno y Librerías ---- ####

# Limpieza del entorno
rm(list=ls())

# Configuración del directorio de trabajo (Ajustar ruta local)
setwd("C:/Users/Tu directorio")

# 1.1. Carga de Librerías Estándar
packs <- c("ggplot2", "dplyr", "stats", "uwot", "factoextra", "BiocManager")
new_packs <- packs[!(packs %in% installed.packages()[,"Package"])]
if(length(new_packs)) install.packages(new_packs)

lapply(packs, library, character.only = TRUE)

# 1.2. Carga de Librería de Bioconductor (RDRToolbox)
# Verifica si RDRToolbox está instalado, si no, lo instala vía BiocManager
if (!requireNamespace("RDRToolbox", quietly = TRUE)) {
  BiocManager::install("RDRToolbox")
}
library(RDRToolbox)

######################################################
#       INFORMACIÓN GENERAL DEL SCRIPT               #
######################################################
# El presente código ejecuta un flujo de trabajo de aprendizaje NO supervisado.
# 1. Preprocesamiento y análisis de ceros (sparseness).
# 2. Reducción de Dimensionalidad Lineal y No Lineal:
#    - Isomap, PCA, MDS, UMAP.
# 3. Clustering (K-means) aplicado sobre la proyección UMAP.
# 4. Validación biológica cruzada con etiquetas reales.
######################################################

#### ---- 2. Carga y Preprocesamiento de Datos ---- ####

# Lectura de datos crudos
data.raw   <- read.csv('data/rna_cancer/data.csv')
labels.raw <- read.csv('data/rna_cancer/labels.csv')

# Convertir todo a numérico primero (quitando la primera columna de nombres)
# Asumiendo que la col 1 es (ID/SAMPLE) y el resto son genes
data_full <- data.frame(sapply(data.raw[, -1], as.numeric)) 

# 2.1. Análisis Exploratorio de Valores Nulos y Ceros
cat("\n--- Análisis de Calidad de Datos ---\n")
cat("¿Existen valores NA?:", anyNA(data_full), "\n")
cat("¿Existen valores 0?: ", any(data_full == 0), "\n")

# Conteo de ceros por columna
zero_counts <- colSums(data_full == 0)

# Visualización de la dispersión (Sparseness)
zero_df <- data.frame(Variable = names(zero_counts), Zeros = as.numeric(zero_counts))

ggplot(zero_df, aes(x = Variable, y = Zeros, fill = Variable)) +
  geom_bar(stat = "identity") +
  labs(title = "Distribución de valores cero por Gen (Sparsity)",
       x = "Genes",
       y = "Conteo de Ceros") +
  theme_minimal() +
  theme(legend.position = "none", axis.text.x = element_blank()) # Se ocultan etiquetas eje X por legibilidad

# FILTRADO POR VARIANZA 
# Calculamos la varianza de cada columna (cada gen)
varianzas <- apply(data_full, 2, var)

# Ordenar los genes de mayor a menor varianza y usamos los índices
# Se puede probar con 5000, 11000 (considerando que el data.raw contiene más de 20000 genes)  
top_genes_index <- order(varianzas, decreasing = TRUE)[1:11000] 

# Crear el dataset final solo con esos genes
data_num <- data_full[, top_genes_index]
cat("Dimensiones finales tras filtrado:", dim(data_num)) 

#cat("¿Existen valores NA?:", anyNA(data_num), "\n")
#cat("¿Existen valores 0?: ", any(data_num == 0), "\n")

#### ---- 3. Reducción de Dimensionalidad ---- ####

# =========================================================
# MÉTODO 1: Isomap (Isometric Mapping)
# =========================================================
# Algoritmo no lineal que conserva la geometría geodésica global.

print("Ejecutando Isomap...")
# Se calculan dimensiones 1 a 15 con k=8 vecinos más cercanos.
isomap.results <- Isomap(data=as.matrix(data_num), dims=1:15, k=8, plotResiduals=TRUE)

# Visualización (Dimensiones 1 y 2)
isomap.df <- data.frame(isomap.results$dim2) 
ggplot(isomap.df, aes(x = X1, y = X2, color = labels.raw$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "Isomap (k=8)", x = 'Dim 1', y = 'Dim 2', color = "Clase Real") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))


# =========================================================
# MÉTODO 2: PCA (Principal Component Analysis)
# =========================================================
# Algoritmo lineal que maximiza la varianza explicada.

print("Ejecutando PCA...")
pca.results <- prcomp(data_num, center=TRUE, scale.=FALSE) 

# Cálculo de varianza explicada acumulada
varianzas <- pca.results$sdev^2
varianza.explicada <- varianzas/sum(varianzas)
varianza.acumulada <- cumsum(varianza.explicada)

# Selección de componentes para el 90% de varianza
n.pc <- min(which(varianza.acumulada > 0.90))
cat("Componentes necesarios para explicar el 90% de la varianza:", n.pc, "\n")

# Visualización
pca.df <- data.frame(pca.results$x)
x_label <- paste0('PC1 (', round(varianza.explicada[1] * 100, 2), '%)')
y_label <- paste0('PC2 (', round(varianza.explicada[2] * 100, 2), '%)')

ggplot(pca.df, aes(x=PC1, y=PC2, color=labels.raw$Class)) +
  geom_point(size=3) +
  scale_color_manual(values=c('red', 'blue', 'green', 'orange', 'purple')) +
  labs(title='PCA - Componentes Principales', x=x_label, y=y_label, color='Clase Real') +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))


# =========================================================
# MÉTODO 3: MDS (Multidimensional Scaling)
# =========================================================
# Busca preservar las distancias por pares (Manhattan) en un espacio de baja dimensión.

print("Ejecutando MDS...")
distances <- dist(data_num, method = 'manhattan')
mds.results <- cmdscale(distances, eig=TRUE, k=2, x.ret=TRUE)

mds.df <- data.frame(mds.results$points)

ggplot(mds.df, aes(x=X1, y=X2, color=labels.raw$Class)) +
  geom_point(size=3) + 
  scale_color_manual(values=c("red", "blue", "green", "orange", "purple")) +
  labs(title="MDS (Distancia Manhattan)", x="Dim 1", y="Dim 2", color = "Clase Real") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))


# =========================================================
# MÉTODO 4: UMAP (Uniform Manifold Approximation and Projection)
# =========================================================
# Técnica moderna de aprendizaje de variedades. Preserva estructura local y global.

print("Ejecutando UMAP...")
# Configuración: 30 vecinos, distancia mínima 0.5 para clusters compactos, métrica Manhattan
umap.results <- umap(data_num, 
                     n_neighbors=30,
                     n_components = 2, 
                     min_dist = 0.5, 
                     metric ="manhattan", 
                     local_connectivity=1, 
                     ret_model = TRUE, 
                     verbose = FALSE) # Verbose False para limpiar consola

umap.df <- data.frame(umap.results$embedding)

ggplot(umap.df, aes(x = X1, y = X2, color = labels.raw$Class)) +
  geom_point(size = 3) +
  scale_color_manual(values = c("red", "blue", "green", "orange", "purple")) +
  labs(title = "UMAP Projection", x = "UMAP 1", y = "UMAP 2", color = "Clase Real") +
  theme_classic() +
  theme(plot.title=element_text(hjust=0.5))


#### ---- 4. Clustering (K-Means sobre UMAP) ---- ####

# Se aplica K-Means sobre las coordenadas generadas por UMAP para validar agrupación automática.

# 4.1. Determinación de K óptimo (Método del Codo - WSS)
set.seed(1234) # Semilla para reproducibilidad del clustering
fviz_nbclust(umap.df, kmeans, method = "wss", k.max = 10, nstart = 50) +
  geom_vline(xintercept = 5, linetype = 2) +
  ggtitle("Método del Codo (Optimal k)") +
  theme_classic()

# 4.2. Ejecución de K-Means (k=5)
# Se selecciona k=5 basado en la gráfica anterior y el conocimiento biológico (5 tipos de cáncer)
kmeans.result <- kmeans(umap.df, centers = 5, iter.max = 100, nstart = 25)

# Visualización de Clusters Calculados
fviz_cluster(kmeans.result, umap.df, 
             geom = "point", 
             ellipse.type = "convex", 
             ggtheme = theme_minimal(),
             main = "Clustering K-Means sobre proyección UMAP") +
  theme(plot.title = element_text(hjust = 0.5))


#### ---- 5. Validación Biológica (Matriz de Confusión) ---- ####

# Comparación cruzada: Clusters Matemáticos (K-means) vs Etiquetas Reales (Biología)
tabla_validacion <- table(Tipo_Real = labels.raw$Class, Cluster_Kmeans = kmeans.result$cluster)

cat("\n--- Validación Cruzada: Clusters vs Tipos de Cáncer ---\n")
print(tabla_validacion)

# INTERPRETACIÓN DE RESULTADOS:
# Se observa una alta correlación entre los clusters generados y los tipos biológicos.
# - BRCA, COAD y PRAD muestran una separación perfecta (100% pureza en cluster).
# - KIRC y LUAD presentan una separación > 98%, con mínimas discrepancias.
# Esto valida la capacidad de UMAP combinado con K-means para distinguir fenotipos de cáncer.

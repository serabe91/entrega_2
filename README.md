# Análisis de Aberraciones Cromosómicas

Este proyecto está diseñado para realizar un análisis exploratorio y visualización de datos sobre anormalidades cromosómicas. Utiliza Python y varias bibliotecas como `pandas`, `matplotlib`, `seaborn` y `reportlab` para cargar, procesar, analizar y visualizar datos relacionados con aberraciones cromosómicas.

## Contenido

1. **Carga y Procesamiento de Datos**
   - Carga archivos CSV y extrae información relevante.
   - Organiza y limpia los datos, eliminando columnas innecesarias y filtrando información específica.

2. **Análisis de Valores Atípicos**
   - Identifica valores atípicos utilizando el método del Rango Intercuartílico (IQR).
   - Genera gráficos de dispersión para comparar los datos originales con los valores atípicos.

3. **Generación de Gráficos**
   - Crea múltiples visualizaciones, incluyendo gráficos de barras, mapas de calor y gráficos de dispersión, para analizar la distribución de anormalidades cromosómicas.

4. **Generación de Reportes en PDF**
   - Crea un informe en PDF que incluye un resumen del análisis, gráficos generados y un índice de anexos.

## Funciones Principales

- `cargar_csv_y_extraer_nombre(file_path)`: Carga un archivo CSV y extrae el nombre base del archivo.
- `organizar_dataframe(df)`: Limpia y organiza el DataFrame, creando un diccionario de morfologías.
- `analizar_outliers_y_transformacion(df_sorted)`: Analiza y visualiza valores atípicos en el DataFrame.
- `generar_graficos(df)`: Genera varios tipos de gráficos para visualizar los datos.
- `crear_reporte_pdf(nombre_archivo, df, df_sorted, path_graficas)`: Genera un informe en PDF con el análisis y visualizaciones.

## Requisitos

- Python 3.12
- Bibliotecas: `pandas`, `matplotlib`, `seaborn`, `reportlab`

## Uso

1. Se deben tener instaladas las bibliotecas necesarias.
2. Colocar el archivo `total_cancer.csv` en el directorio de trabajo.
3. Ejecutar el script para realizar el análisis y generar el reporte.

## Licencia

Proyecto bajo la Licencia MIT. Consultar el archivo `LICENSE` para más detalles.

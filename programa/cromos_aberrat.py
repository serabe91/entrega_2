#===================================================================================================================#
#====================================================LIBRERIAS======================================================#
#===================================================================================================================#

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
#===================================================================================================================#
from reportlab.pdfgen import canvas
from reportlab.lib.units import inch
from reportlab.platypus import Frame
from reportlab.platypus import Image
from reportlab.platypus import Paragraph
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.styles import ParagraphStyle

#===================================================================================================================#
#====================================================FUNCIONES======================================================#
#===================================================================================================================#

def cargar_csv_y_extraer_nombre(file_path):
    """
    Esta función lee un archivo CSV desde la ruta proporcionada y lo carga en un
    DataFrame de pandas, extrayendo además el nombre base del archivo sin su extensión.
    Adicionalmente, aísla la primera parte del nombre del archivo antes del carácter
    guion bajo. Permite tanto la extracción de datos como el procesamiento básico del
    nombre del archivo en una sola función.

    Parámetros:
        file_path (str): La ruta completa del archivo CSV que se cargará y procesará.

    Retorna:
        Una tupla que contiene:
            - pandas.DataFrame: Los datos cargados desde el archivo CSV.
            - str: El nombre base extraído del archivo antes del carácter guion bajo.
"""

    df = pd.read_csv(file_path, sep = ",")    # Cargar los archivos CSV en un DataFrame
    file_name = os.path.splitext(os.path.basename(file_path))[0]    # Obtener el nombre base del archivo sin la ext.
    file_name = file_name.split("_")[0]     # Obtener solo la primera palabra antes del guion bajo
    return df, file_name     # Retornar el DataFrame y el nombre extraído

file_path = "total_cancer.csv"  # Especificar la ruta del archivo
df, nombre_archivo = cargar_csv_y_extraer_nombre(file_path)

#===================================================================================================================#

def organizar_dataframe(df):
    """
    Organiza y procesa un DataFrame de Pandas para limpiar, transformar y extraer
    información específica sobre anormalidades, morfologías y datos cromosómicos,
    mientras crea un diccionario basado en códigos morfológicos.

    Argumentos:
        df (pandas.DataFrame): DataFrame de entrada que contiene las columnas
        'Abnormality', 'Topo', 'Gene', 'Chromosome', 'Morph', 'MorphName', 'Type' y
        'TopoName'. Se espera que el DataFrame tenga datos relevantes en estos campos
        para que la función realice las transformaciones necesarias.

    Retorna:
        Una tupla que contiene:
            - pandas.DataFrame: Un DataFrame limpio y ordenado donde se eliminan columnas
            innecesarias, se procesa la información morfológica, y la columna 'Chromosome'
            es numérica y compatible para análisis y ordenamientos adicionales.
            - dict: Un diccionario donde las claves son códigos morfológicos únicos
            ('Morph_code') y los valores son las descripciones morfológicas correspondientes
            ('morfo'), con duplicados eliminados y claves ordenadas en orden numérico.
    """

    # Organizar anormalidades solo para quedar con el texto indicativo
    df["Abnormality"] = (
        df["Abnormality"]  # Selección de la columna
        .str.split("(")  # Dividir usando el paréntesis como separador
        .str[0])  # Tomar la primera parte (antes del primer paréntesis)
    # Eliminación de las columnas de topografía y genes
    df = df.drop("Topo", axis = 1)
    df = df.drop("Gene", axis = 1)
    # Filtrar filas donde la columna "Chromosome" no sea "X" o "Y"
        # isin([]) verifica si los valores en la columna "Chromosome" son "X" o "Y".
        # ~ invierte la condición, de modo que selecciona las filas donde "Chromosome" no es "X" ni "Y".
        # El DataFrame resultante deja por fuera las filas que contienen "X" o "Y" en la columna "Chromosome".
    df = df[~df["Chromosome"].isin(["X", "Y"])]
    # Asegura que las columnas son strings y reemplazar NaN con cadenas vacías
    df["Morph"] = df["Morph"].fillna("").astype(str)
    df["MorphName"] = df["MorphName"].fillna("").astype(str)
    # Concatenar las columnas y eliminar las innecesarias
    df["Morphology"] = df["Morph"] + ": " + df["MorphName"]
    df = df.drop("Morph", axis = 1)
    # Función para dividir la columna en número y morfo
    def extraer_num_morfo(texto):
        numero, morfo = texto.split(": ")
        return int(numero), morfo.strip()
    # Aplicar la función y crear una nueva columna "numero" y "morfo"
        # La función apply aplica la función a cada elemento de la columna "Morphology"
        # La función extraer_num_morfo(x) divide la columna en número y morfo
        # pd.Series() convierte los valores en una serie de Pandas
    df[["Morph_code", "morfo"]] = df["Morphology"].apply(lambda x: pd.Series(extraer_num_morfo(x)))
    # Eliminar las columnas innecesarias del DataFrame
    df = df.drop(columns = "Morphology")
    df = df.drop(columns = "MorphName")
    # Crear el diccionario usando "numero" como clave y "morfo" como valor, eliminando duplicados
    diccionario_morfo = df.drop_duplicates("Morph_code").set_index("Morph_code")["morfo"].to_dict()
    # Ordenar el diccionario basado en orden numérico de key
    diccionario_morfo = dict(sorted(diccionario_morfo.items()))
    # Reemplazar los valores "U" y "B" con su palabra correspondiente
    df["Type"] = df["Type"].replace({"U": "unbalanced", "B": "balanced"})
    # Eliminar columnas no necesarias
    df = df.drop(columns = "morfo")
    df = df.drop(columns = "TopoName")
    # Asegurarse de que "Chromosome" sea numérico
        # La función pd.to_numeric() intenta convertir la columna "Chromosome" en valores numéricos.
        # Si no puede convertir algún valor (por ejemplo, texto no numérico), lo reemplaza con NaN gracias al parámetro errors="coerce"
    df["Chromosome"] = pd.to_numeric(df["Chromosome"], errors = "coerce")
    # Ordenar después de la conversión
    df_sorted = df.sort_values(by = "Chromosome")
    df.rename(columns={"": "Morphology"}, inplace = True)
    return df_sorted, diccionario_morfo

df_sorted, dict_morfo = organizar_dataframe(df)

#===================================================================================================================#

def analizar_outliers_y_transformacion(df_sorted):
    """
    Analiza valores atípicos en un DataFrame ordenado por 'CaseCount', los identifica
    usando el método IQR, los compara con los datos originales mediante gráficos de
    dispersión y genera representaciones visuales de los valores atípicos por cromosoma.
    Devuelve el DataFrame modificado junto con los datos de los valores atípicos por separado.

    Argumentos:
        df_sorted (pandas.DataFrame): Un DataFrame ordenado por la columna 'CaseCount', que
        contiene las columnas 'CaseCount' y 'Chromosome'.

    Retorna:
        Una tupla que contiene:
        - pandas.DataFrame: El DataFrame de entrada con las modificaciones opcionales aplicadas.
        - pandas.DataFrame: Un DataFrame que contiene los valores atípicos detectados.
    """

    # IQR para CaseCount
    # Calcular Q1, Q3 e IQR
    Q1 = df_sorted["CaseCount"].quantile(0.25)
    Q3 = df_sorted["CaseCount"].quantile(0.75)
    IQR = Q3 - Q1
    # Definir límites de outliers
    lower_limit = Q1 - 1.5 * IQR
    upper_limit = Q3 + 1.5 * IQR
    # Filtrar outliers y generar tabla
    outliers_iqr = df_sorted[(df_sorted["CaseCount"] # DataFrame ordenado por la columna CaseCount
                              < lower_limit) | (df_sorted["CaseCount"] > upper_limit)] # Detección de valores menones y mayores a los límites

    # Comparación outliers vs datos originales en scatterplot
    plt.figure(figsize = (10, 6)) # Configura el tamaño de la figura
    sns.set_style("whitegrid") # Establece el estilo del gráfico
    sns.scatterplot(x = "Chromosome", # Fuente datos en x
                    y = "CaseCount", # Fuente datos en y
                    data = outliers_iqr, # DF fuente de los datos
                    color = "red", # Color de los puntos
                    label = "Outliers") # Etiqueta de los datos
    sns.scatterplot(x = "Chromosome", # Fuente datos en x
                    y = "CaseCount", # Fuente datos en y
                    data = df_sorted, # DF fuente de los datos
                    alpha = 0.5, # Transparencia de los puntos
                    label = "Datos originales") # Etiqueta de los datos
    plt.title("Comparación de Outliers vs. Datos Originales") # Asignación de título del gráfico
    plt.legend() # Recuadro de etiquetas
    plt.savefig("grafica_1.png", dpi=800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente



    # Conteo de outliers por Chromosome
    outliers_by_chromosome = outliers_iqr["Chromosome"].value_counts() # Cuenta cuántas veces aparece cada valor único en la columna Chromosome
    # Gráfico
    plt.figure(figsize = (10, 6)) # Configura el tamaño de la figura
    sns.countplot(x = "Chromosome", # Fuente datos
                  data = outliers_iqr, # DF fuente de los datos
                  palette = "viridis", # Selección de la paleta de colores
                  order = outliers_by_chromosome.index) # Establece el orden de los elementos del gráfico
    plt.title("Distribución de Outliers por Chromosome") # Asignación de título del gráfico
    plt.savefig("grafica_2.png", dpi = 800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente

    # Devolver el DataFrame modificado con la nueva columna
    return df_sorted, outliers_iqr

# Llamado de la función
df_modificado, outliers_iqr = analizar_outliers_y_transformacion(df)

#===================================================================================================================#
#=====================================================GRAFICAS======================================================#
#===================================================================================================================#

def generar_graficos(df):
    """
    Función con la que generamos varios tipos de gráficos y visualizaciones a partir de un DataFrame dado.
    Esta función genera gráficos de barras, dispersión, mapas de calor y gráficos de tiras para analizar anormalidades cromosómicas,
    su frecuencia, tipo y distribución basándose en el conjunto de datos proporcionado.

    Parámetros:
        df (pd.DataFrame): Los datos de entrada que contienen información sobre anormalidades cromosómicas.
        Se espera que el DataFrame incluya las siguientes columnas:
        - Chromosome: El identificador del cromosoma.
        - Arm: El brazo específico del cromosoma.
        - Abnormality: Tipo de anormalidad cromosómica.
        - CaseCount: La cantidad de casos para cada anormalidad.
        - Type: El tipo de anormalidad cromosómica.
        - Band: El identificador de la banda cromosómica.

    Efectos secundarios:
        - Guarda varias imágenes de gráficos (grafica_3.png, grafica_4.png, grafica_5.png, grafica_6.png, grafica_7.png,
         grafica_8.png) en el directorio de trabajo actual. Cada imagen corresponde a una visualización específica generada.

    Descripción de gráficos generados:
        1. Un gráfico de barras que muestra la distribución de anormalidades por cromosoma y brazo (grafica_3.png).
        2. Un gráfico de barras que presenta la frecuencia de anormalidades ordenadas por la cantidad de casos (grafica_4.png).
        3. Un mapa de calor que visualiza los tipos de anormalidades en los cromosomas (grafica_5.png).
        4. Un gráfico de dispersión que muestra la distribución de casos por banda cromosómica (grafica_6.png).
        5. Un gráfico de tiras que ilustra los casos por tipo de anormalidad (grafica_7.png).
        6. Un mapa de calor con indexación avanzada, que revela las frecuencias de anormalidades por banda cromosómica y tipo (grafica_8.png).

    Advertencias:
        - Nos debemos asegurar de que el DataFrame de entrada contenga datos válidos y completos en el formato requerido para una correcta generación de gráficos.
        - Los nombres de archivo generados sobrescribirán archivos existentes con el mismo nombre en el directorio.
"""

    # 1. Conteo de anormalidades por cromosoma y brazo
        # groupby agrupa los datos por las columnas "Chromosome" y "Arm".
        # [].count() cuenta cuántas entradas (no nulas) hay en la columna "Abnormality" para cada grupo.
        # reset_index() restablece el índice del DataFrame, creando un formato más limpio y accesible para trabajar con los resultados.
    resumen = df.groupby(["Chromosome", "Arm"])["Abnormality"].count().reset_index()
    # Gráfico respectivo:
    sns.barplot(x = "Chromosome", # Fuente datos en x
                y = "Abnormality", # Fuente datos en y
                hue = "Arm", # Para distinguir la variable categórica en el gráfico
                data = resumen) # DF fuente de los datos
    sns.color_palette("tab10") # Selección de la paleta de colores
    plt.grid(True) # Activa la cuadrícula del gráfico
    plt.title("Distribución de anormalidades por cromosoma y brazo") # Asignación de título del gráfico
    plt.savefig("grafica_3.png", dpi = 800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente



    # 2. Agrupar por anormalidad y sumar los casos
        # groupby() agrupa los datos por la columna "Abnormality".
        # [].sum() suma los valores de la columna "CaseCount" dentro de cada grupo de anomalías.
        # reset_index() restablece el índice, convirtiendo la columna "Abnormality" en una columna normal del DataFrame
    freq_abnormalidades = df.groupby("Abnormality")["CaseCount"].sum().reset_index()
    # Ordenar y graficar
    sns.barplot(x = "CaseCount", # Fuente datos en x
                y = "Abnormality", # Fuente datos en y
                data = freq_abnormalidades.sort_values("CaseCount", ascending = False), # sort_values() ordena el DataFrame
                                                                                        # freq_abnormalidades según la columna "CaseCount",
                                                                                        # de mayor a menor.
                color = "red") # Color de las barras
    plt.title("Frecuencia de anormalidades") # Asignación de título del gráfico
    plt.grid(True) # Activa la cuadrícula del gráfico
    plt.savefig("grafica_4.png", dpi = 800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente



    # 3. Tabla cruzada
    tabla_tipo = pd.crosstab(df["Chromosome"], df["Type"]) # Tabla de contingencia entre dos columnas del DataFrame.
    # Gráfico de calor
    sns.heatmap(tabla_tipo,
                cmap = "YlOrBr", # define el esquema de colores
                annot = True, # muestra los valores numéricos en las celdas.
                fmt = ".0f")  # fmt=".0f" muestra números enteros
    plt.title("Tipo de anormalidad por cromosoma") # Asignación de título del gráfico
    plt.savefig("grafica_5.png", dpi = 800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente



    # 4. Crear un scatterplot
    plt.figure(figsize = (12, 6)) # Configura el tamaño de la figura
    sns.scatterplot(x = "Band", # Fuente datos en x
                    y = "CaseCount", # Fuente datos en y
                    hue = "Type", # Para distinguir la variable categórica en el gráfico
                    alpha = 0.5, # Transparencia de los puntos
                    data = df, # DF fuente de los datos
                    palette = "coolwarm", # Selección de la paleta de colores
                    s = 200) # Ajusta el tamaño de los puntos en el gráfico
    plt.title("Distribución de casos por banda cromosómica") # Asignación de título del gráfico
    plt.xticks(rotation = 90) #  Rota las etiquetas del eje x 90 grados en sentido antihorario
    plt.savefig("grafica_6.png", dpi = 800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente



    # 5. Stripplot para casos por abnormality
    plt.figure(figsize = (10, 6)) # Configura el tamaño de la figura
    sns.stripplot(x = "Abnormality", # Fuente datos en x
                  y = "CaseCount", # Fuente datos en y
                  hue = "Type", # Para distinguir la variable categórica en el gráfico
                  data = df, # DF fuente de los datos
                  jitter = True, # Desplaza aleatoriamente los puntos para evitar superposiciones
                  alpha = 0.2, # Transparencia de los puntos
                  dodge = True, # Desplaza las categorías de los puntos para que no se superpongan en el gráfico
                  palette = "deep", # Selección de la paleta de colores
                  s = 10) # Ajusta el tamaño de los puntos en el gráfico
    plt.title("Casos por tipo de anormalidad") # Asignación de título del gráfico
    plt.xticks(rotation = 45) #  Rota las etiquetas del eje x 45 grados en sentido antihorario
    plt.savefig("grafica_7.png", dpi = 800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente


    # 6. Crear tabla con combinaciones
        # index = [df["Chromosome"], df["Band"]]: Define el índice de la tabla como una combinación de las columnas Chromosome y Band.
        # Creando un índice jerárquico o multiíndice, donde cada fila está identificada por un par de valores (un cromosoma y una banda)

        # columns = df["Type"]: Define las columnas de la tabla basándose en la variable Type.
        # Cada categoría de Type se convierte en una columna separada
    heatmap_avanzado = pd.crosstab(index = [df["Chromosome"], df["Band"]], columns = df["Type"])
    # Gráfico heatmap
    plt.figure(figsize = (12, 10)) # Configura el tamaño de la figura
    sns.heatmap(heatmap_avanzado,
                cmap = "viridis", # Define el esquema de colores
                annot = False, # Indica que no se muestren anotaciones numéricas en las celdas
                cbar_kws = {"label": "Frecuencia"}) # Configura las opciones del color bar que aparece al lado del mapa de calor.
    plt.title("Frecuencia de anormalidades por banda y tipo") # Asignación de título del gráfico
    plt.savefig("grafica_8.png", dpi = 800) # Guardado de gráfica generada
    plt.close() # Cerrar personalización del gráfico correspondiente


# Llamar a la función con el DataFrame df
generar_graficos(df_sorted)


#===================================================================================================================#
#===================================================================================================================#
#=================================================REPORTE PDF=======================================================#
#===================================================================================================================#
#===================================================================================================================#

# Función para crear el PDF
def crear_reporte_pdf(nombre_archivo, df, df_sorted, path_graficas):
    """
    Genera un informe en PDF con información y visualizaciones relacionadas con el análisis de aberraciones cromosómicas.
    La función utiliza la biblioteca ReportLab para crear un documento PDF de varias páginas que incluye elementos como:
    - Página de título
    - Sección de introducción
    - Tabla de contenidos
    - Índice de anexos
    - Índice de gráficos
    - Gráficos de visualización

    Argumentos:
        nombre_archivo (str): El nombre del archivo PDF resultante que se creará.
        df (DataFrame): El DataFrame que contiene información del conjunto de datos original.
        df_sorted (DataFrame): Datos preprocesados y organizados para análisis y ordenamiento.
        path_graficas (str): Ruta del directorio donde se almacenan las imágenes de gráficos de visualización.
    """

    # Crear un objeto canvas con tamaño de página personalizado (horizontal)
    width, height = 11 * inch, 8.5 * inch  # Tamaño de página en formato horizontal
    c = canvas.Canvas(nombre_archivo, pagesize = (width, height))  # Tamaño de página personalizado


    # 1. Página de presentación
    c.setFont("Helvetica-Bold", 16)
    titulo = "Reporte de Análisis de Datos de Aberraciones Cromosómicas"
    autor = "Sergio Daniel Abello Hernández"
    universidad = "Universidad Nacional de Colombia Sede Medellín"
    curso = "Fundamentos de Programación para Ciencias Biológicas"
    semestre = "2024-2S"
    # Cálculos para centrar texto horizontalmente
    def centrar_texto(texto, font, size, y):
        c.setFont(font, size)
        ancho_texto = c.stringWidth(texto, font, size)
        x = (width - ancho_texto) / 2  # Centrado horizontal
        c.drawString(x, y, texto)
    # Coordenadas de inicio (centrar verticalmente)
    espacio_entre_lineas = 115
    y_inicial = height / 2 + 2 * espacio_entre_lineas
    # Dibujar el contenido centrado
    centrar_texto(titulo, "Helvetica-Bold", 16, y_inicial)
    centrar_texto(autor, "Helvetica", 12, y_inicial - espacio_entre_lineas)
    centrar_texto(universidad, "Helvetica", 12, y_inicial - 2 * espacio_entre_lineas)
    centrar_texto(curso, "Helvetica", 12, y_inicial - 3 * espacio_entre_lineas)
    centrar_texto(semestre, "Helvetica", 12, y_inicial - 4 * espacio_entre_lineas)
    # Finalizar la página de presentación
    c.showPage()

    # 2. Página de Introducción
    # Definir estilos de texto
    styles = getSampleStyleSheet()
    title_style = styles["Title"]
    heading_style = styles["Heading2"]
    # Crear estilo para texto justificado
    justified_body_style = ParagraphStyle(
        "JustifiedBody",
        parent = styles["BodyText"],
        alignment = TA_JUSTIFY,  # Justificar texto
        fontSize = 10,
        leading = 14  # Espaciado entre líneas
    )
    # Crear el contenido como párrafos
    contenido = [
        Paragraph("<b>Introducción</b>", title_style),
        Paragraph("<b>Propósito general del código:</b>", heading_style),
        Paragraph(
            "El código está diseñado para realizar análisis exploratorio y visualización de datos sobre anormalidades cromosómicas, "
            "organizando, limpiando y representando gráficamente la información en múltiples formas. Combina procesamiento de datos, "
            "detección de valores atípicos y generación de gráficos para identificar patrones, tendencias y distribuciones relacionadas "
            "con anomalías genéticas, cromosómicas y otros atributos asociados.", justified_body_style),
        Paragraph("<b>Propósitos específicos:</b>", heading_style),
        Paragraph(
            "<b>Procesamiento de datos genéticos:</b> Se enfoca en limpiar y organizar datos cromosómicos y morfológicos para facilitar su análisis. "
            "Y filtra información específica (por ej: elimina cromosomas sexuales X e Y, organiza anormalidades y códigos morfológicos).",
            justified_body_style),
        Paragraph(
            "<b>Análisis Estadístico y Exploratorio:</b> Identifica valores atípicos (outliers) en las cuentas de casos utilizando el IQR (Rango Intercuartílico), "
            "aplica transformaciones logarítmicas para normalizar datos y facilitar su interpretación visual, y permite detectar qué cromosomas o bandas tienen "
            "mayores anomalías o concentraciones de casos.", justified_body_style),
        Paragraph(
            "<b>Visualización de Datos:</b> Genera gráficos para comprender la distribución de anormalidades y otros parámetros genéticos, como lo son: "
            "Gráficos de barras, Mapas de calor, Scatterplots y FacetGrid, facilitando el análisis visual de los datos cromosómicos.",
            justified_body_style),
        Paragraph(
            "<b>Creación de un Diccionario de Morfologías:</b> Extrae y organiza información de códigos morfológicos y nombres de anomalías, "
            "creando un diccionario ordenado que puede servir como referencia para análisis posteriores.",
            justified_body_style),
        Paragraph("<b>Fuente de información:</b>", heading_style),
        Paragraph(
            "Toda la información fue extraída de: Mitelman Database of Chromosome Aberrations and Gene Fusions in Cancer."
            " La información de la base de datos Mitelman sobre aberraciones cromosómicas y fusiones de genes en el cáncer relaciona los cambios citogenéticos"
            " y sus consecuencias genómicas, en particular las fusiones de genes, con las características de los tumores, ya sea en casos individuales o en asociaciones."
            " Todos los datos han sido seleccionados manualmente de la literatura por Felix Mitelman en colaboración con Bertil Johansson y Fredrik Mertens. "
            "La base de datos Mitelman cuenta con el apoyo del Instituto Nacional del Cáncer, la Sociedad Sueca contra el Cáncer y la Fundación Sueca contra el "
            "Cáncer Infantil. La base de datos se actualiza trimestralmente en enero, abril, julio y octubre."
            " Especificamente, el DataFrame proviene de la sección de Aberraciones Cromosómicas Recurrentes y en ella, solo se bajaron las aberraciones estructurales"
            "para ser analizadas, incluyendo todos los genes, topografías y morfologías y disponibles en la base de datos.", justified_body_style),
    ]
    # Crear un marco para organizar los párrafos
    frame = Frame(1 * inch, 1 * inch, width - 2 * inch, height - 2 * inch, showBoundary=0)
    # Dibujar el contenido dentro del marco
    frame.addFromList(contenido, c)
    # Guardar página
    c.showPage()


    # 3. Página de contenido y anexos
    c.setFont("Helvetica-Bold", 18)
    c.drawString(1 * inch, height - 1 * inch, "Contenido")
    c.setFont("Helvetica", 12)
    c.drawString(1 * inch, height - 1.5 * inch, "1.  Página de presentación")
    c.drawString(1 * inch, height - 2 * inch, "2.  Introducción")
    c.drawString(1 * inch, height - 2.5 * inch, "3.  Contenido y anexos")
    c.drawString(1 * inch, height - 3 * inch, "4.  Índice de gráficas")
    c.drawString(1 * inch, height - 3.5 * inch, "5.  Gráficas finales")

    c.setFont("Helvetica-Bold", 18)
    c.drawString(1 * inch, height - 4.5 * inch, "Índice de Anexos")
    c.setFont("Helvetica", 12)
    c.drawString(1 * inch, height - 5 * inch, "1. DataFrame Original")
    c.drawString(1 * inch, height - 5.5 * inch, "2. DataFrame Organizado")
    c.drawString(1 * inch, height - 6 * inch, "3. DataFrame Outliers")
    c.drawString(1 * inch, height - 6.5 * inch, "4. Diccionario de morfologías")
    c.showPage()


    # 4. Índice de gráficas
    c.setFont("Helvetica-Bold", 18)
    c.drawString(1 * inch, height - 1 * inch, "Índice de Gráficas")
    c.setFont("Helvetica", 12)
    c.drawString(1 * inch, height - 1.5 * inch, "1. Scatterplot: Gráfico de dispersión de outliers")
    c.drawString(1 * inch, height - 2 * inch, "2. Countplot: Gráfico de conteo de outliers por cromosoma")
    c.drawString(1 * inch, height - 2.5 * inch, "3. Barplot: Conteo de anormalidades por cromosoma y brazo")
    c.drawString(1 * inch, height - 3 * inch, "4. Barplot: Frecuencia de anormalidades")
    c.drawString(1 * inch, height - 3.5 * inch, "5. Heatmap: Tipo de anormalidad por cromosoma")
    c.drawString(1 * inch, height - 4 * inch, "6. Scatterplot: Distribución de casos por banda cromosómica")
    c.drawString(1 * inch, height - 4.5 * inch, "7. StripPLot: Casos por tipo de anormalidad")
    c.drawString(1 * inch, height - 5 * inch, "8. Heatmap: Frecuencia de anormalidades por banda y tipo")
    c.showPage()

    # 5. Gráficas finales (agregar cada gráfica en una página nueva)
    def agregar_graficas(c, path_graficas):
        for i in range(1, 9):  # Asumiendo que tienes 8 gráficas
            # Crear nombre de archivo de la gráfica
            grafica_path = os.path.join(path_graficas, f"grafica_{i}.png")
            # Verificar si la gráfica existe
            if os.path.exists(grafica_path):
                c.showPage()  # Crear una nueva página con el tamaño definido
                c.setFont("Helvetica-Bold", 14)
                c.drawString(1 * inch, height - 1 * inch, f"Gráfica {i}")

                # Ajustar el tamaño de la imagen
                img_width = 750  # Ancho de la imagen en puntos
                img_height = 400  # Alto de la imagen en puntos

                # Crear objeto de imagen
                img = Image(grafica_path, width = img_width, height = img_height)

                # Calcular la posición centrada en la página
                x_pos = (width - img_width) / 2  # Centrar horizontalmente
                y_pos = (height - img_height) / 2  # Centrar verticalmente

                # Dibujar la imagen en el PDF
                img.drawOn(c, x_pos, y_pos)

    # Llamar la función para agregar las gráficas
    agregar_graficas(c, path_graficas)

    # Finalizar el PDF
    c.save()


#===================================================================================================================#
#===================================================================================================================#
#====================================================COMANDOS=======================================================#
#===================================================================================================================#
#===================================================================================================================#

# Llamar la función para crear el reporte
crear_reporte_pdf("reporte_analisis_datos.pdf", df, df_sorted, path_graficas="./")

# Impresión de anexos como archivos .txt en carpeta
df.to_csv("DataFrame_original.txt", sep = ",", index = False)

df_sorted.to_csv("DataFrame_sorted.txt", sep=",", index = False)

outliers_iqr.to_csv("DataFrame_Outliers.txt", sep=",", index = False)

with open("dictmorfo_lxl.txt", "w") as file:
    for key, value in dict_morfo.items():
        file.write(f"{key}: {value}\n")



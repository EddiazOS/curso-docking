# Clase 1: Introducción al Docking Molecular

> **Semillero KAIROS** · IA, Ciencia de Datos y Bioinformática  
> Universidad de Cartagena · Facultad de Ciencias Exactas y Naturales  
> Código: KAIROS-2026-DOCK-01 | Fecha: 2026-03-15

---

## Tabla de Contenidos

1. [¿Qué es el Docking Molecular?](#1-qué-es-el-docking-molecular)
2. [Algoritmos de Docking](#2-algoritmos-de-docking)
3. [Proteínas: Estructura 3D y Actividad Biológica](#3-proteínas-estructura-3d-y-actividad-biológica)
4. [La Base de Datos PDB](#4-la-base-de-datos-pdb)
5. [Cómo Navegar el PDB: Información Relevante](#5-cómo-navegar-el-pdb-información-relevante)
6. [El Archivo `.pdb`: Anatomía y Lectura](#6-el-archivo-pdb-anatomía-y-lectura)
7. [Herramientas del Curso](#7-herramientas-del-curso)
8. [Interpretación de Resultados y Falsos Positivos](#8-interpretación-de-resultados-y-falsos-positivos)
9. [Tipos de Estudios de Docking](#9-tipos-de-estudios-de-docking)

---

## 1. ¿Qué es el Docking Molecular?

El **docking molecular** (o acoplamiento molecular) es una metodología computacional que predice la orientación preferida y la afinidad de unión de una molécula pequeña (ligando) dentro del sitio activo de una macromolécula receptora (generalmente una proteína), cuando ambas forman un complejo estable.

### Principios Básicos

| Principio | Descripción |
|-----------|-------------|
| **Reconocimiento molecular** | La interacción ligando–receptor es altamente específica, gobernada por complementariedad estérica y electrostática. |
| **Función de puntuación (*scoring function*)** | Estima la energía libre de unión (ΔG) del complejo. A menor ΔG (más negativo), mayor afinidad predicha. |
| **Flexibilidad conformacional** | El receptor y/o el ligando pueden tratarse como rígidos o flexibles durante la búsqueda. |
| **Espacio de búsqueda conformacional** | El algoritmo explora las posibles poses del ligando dentro de un volumen definido (caja de docking). |

### Aplicaciones principales

- **Diseño de fármacos asistido por computadora (CADD):** identificar y optimizar candidatos terapéuticos.
- **Virtual screening:** tamizar millones de compuestos en *in silico* antes de síntesis o ensayos biológicos.
- **Comprensión mecanística:** estudiar cómo un inhibidor bloquea un sitio activo enzimático.
- **Re-propósito de fármacos:** evaluar si medicamentos aprobados pueden actuar sobre nuevas dianas.

> **Referencia fundamental:** Morris, G. M. & Lim-Wilby, M. (2008). Molecular docking. *Methods in Molecular Biology*, 443, 365–382. https://doi.org/10.1007/978-1-59745-177-2_19

---

## 2. Algoritmos de Docking

Los algoritmos de docking tienen el reto de explorar eficientemente un espacio conformacional enorme. Se clasifican principalmente según la estrategia de búsqueda.

### 2.1 Clasificación general de algoritmos

| Categoría | Algoritmo | Característica |
|-----------|-----------|----------------|
| **Búsqueda sistemática** | Exhaustive search, fragmentación (DOCK) | Evalúa posiciones en una grilla, costoso computacionalmente |
| **Estocásticos / Evolutivos** | Algoritmos genéticos (AutoDock, AutoDock Vina) | Optimización inspirada en evolución natural; balance velocidad/precisión |
| **Monte Carlo** | GLIDE, ICM | Perturbaciones aleatorias aceptadas según criterio de Boltzmann |
| **Aprendizaje profundo** | DiffDock, EquiBind | Modelos de ML que predicen poses directamente sin función de energía clásica |

### 2.2 Software de pago vs. código abierto

**Herramientas de pago** como **Schrödinger Glide** o **MOE (Chemical Computing Group)** utilizan campos de fuerza propietarios (OPLS4, AMBER), funciones de puntuación altamente calibradas y flujos de trabajo gráficos avanzados. Están optimizadas para rendimiento industrial pero requieren licencias costosas.

**Herramientas open source** como **AutoDock Vina**, **GNINA** o **smina** utilizan algoritmos genéticos o gradiente de descenso con funciones de puntuación basadas en energía o ML, son completamente gratuitas y reproducibles, y son el estándar en entornos académicos.

> **Referencia:** Trott, O. & Olson, A. J. (2010). AutoDock Vina: Improving the speed and accuracy of docking. *Journal of Computational Chemistry*, 31(2), 455–461. https://doi.org/10.1002/jcc.21334

---

## 3. Proteínas: Estructura 3D y Actividad Biológica

Las proteínas son **polímeros de aminoácidos** que ejecutan prácticamente todas las funciones celulares: catálisis enzimática, señalización, transporte, defensa inmune y más. Su función está determinada de manera directa por su estructura tridimensional.

### Niveles de organización estructural

```
Secuencia de aa → Estructura primaria
Hélices α / Láminas β → Estructura secundaria
Plegamiento global → Estructura terciaria
Asociación de subunidades → Estructura cuaternaria
```

La estructura terciaria crea cavidades, surcos y bolsillos en la superficie proteica: estos son los **sitios de unión a ligandos** (sitios activos en enzimas, sitios alostéricos, etc.). La geometría y las propiedades fisicoquímicas de estos bolsillos determinan qué moléculas pueden unirse con alta afinidad, principio conocido como **complementariedad estérica y electrostática** (modelo "llave-cerradura" de Emil Fischer, 1894, y su extensión al modelo de "ajuste inducido" de Koshland, 1958).

> **Referencia estructural:** Branden, C. & Tooze, J. (1999). *Introduction to Protein Structure* (2nd ed.). Garland Publishing.

---

## 4. La Base de Datos PDB

El **Protein Data Bank (PDB)** es el repositorio global y de acceso libre para estructuras tridimensionales de macromoléculas biológicas determinadas experimentalmente. Fue fundado en 1971 y actualmente aloja más de **220,000 estructuras** (proteínas, ácidos nucleicos y complejos).

### Acceso

- **Portal principal:** [https://www.rcsb.org](https://www.rcsb.org)
- **Espejo europeo (PDBe):** [https://www.ebi.ac.uk/pdbe/](https://www.ebi.ac.uk/pdbe/)
- **Espejo japonés (PDBj):** [https://pdbj.org](https://pdbj.org)

Cada estructura se identifica con un código único de **4 caracteres alfanuméricos** (ej. `1HSG`, `6LU7`, `3HTB`).

> **Referencia:** Berman, H. M. et al. (2000). The Protein Data Bank. *Nucleic Acids Research*, 28(1), 235–242. https://doi.org/10.1093/nar/28.1.235

---

## 5. Cómo Navegar el PDB: Información Relevante

Al buscar una estructura en [rcsb.org](https://www.rcsb.org) e ingresar a su página, encontrarás información crítica organizada en varias pestañas. A continuación se describe qué buscar y cómo interpretarlo.

### 5.1 Métricas de calidad estructural

La calidad de una estructura cristalográfica de rayos X se evalúa con:

| Métrica | Valor aceptable | Qué indica |
|---------|-----------------|------------|
| **Resolución (Å)** | ≤ 2.5 Å (ideal < 2.0 Å) | Detalle atómico; menor valor = mejor resolución |
| **R-factor** | < 0.25 | Ajuste del modelo a los datos de difracción |
| **R-free** | < 0.30 (y próximo al R-factor) | Validación cruzada; diferencia grande indica sobreajuste |
| **Ramachandran plot** | > 95% residuos en región favorable | Geometría del esqueleto peptídico |
| **MolProbity score** | < 2.0 | Calidad global; escala tipo percentil |

Para estructuras de RM o crio-EM, los criterios difieren. En crio-EM se usa la **resolución FSC (0.143)** como indicador principal.

### 5.2 Información del sitio de unión y ligandos

En la pestaña **"Ligand"** o **"Small Molecules"** del PDB encontrarás:
- **Código HET** del ligando (ej. `LIG`, `ATP`, `NAG`): identificador de 3 caracteres.
- **Descripción química** y enlace a bases de datos de ligandos (ChEBI, PubChem, DrugBank).
- **Interacciones** con residuos del bolsillo de unión.

En la sección **"Binding Site"** o usando la herramienta integrada **"3D View" (Mol\*)** puedes visualizar el sitio activo, los residuos clave y medir distancias de interacción (puentes de H, interacciones hidrofóbicas, π-stacking).

### 5.3 Otros elementos a verificar

- **Organismos fuente:** asegúrate de que la proteína sea del organismo de tu interés (ej. *Homo sapiens* vs. proteína bacteriana homóloga).
- **Mutaciones o construcciones truncadas:** algunas estructuras usan proteínas mutadas para cristalización; esto puede afectar el sitio de unión.
- **Moléculas de agua y iones:** pueden ser importantes para el docking (aguas conservadas en el sitio activo).
- **Número de cadenas:** verificar si es monómero, dímero, etc. y cuál cadena usar.

---

## 6. El Archivo `.pdb`: Anatomía y Lectura

El formato de archivo **`.pdb`** es texto plano estructurado en columnas. Cada línea comienza con un **tipo de registro** de 6 caracteres que indica qué tipo de información contiene.

### Estructura general del archivo

```
HEADER    HYDROLASE                               ...
TITLE     CRYSTAL STRUCTURE OF HIV-1 PROTEASE...
REMARK 2  RESOLUTION. 2.00 ANGSTROMS.
SEQRES  1 A  99  MET GLN ILE THR LEU TRP ...
HETATM 300  C1  LIG A 201      10.123  20.456  30.789  ...
ATOM      1  N   MET A   1      ...
ATOM      2  CA  MET A   1      ...
CONECT  300  301  302
END
```

### Registros más importantes

| Registro | Descripción |
|----------|-------------|
| `HEADER` | Clasificación y fecha de depósito |
| `REMARK` | Notas diversas: resolución, método experimental, referencias |
| `SEQRES` | Secuencia de aminoácidos de cada cadena |
| `ATOM` | Coordenadas XYZ de átomos de la cadena polipeptídica |
| `HETATM` | Coordenadas de heteroátomos: ligandos, cofactores, agua, iones |
| `CONECT` | Conectividad entre átomos (importante para ligandos) |
| `END` | Marca el fin del archivo |

### Columnas del registro `ATOM` / `HETATM`

```
Columnas  1- 6   Tipo de registro (ATOM o HETATM)
Columnas  7-11   Número de átomo
Columnas 13-16   Nombre del átomo (CA = carbono alfa)
Columnas 18-20   Nombre del residuo (aminoácido o ligando)
Columna  22      ID de cadena (A, B, C...)
Columnas 23-26   Número de residuo
Columnas 31-38   Coordenada X (Å)
Columnas 39-46   Coordenada Y (Å)
Columnas 47-54   Coordenada Z (Å)
Columnas 55-60   Ocupancia (normalmente 1.00)
Columnas 61-66   Factor B (temperatura, flexibilidad atómica)
Columna  77-78   Símbolo del elemento
```

> **Para qué sirve:** Los archivos `.pdb` son el formato de entrada estándar para la mayoría de los programas de docking, visualización (PyMOL, UCSF Chimera) y preparación de receptores (AutoDockTools, Protein Preparation Wizard).

---

## 7. Herramientas del Curso

Durante el curso usaremos una combinación de herramientas open source, gratuitas y reproducibles. Aquí un primer vistazo:

| Herramienta | Función principal | Acceso |
|-------------|------------------|--------|
| **AutoDock Vina** | Motor de docking molecular | [vina.scripps.edu](https://vina.scripps.edu) |
| **AutoDockTools / MGLTools** | Preparación de receptor y ligando, generación de archivos `.pdbqt` | [mgltools.scripps.edu](http://mgltools.scripps.edu) |
| **PyMOL** | Visualización molecular 3D, inspección de estructuras | [pymol.org](https://pymol.org) (versión educativa gratuita) |
| **UCSF Chimera / ChimeraX** | Visualización y preparación avanzada | [cgl.ucsf.edu/chimera](https://www.cgl.ucsf.edu/chimera/) |
| **Open Babel** | Conversión entre formatos moleculares | [openbabel.org](http://openbabel.org) |
| **RCSB PDB** | Obtención de estructuras proteicas | [rcsb.org](https://www.rcsb.org) |
| **PubChem / ChEMBL** | Obtención de ligandos y datos de actividad biológica | [pubchem.ncbi.nlm.nih.gov](https://pubchem.ncbi.nlm.nih.gov) |

> Las instrucciones de instalación y configuración se detallarán en la sesión práctica correspondiente a cada herramienta.

---

## 8. Interpretación de Resultados y Falsos Positivos

Al terminar una corrida de docking, el software genera un conjunto de **poses** (orientaciones predichas del ligando) ordenadas por su **energía de unión estimada** (kcal/mol). El análisis correcto de estos resultados es tan importante como la corrida misma.

### 8.1 Parámetros clave de evaluación

| Parámetro | Descripción | Interpretación |
|-----------|-------------|----------------|
| **Energía de unión (ΔG predicho)** | Valor en kcal/mol; más negativo = mayor afinidad | Referencia general: < −7 kcal/mol indica unión potencialmente relevante |
| **RMSD entre poses** | Diferencia en Å entre la mejor pose y las demás | RMSD < 2 Å entre las mejores poses sugiere convergencia y confiabilidad |
| **Interacciones con residuos clave** | Puentes de hidrógeno, contactos hidrofóbicos, π-stacking con residuos catalíticos | Valida que el ligando ocupa el sitio de unión correcto |
| **Modo de unión visual** | Inspección en PyMOL o Chimera | Verificar que el ligando NO esté "saliendo" del bolsillo o en conformación imposible |

### 8.2 Identificación de Falsos Positivos

Un **falso positivo** es un compuesto que obtiene una buena puntuación de docking pero que en realidad **no tiene actividad biológica** real. Son el principal problema del virtual screening. Sus causas más frecuentes son:

- **Agregadores:** moléculas que forman agregados coloidales y "secuestran" proteínas de forma inespecífica. Revisar el perfil PAINS (Pan-Assay Interference Compounds) con [FAF-Drugs4](https://fafdrugs4.rpbs.univ-paris-diderot.fr/) o el filtro PAINS de RDKit.
- **Artefactos de la función de scoring:** las funciones de puntuación son aproximaciones; pueden sobreestimar la afinidad de ciertos grupos químicos (catecoles, quinonas, aldehídos reactivos).
- **Poses en bolsillos erróneos:** el ligando puede obtener buena puntuación en una cavidad que no es el sitio activo. Siempre verificar visualmente y comparar con el cristal co-cristalizado si existe.
- **Redocking como control:** si la proteína tiene un ligando co-cristalizado, re-dockearlo y verificar que el RMSD de la pose predicha respecto al cristal sea < 2 Å es el **estándar de validación** más usado.

### 8.3 Estrategias de validación

1. **Re-docking del ligando nativo:** RMSD < 2.0 Å es criterio de éxito del protocolo.
2. **Enriquecimiento con señuelos (DUD-E benchmark):** calcular métricas como AUC-ROC y factor de enriquecimiento.
3. **Análisis de interacciones farmacofóricas:** comparar interacciones del hit con el farmacoporo conocido del sitio activo.
4. **Filtros de druglikeness:** regla de Lipinski (Ro5), ADMET básico con SwissADME.

> **Referencia:** Shoichet, B. K. (2004). Virtual screening of chemical libraries. *Nature*, 432, 862–865. https://doi.org/10.1038/nature03197

---

## 9. Tipos de Estudios de Docking

El docking no es una metodología única; existen distintos diseños de estudio según el objetivo científico. A continuación se describen los más relevantes.

### Virtual Screening (VS)

Se busca identificar compuestos activos dentro de una **librería de miles o millones de moléculas** frente a una diana conocida. El objetivo es reducir el espacio químico a un conjunto de candidatos para ensayos biológicos posteriores. Es el uso más extendido del docking en drug discovery.

*Flujo típico:* librería filtrada por Ro5 → docking → ranking por score → inspección visual del top 1% → ensayos in vitro.

### Re-docking

Se utiliza como **control de calidad del protocolo**. Consiste en extraer el ligando co-cristalizado de su estructura PDB y volver a dockearlo en el mismo receptor. Si el programa reproduce la pose experimental (RMSD < 2 Å), el protocolo es confiable para ese sistema.

### Docking Cruzado (Cross-docking)

El ligando de una estructura PDB se dockea en el receptor de **otra estructura PDB del mismo blanco**, generalmente con una conformación diferente. Evalúa la robustez del método ante variaciones conformacionales del receptor.

### Docking Inverso (Reverse Docking / Target Fishing)

En lugar de buscar ligandos para un receptor fijo, se toma **un ligando de interés y se lo dockea contra una biblioteca de proteínas** para identificar posibles dianas. Es útil para estudiar polimorfismo farmacológico o re-propósito de fármacos.

### Docking con Receptor Flexible (Ensemble Docking)

Se generan **múltiples conformaciones del receptor** (desde DM, NMR o modelos) y se realiza docking en cada una. Captura el efecto del ajuste inducido y mejora la precisión en proteínas altamente flexibles.

### Docking Covalente

Modela la formación de un **enlace covalente** entre el ligando y un residuo específico del receptor (ej. Cys, Ser). Es crucial para inhibidores covalentes como los usados en KRAS G12C o EGFR.

| Tipo de Docking | Objetivo | Escala | Herramienta recomendada |
|-----------------|----------|--------|------------------------|
| Virtual Screening | Identificar hits en biblioteca | Millones de compuestos | Vina, GNINA, Glide (pago) |
| Re-docking | Validar protocolo | 1 ligando | AutoDock Vina |
| Cross-docking | Robustez del método | Pocas estructuras | AutoDock Vina |
| Docking Inverso | Identificar dianas (target fishing) | Múltiples proteínas | idTarget, PharmMapper |
| Ensemble Docking | Flexibilidad del receptor | Múltiples conformaciones | AutoDock Vina + clustering |
| Docking Covalente | Inhibidores covalentes | 1 ligando | CovDock (Schrödinger), AutoDock |

> **Referencia:** Pinzi, L. & Rastelli, G. (2019). Molecular docking: Shifting paradigms in drug discovery. *International Journal of Molecular Sciences*, 20(18), 4331. https://doi.org/10.3390/ijms20184331

---

## Material Complementario

- 📄 [AutoDock Vina — Documentación oficial](https://autodock-vina.readthedocs.io/en/latest/)
- 🗄️ [RCSB PDB — Tutorial de búsqueda avanzada](https://www.rcsb.org/docs/search-and-browse/advanced-search/introduction-to-advanced-search)
- 📘 [Introduction to Molecular Docking (EMBL-EBI Training)](https://www.ebi.ac.uk/training/online/courses/molecular-docking/)
- 🛠️ [SwissADME — Propiedades farmacoquímicas](http://www.swissadme.ch/)
- 🔬 [FAF-Drugs4 — Filtro PAINS](https://fafdrugs4.rpbs.univ-paris-diderot.fr/)

---

*Documento generado para el Semillero KAIROS · KAIROS-2026-DOCK-01 · 2026-03-15*  
*Dr. Ricardo Vivas Reyes (Director) · Cesar Tovio Gracia (Codirector)*

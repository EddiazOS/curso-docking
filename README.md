# Docking Molecular — Semillero de Química

> Curso de introducción a la química computacional mediante docking molecular.  
> 6 sesiones · 2 horas/sesión · Enfoque: terminal, Python y análisis crítico de resultados.

---

## Descripción general

Este repositorio contiene la información fundamental para hacer seguimiento a las clases de introducción a las técnicas de biología computacional a partir de docking molecular. Dentro de este se alojan todos los scripts y programas ejecutados para realizar tanto los ejercicios de demostración presentados en clase como los ejercicios propuestos para trabajo individual. Adicionalmente, toda la información base con la que se ha desarrollado el curso se encuentra citada en la sección de [Fuentes y recursos complementarios](#fuentes-y-recursos-complementarios).

El curso está desarrollado en su mayoría en base al lenguaje de programación Python, además de incluir frameworks y programas de línea de comando. Por lo tanto, se recomienda desarrollar los ejercicios dentro de un ordenador con sistema operativo Linux/macOS. Para el caso de los usuarios de Windows, lo más recomendable es instalar el subsistema de línea de comandos mediante **WSL2**.

Adicional a las herramientas de programación, se propone el uso de herramientas agénticas de IA conectadas a bases de datos y servicios científicos, con el fin de plantear flujos de trabajo más eficientes. El uso e interpretación de estos servicios no es obligatorio y queda reservado como recurso complementario.

---

## Objetivos

### Objetivo general

Capacitar a los estudiantes del semillero en el uso fundamentado de la metodología de docking molecular como herramienta de exploración estructural, desarrollando tanto la competencia técnica para ejecutar flujos de trabajo reproducibles como el criterio analítico para interpretar, validar y comunicar resultados en el contexto de una investigación química real.

### Objetivos específicos

- **Conceptuales:** Introducir los principios fisicoquímicos que sustentan el docking molecular: reconocimiento ligando-receptor, búsqueda conformacional y función de scoring como aproximación a la energía libre de unión. Identificar el alcance real de la técnica, sus supuestos simplificadores y las condiciones bajo las cuales sus resultados son confiables o deben complementarse con otros métodos.
- **Computacionales:** Configurar y gestionar un entorno computacional reproducible usando conda y Python, capaz de sostener un flujo de trabajo de docking completo, tanto desde la terminal como mediante scripts reutilizables.
- **Analíticos:** Establecer bases para evaluar críticamente los resultados de un docking más allá del valor numérico del score: seleccionar poses con criterio estructural, validar el protocolo mediante re-docking, analizar interacciones proteína-ligando e integrar propiedades fisicoquímicas del ligando en la interpretación.
- **Investigativos:** Diseñar y ejecutar un mini-screening virtual sobre una diana de interés propio, comunicando los resultados en formato de reporte científico corto con rigor metodológico.

---

## Estructura del curso

| Sesión | Título | Eje central | Reto de análisis |
|:---:|---|---|---|
| **S1** | Fundamentos: ¿qué predice realmente el docking? | Función de scoring, ΔG de unión, búsqueda conformacional, modelos de reconocimiento molecular, alcance y limitaciones | Buscar una proteína de interés en el PDB y evaluar si su estructura es adecuada para un estudio de docking |
| **S2** | Entorno computacional: terminal, Python y estructura de proyecto | Reproducibilidad computacional, gestión de entornos conda, organización de proyectos científicos | Usar RDKit para calcular propiedades básicas de un ligando y conectarlas con las limitaciones del docking discutidas en S1 |
| **S3** | Preparación del receptor, ligando y primer docking (GUI/terminal) | Inspección de PDB, aguas, ocupancias, cargas Gasteiger, torsiones rotables, definición del grid box | Dos criterios distintos de preparación del receptor producen scores diferentes — ¿cuál es correcto y cómo se documenta? |
| **S4** | Docking con Python: scripting, APIs y automatización | Python bindings de Vina, acceso programático a PDB/PubChem, flujo reproducible de punta a punta | El script falla silenciosamente en un ligando (poses fuera del grid) — proponer una validación programática |
| **S5** | Análisis de interacciones y validación del protocolo | Re-docking, RMSD, tipos de interacciones proteína-ligando, diagnóstico de fallos de protocolo | Dos ligandos con scores similares pero perfiles de interacción distintos — ¿cuál priorizar para ensayo experimental? |
| **S6** | Mini virtual screening y reporte científico | Ligand efficiency, Regla de Lipinski, priorización de candidatos, estructura de comunicación científica | El mejor score pertenece a un compuesto con bajo drug-likeness — integrar score, Lipinski y ligand efficiency para priorizar |

---

## Estructura del repositorio

```
curso-docking/
├── README.md
├── sesiones/
│   ├── S1_fundamentos/
│   ├── S2_entorno/
│   ├── S3_preparacion_docking_GUI/
│   ├── S4_docking_python/
│   ├── S5_analisis_validacion/
│   └── S6_virtual_screening/
├── datos/
│   ├── receptores/
│   └── ligandos/
└── entorno/
    └── environment.yml
```

Cada carpeta de sesión contiene:
- `README.md` con el material teórico, instrucciones del ejercicio de demostración y del ejercicio individual.
- Scripts de demostración (`.py` o `.sh`).
- Carpeta `ejercicio_propuesto/` con el enunciado y espacio para la solución del estudiante.

---

## Requerimientos técnicos

### Sistema operativo
- Linux o macOS (recomendado)
- Windows con **WSL2** habilitado

### Gestor de entornos
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) o [Mambaforge](https://github.com/conda-forge/miniforge)

### Entorno conda

Para crear el entorno del curso con todas las dependencias:

```bash
conda env create -f entorno/environment.yml
conda activate docking-curso
```

### Herramientas principales

| Herramienta | Rol | Instalación |
|---|---|---|
| [AutoDock Vina 1.2.0](https://vina.scripps.edu/) | Motor de docking | `conda install -c conda-forge autodock-vina` |
| [Meeko](https://github.com/forlilab/Meeko) | Preparación de receptor y ligando | `pip install meeko` |
| [Open Babel](https://openbabel.org/) | Conversión de formatos moleculares | `conda install -c conda-forge openbabel` |
| [RDKit](https://www.rdkit.org/) | Quimioinformática y cálculo de propiedades | `conda install -c conda-forge rdkit` |
| [PLIP](https://github.com/pharmai/plip) | Análisis de interacciones proteína-ligando | `pip install plip` |
| [PubChemPy](https://pubchempy.readthedocs.io/) | Acceso programático a PubChem | `pip install pubchempy` |
| [PyPDB](https://github.com/williamgilpin/pypdb) | Acceso programático al PDB | `pip install pypdb` |
| [PyMOL](https://pymol.org/) | Visualización molecular (GUI) | Descarga educativa gratuita |

### Bases de datos

- [RCSB Protein Data Bank (PDB)](https://www.rcsb.org/) — estructuras cristalográficas de proteínas
- [PubChem](https://pubchem.ncbi.nlm.nih.gov/) — estructuras y propiedades de ligandos
- [ChEMBL](https://www.ebi.ac.uk/chembl/) — datos de actividad biológica de compuestos
- [UniProt](https://www.uniprot.org/) — información funcional de proteínas

---

## Fuentes y recursos complementarios

### Referencias metodológicas fundamentales

- Eberhardt, J. et al. (2021). **AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings**. *Journal of Chemical Information and Modeling*, 61(8), 3891–3898. https://doi.org/10.1021/acs.jcim.1c00203
- Trott, O. & Olson, A. J. (2010). **AutoDock Vina: Improving the speed and accuracy of docking**. *Journal of Computational Chemistry*, 31(2), 455–461. https://doi.org/10.1002/jcc.21334
- Pagadala, N. S. et al. (2017). **Software for molecular docking: a review**. *Biophysical Reviews*, 9, 91–102. https://doi.org/10.1007/s12551-016-0247-1
- Torres, P. H. M. et al. (2019). **Ten quick tips to perform meaningful and reproducible molecular docking**. *PLOS Computational Biology*, 15(4), e1006829. https://doi.org/10.1371/journal.pcbi.1006829
- Meng, X.-Y. et al. (2011). **Molecular Docking: A Powerful Approach for Structure-Based Drug Discovery**. *Current Computer-Aided Drug Design*, 7(2), 146–157.

### Documentación técnica

- [AutoDock Vina — Documentación oficial](https://autodock-vina.readthedocs.io/)
- [RDKit — Documentación oficial](https://www.rdkit.org/docs/)
- [Meeko — GitHub](https://github.com/forlilab/Meeko)
- [PLIP — GitHub](https://github.com/pharmai/plip)

### Recursos educativos

- Dallakyan, S. & Olson, A. J. (2015). **Demonstration of AutoDock as an Educational Tool for Drug Discovery**. *Methods in Molecular Biology*, 1263, 243–250.
- Sapay, N. & Bhattacharjee, N. (2020). **Molecular Docking Using Chimera and Autodock Vina for Non-bioinformaticians**. *JMIR Bioinformatics and Biotechnology*. https://doi.org/10.2196/14232

---

## Licencia

Este material es de uso educativo para el semillero de química. Los scripts y materiales originales de este repositorio se distribuyen bajo licencia [MIT](LICENSE).

# MSc Thesis: Semantic Standardisation and Integration of Non-Coding RNA Data

## Abstract 

Most of the human genome consists of non-coding RNA molecules (ncRNA), making their study
and understanding fundamental in modern molecular biology. These ncRNA play a crucial role in gene
regulation and expression, becoming essential elements for diagnosing, developing, and treating various
anomalies in fields such as physiology and genetics.
In recent years, sequencing technologies have facilitated the identification of a wide diversity of these
molecules. However, the integration of these data has not kept pace with their discovery, hindering their
effective management and analysis.
Although there are numerous resources and databases, such as BioGateway, that provide relevant
information on genetic elements, many lack comprehensive information on ncRNA, which are crucial for
gene regulation. To address this gap, the aim is to develop a ncRNA knowledge graph that is interoperable and integrable with BioGateway and other resources, ensuring data standardization and facilitating
interoperability across different sources.
In this work, a schema has been designed to model the available information on these molecules, inte grating five datasets. This schema will describe the data of the subsequently created graph. Additionally,
to promote interoperability, integration with the BioGateway knowledge networks is sought, using their
specific ontologies and annotation resources, as well as with Biolink, adopting its data model to structure
and relate information with other biological data. Finally, different RDFs are generated and loaded into
Virtuoso, enabling queries of greater biological relevance than those permitted by the source databases.

**Keywords**: Non-Coding RNA (ncRNA), Gene Regulation, Data Integration, Data Semantics, Know-
ledge Graphs, Resource Description Framework (RDF), Ontologies, BioGateway.
## Workflow

![Copia de PASO 1 (1)](https://github.com/user-attachments/assets/67916b99-56f0-4d63-8247-cb062f999719)


## Tools

- Employing semantic web tools to ensure standardized information representation.
- Using the RDFLib library in Python to generate RDF files.
- Utilizing RDF and SPARQL for effective data modeling and querying.
- Leveraging platforms such as Virtuoso for managing knowledge graphs.


## Directories Description

- [Code](./Code): This folder contains source code files used in the project, developed to standarize the data ([**R**](./Code/R)) and to generate the RDF files ([**RDFLib**](./Code/RDFLib)).

- [DB](./DB): This folder contains the datasets of the different databases used in the work: **Ensembl**, **RefSeq**, **LNCipedia**, **miRBase** and **miRCancer**.

- [RDF](./RDF): This folder contains the RDF files generated, classified according to the database to which they belong.

- [Turtle](./Turtle): This folder contains the Turtle RDF files that define the properties and classes used.

## Contact

For questions, collaborations, or suggestions, please contact the project's lead author via [noeliamoreng@gmail.com].



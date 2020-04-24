This Github repository contains the code and data to used for in the published study: "Mutational signature in colorectal cancer caused by genotoxic *pks+ E.coli*", C. Pleguezuelos-Manzano, J. Puschhof and A. Rosendahl Huber et al. 

https://www.nature.com/articles/s41586-020-2080-8
DOI: https://doi.org/10.1038/s41586-020-2080-8

For questions or suggestions: Please contact a.k.m.rosendahlhuber@prinsesmaximacentrum.nl

Please cite this paper when using this code for your research. BAM files containing the raw sequencing dat have been deposited at the European Genome-phenome Archive (https://ega-archive.org/) under the accesion code: EGAS00001003934.

SBS-*pks* and ID-*pks* signatures can be found in the folder Output. 

### ---- Usage ----- 
Workflow: For the most easily deployment, unzip Genotoxic_Ecoli.zip folder. Set working directory in script #1 and execute R scripts following numbering 1-9.
Dependencies: R version 3.6.0

1. Load data and analyze single base substitution load

2. Analysis of indels 

3. Transcriptional strand bias

4. Wider single base substitution and indel context analysis

5. Analyze data of second exposed organoid line 

6. Analysis of organoids exposed to recomplemented *E.coli* strain

7. Presence of *pks*-patterns in > 1bp deletions 

NOTE: scripts 8-9 use patient-level somatic variant and clinical data have been obtained from the Hartwig Medical Foundation under the data request number DR-084. Somatic variant and clinical data are freely available for academic use from the Hartwig Medical Foundation through standardized procedures. Privacy and publication policies, including co-authorship policies, can be retrieved from: https://www.hartwigmedicalfoundation.nl/en/data-policy/. 
Data request forms can be downloaded from https://www.hartwigmedicalfoundation.nl/en/applying-for-data/.
To gain access to the data, this data request form should be emailed to info@hartwigmedicalfoundation.nl., upon which it will be evaluated within 6 weeks by the HMF Scientific Council and an independent Data Access Board.
When access is granted, the requested data become available through a download link provided by HMF.

8. SBS and indel refitting of mutational signatures in HMF data

9. Presence of *pks*-motifs in > 1bp deletions and flanking bases in HMF data. 
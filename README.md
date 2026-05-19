# rna-central-pipeline
Import RNACentral IDs as 'external database ids' for all non-coding genes in RGD to generate links to RNACentral database.

LOGIC:

 1. download file with RNACentral RefSeq mappings:
   https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/refseq.tsv

 2. match RefSeq ids against RGD database to find gene rgd id
   - first match against transcripts table and look for RefSeq accession
   - if not found, try to look for RefSeq nucleotide external database id
     (for genes that have no transcripts loaded, f.e. because they are not in current annotation release)
   - if multiple genes are matched, report that

 3. generate RNACentral external database ids for the matching gene rgd ids

 4. for rat only, download file with RNACentral RGD mappings:
   https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/rgd.tsv
   and match by RGD ID directly

 5. download file with RNACentral Ensembl mappings:
   https://ftp.ebi.ac.uk/pub/databases/RNAcentral/current_release/id_mapping/database_mappings/ensembl.tsv
   and match by Ensembl gene id (stripping any version suffix)

 6. generate RNACentral external database ids for all matching gene rgd ids

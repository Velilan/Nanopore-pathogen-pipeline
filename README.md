# Citrus Pathogen Detection Pipeline (Nanopore)

Pipeline for detection and identification of major citrus pathogens from Oxford Nanopore sequencing data.

## Target pathogens

- Candidatus Liberibacter asiaticus (CaLas)
- Candidatus Liberibacter africanus (CaLaf)
- Candidatus Liberibacter americanus (CaLam)
- Candidatus Liberibacter solanacearum (CaLsol)
- Liberibacter crescens
- Xylella fastidiosa
- Xanthomonas citri
- Phyllosticta citricarpa

## Workflow

Nanopore FASTQ files
↓
FastpLong quality filtering
↓
Kraken2 taxonomic classification
↓
Extraction of pathogen-associated reads
↓
Minimap2 mapping against custom reference database
↓
Coverage, depth and identity calculations
↓
Automated pathogen detection report

## Dependencies

- Python 3
- fastplong
- kraken2
- minimap2
- samtools

## Usege

python pipeline_patogenosV16.py \
--run-dir datos \
--work-dir resultados \
--kraken-db base_datos \
--threads 16

## Read filtering

Reads are filtered using FastpLong with:

- Minimum quality score: Q10
- Maximum unqualified bases: 60%

Outputs:

- filtered.fastq
- fastplong.html
- fastplong.json

## Detection strategy

### Liberibacter

1. Kraken2 identifies candidate Liberibacter reads.
2. Reads are mapped against a custom set of Liberibacter references.
3. Species assignment is based primarily on:
- Number of reads with ≥95% identity
- Mean alignment identity
- Reads with ≥99% identity

### Xylella fastidiosa

Positive detection:

- ≥10 mapped reads
- ≥0.1% genome coverage

### Xanthomonas citri

Positive detection:

- ≥10 mapped reads
- ≥0.1% genome coverage

### Phyllosticta citricarpa

Positive detection:

- ≥10 mapped reads
- ≥0.1% genome coverage

Identity statistics include:

- Reads ≥80%
- Reads ≥85%
- Reads ≥90%
- Reads ≥95%
- Reads ≥97%
- Reads ≥99%
- Mean identity
- Median identity

## Output files

For each sample/barcode:

- merged.fastq
- filtered.fastq
- fastplong.html
- fastplong.json
- kraken.out
- kraken.report
- QC_summary.txt
- Pathogen_summary.txt
- Pathogen_summary.html

Additional BAM files are generated for each pathogen reference.

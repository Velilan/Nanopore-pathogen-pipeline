Citrus Pathogen Detection Pipeline (Nanopore)

Pipeline for detection and identification of major citrus pathogens from Oxford Nanopore sequencing data using quality filtering, taxonomic classification, reference mapping and alignment identity-based validation.

Target pathogens
Candidatus Liberibacter asiaticus (CaLas)
Candidatus Liberibacter africanus (CaLaf)
Candidatus Liberibacter americanus (CaLam)
Candidatus Liberibacter solanacearum (CaLsol)
Liberibacter crescens (BT-0, BT-1)
Xylella fastidiosa
Xanthomonas citri
Phyllosticta citricarpa

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

## Usage
conda env create -f environment.yml
conda activate EMERGENOW

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

Candidate reads are identified by Kraken2 and mapped against a custom database containing:

-CaLas
-CaLaf
-CaLam
-CaLsol
-Liberibacter crescens BT-0
-Liberibacter crescens BT-1

Species assignment is ranked using:

-Reads >99% identity
-Reads >97% identity
-Reads >95% identity
-Mean alignment identity

Detection states:

-POSITIVE
-WEAK DETECTION (INCONCLUSIVE)
-NEGATIVE

### Xylella fastidiosa

Subspecies assignment is ranked using:

-Reads >99% identity
-Reads >97% identity
-Reads >95% identity
-Mean alignment identity

### Xanthomonas citri

Positive detection:

Subspecies assignment is ranked using:

-Reads >99% identity
-Reads >97% identity
-Reads >95% identity
-Mean alignment identity

### Phyllosticta citricarpa

Reference assignment is ranked using:

-Reads >99% identity
-Reads >97% identity
-Reads >95% identity
-Mean alignment identity

For P. citricarpa, additional identity ratios are reported:
95% ratio = reads95 / mapped reads
97% ratio = reads97 / mapped reads
99% ratio = reads99 / mapped reads
These metrics help distinguish true matches from large numbers of low-quality alignments.

Alignment statistics

The following metrics are calculated for every reference:

-Mapped reads
-Genome coverage (%)
-Mean depth (X)
-Mean identity (%)
-Reads >80% identity
-Reads >85% identity
-Reads >90% identity
-Reads >95% identity
-Reads >97% identity
-Reads >99% identity

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
Notes

This pipeline was developed for rapid detection of citrus bacterial and fungal pathogens using Oxford Nanopore sequencing and a custom Kraken2 plus Minimap2 validation workflow. Species calls are based primarily on high-identity read support rather than raw mapping counts alone.

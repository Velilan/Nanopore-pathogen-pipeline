#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import logging
import subprocess
import shutil
from pathlib import Path
from collections import defaultdict
import statistics
import re

# ---------------- LOG ----------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

log = logging.getLogger(__name__)

# ---------------- DEPENDENCIAS ----------------

def which(cmd):
    return shutil.which(cmd)

def check_dependencies():

    required = [
        "fastplong",
        "kraken2",
        "minimap2",
        "samtools"
    ]

    missing = [x for x in required if which(x) is None]

    if missing:
        raise RuntimeError(
            f"Program missing: {', '.join(missing)}"
        )

# ---------------- CMD ----------------

def run_cmd(cmd):

    log.info(" ".join(cmd))
    subprocess.run(cmd, check=True)

# ---------------- BARCODES ----------------

def find_barcodes(run_dir):

    return [
        p for p in Path(run_dir).glob("barcode*")
        if p.is_dir()
    ]

# ---------------- MERGE FASTQ ----------------

def merge_fastqs(files, output):
    log.info(f"Mergeando {len(files)} FASTQ")

    with open(output, "wb") as out:

        for f in files:

            if str(f).endswith(".gz"):

                subprocess.run(
                    ["gunzip", "-c", str(f)],
                    stdout=out,
                    check=True
                )

            else:

                with open(f, "rb") as infile:
                    out.write(infile.read())

def count_fastq_reads(fastq):
    n = 0

    with open(fastq) as f:
        for _ in f:
            n += 1

    return n // 4


def calculate_mean_qscore(fastq):

    total_q = 0
    total_bases = 0

    with open(fastq) as f:
        while True:

            h = f.readline()
            if not h:
                break

            seq = f.readline()
            f.readline()
            qual = f.readline().strip()

            for c in qual:
                total_q += ord(c) - 33
                total_bases += 1

    if total_bases == 0:
        return 0

    return round(total_q / total_bases, 2)



# ---------------- TAXIDS ----------------

TARGETS = {
    "CaLas": {34021},
    "CaLaf": {34020},
    "CaLam": {309868, 1261131},
    "CaLsol": {556287, 658172},
    "Liberibacter_crescens": {1215346},
    "Liberibacter_genus": {34019},
    "Xanthomonas_citri": {346},
    "Xylella_fastidiosa": {2371},
    "P_citricarpa": {55181}
}

LIBERIBACTER_REFS = {
    "CaLas":
        "validation_refs/liberibacter/CaLas.fasta",

    "CaLaf":
        "validation_refs/liberibacter/CaLaf.fasta",

    "CaLam":
        "validation_refs/liberibacter/CaLam.fasta",

    "CaLsol":
        "validation_refs/liberibacter/CaLsol.fasta",

    "BT0":
        "validation_refs/liberibacter/Liberibacter crescens BT-0.fasta",

    "BT1":
        "validation_refs/liberibacter/Liberibacter crescens BT-1.fasta"
}

XYLELLA_REFS = {
    "multiplex":
        "validation_refs/xylella/multiplex.fasta",

    "fastidiosa":
        "validation_refs/xylella/Xylella fastidiosa subsp. fastidiosa.fasta",

    "pauca":
        "validation_refs/xylella/Xylella fastidiosa subsp. pauca.fasta",

    "sandyi":
        "validation_refs/xylella/Xylella fastidiosa subsp. sandyi.fasta",

    "morus":
        "validation_refs/xylella/Xylella fastidiosa subsp. morus.fasta"
}

XANTHOMONAS_REFS = {
    "citri":
        "validation_refs/xanthomonas/Xanthomonas citri pv. citri strain MN12.fasta",

    "aurantifolii":
        "validation_refs/xanthomonas/Xanthomonas citri pv. aurantifolii strain 1622.fasta",

    "anacardii":
        "validation_refs/xanthomonas/Xanthomonas citri pv. anacardii strain T51.fasta",

    "bilvae":
        "validation_refs/xanthomonas/Xanthomonas citri pv. bilvae strain NCPPB 3213.fasta",

    "durantae":
        "validation_refs/xanthomonas/Xanthomonas citri pv. durantae strain LMG696.fasta",

    "fuscans":
        "validation_refs/xanthomonas/Xanthomonas citri pv. fuscans strain ISO118C1.fasta",

    "glycines":
        "validation_refs/xanthomonas/Xanthomonas citri pv. glycines strain ICMP5732.fasta",

    "malvacearum":
        "validation_refs/xanthomonas/Xanthomonas citri pv. malvacearum strain CFBP 2036.fasta",

    "mangiferaeindicae":
        "validation_refs/xanthomonas/Xanthomonas citri pv. mangiferaeindicae strain T11.fasta",

    "punicae":
        "validation_refs/xanthomonas/Xanthomonas citri pv. punicae strain 119.fasta",

    "vignicola":
        "validation_refs/xanthomonas/Xanthomonas citri pv. vignicola strain CFBP7111.fasta"
}

PCIT_REFS = {
    "CBS": "validation_refs/phyllosticta_multi/CBS.fna",
    "CPC": "validation_refs/phyllosticta_multi/CPC.fna",
    "Gc12": "validation_refs/phyllosticta_multi/Gc12.fna",
    "LGMF06": "validation_refs/phyllosticta_multi/LGMF06.fna",
    "LSVM1123": "validation_refs/phyllosticta_multi/LSVM1123.fna",
    "LSVM359": "validation_refs/phyllosticta_multi/LSVM359.fna",
    "USA106": "validation_refs/phyllosticta_multi/USA106.fna",
    "USA109": "validation_refs/phyllosticta_multi/USA109.fna",
    "USA110": "validation_refs/phyllosticta_multi/USA110.fna",
    "USA111": "validation_refs/phyllosticta_multi/USA111.fna",
    "USA119": "validation_refs/phyllosticta_multi/USA119.fna",
    "USA120": "validation_refs/phyllosticta_multi/USA120.fna",
    "USA121": "validation_refs/phyllosticta_multi/USA121.fna",
    "USA122": "validation_refs/phyllosticta_multi/USA122.fna",
    "USA14": "validation_refs/phyllosticta_multi/USA14.fna",
    "USA27": "validation_refs/phyllosticta_multi/USA27.fna",
    "USA28": "validation_refs/phyllosticta_multi/USA28.fna",
    "USA43": "validation_refs/phyllosticta_multi/USA43.fna",
    "USA48": "validation_refs/phyllosticta_multi/USA48.fna",
    "USA66": "validation_refs/phyllosticta_multi/USA66.fna",
    "USA68": "validation_refs/phyllosticta_multi/USA68.fna",
    "USA73": "validation_refs/phyllosticta_multi/USA73.fna",
    "USA74": "validation_refs/phyllosticta_multi/USA74.fna"
    }

LIBERIBACTER_NAMES = {
    "CaLas": "Candidatus Liberibacter asiaticus",
    "CaLaf": "Candidatus Liberibacter africanus",
    "CaLam": "Candidatus Liberibacter americanus",
    "CaLsol": "Candidatus Liberibacter solanacearum",
    "BT0": "Liberibacter crescens BT-0",
    "BT1": "Liberibacter crescens BT-1"
}


def calculate_identity_stats(bam_file):
    cmd = ["samtools", "view", "-F", "4","-q", "20", str(bam_file)]

    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        check=True
    )

    identities = []

    reads95 = set()
    reads90 = set()
    reads85 = set()
    reads80 = set()

    for line in result.stdout.splitlines():

        fields = line.split("\t")

        read_id = fields[0]

        cigar = fields[5]

        nm = None

        for tag in fields[11:]:

            if tag.startswith("NM:i:"):
                nm = int(tag.split(":")[-1])
                break

        if nm is None:
            continue

        aligned_len = 0

        for length, op in re.findall(r'(\d+)([MIDNSHP=X])', cigar):

            if op in ["M", "=", "X", "I", "D"]:
                aligned_len += int(length)

        if aligned_len == 0:
            continue

        identity = (
            (aligned_len - nm)
            / aligned_len
        ) * 100

        identities.append(identity)

        if identity >= 80:
            reads80.add(read_id)

        if identity >= 85:
            reads85.add(read_id)

        if identity >= 90:
            reads90.add(read_id)

        if identity >= 95:
            reads95.add(read_id)

    if not identities:

        return {
            "mean_identity": 0,
            "median_identity": 0,
            "reads99": 0,
            "reads97": 0,
            "reads95": 0,
            "reads90": 0,
            "reads85": 0,
            "reads80": 0
        }

    print("\n===== IDENTITY DISTRIBUTION =====")
    print("Total mapped:", len(identities))
    print(">=80% :", len([x for x in identities if x >= 80]))
    print(">=85% :", len([x for x in identities if x >= 85]))
    print(">=90% :", len([x for x in identities if x >= 90]))
    print(">=95% :", len([x for x in identities if x >= 95]))
    print(">=97% :", len([x for x in identities if x >= 97]))
    print(">=99% :", len([x for x in identities if x >= 99]))
    print("Mean identity:", round(statistics.mean(identities), 2))
    print("Median identity:", round(statistics.median(identities), 2))
    print("===============================\n")


    return {
        "mean_identity": round(
            statistics.mean(identities), 2
        ),
        "median_identity": round(
            statistics.median(identities), 2
        ),
        "reads99": len([x for x in identities if x >= 99]),
        "reads97": len([x for x in identities if x >= 97]),
        "reads95": len(reads95),
        "reads90": len(reads90),
        "reads85": len(reads85),
        "reads80": len(reads80)
    }

def map_against_references(
    fastq,
    refs,
    output_dir,
    prefix,
    db_dir
):

    resultados = {}

    if not fastq.exists():
        return resultados

    for nombre, ref in refs.items():

        ref = Path(db_dir) / ref

        bam = output_dir / f"{prefix}_{nombre}.bam"

        cmd = (
            f'minimap2 -ax map-ont "{ref}" "{fastq}" '
            f'| samtools sort -o "{bam}"'
        )

        subprocess.run(
            cmd,
            shell=True,
            check=True
        )

        cmd_reads = f"samtools view -F 4 -q 20 '{bam}' | cut -f1 | sort -u | wc -l"

        mapped_reads = int(
            subprocess.check_output(
                cmd_reads,
                shell=True,
                text=True
            ).strip()
        )

        cov_result = subprocess.run(
            [
                "samtools",
                "coverage",
                str(bam)
            ],
            capture_output=True,
            text=True,
            check=True
        )
        coverage_pct = 0.0
        mean_depth = 0.0

        lines = cov_result.stdout.strip().splitlines()

        if len(lines) > 1:
            genome_size_total = 0
            covbases_total = 0
            depth_weighted = 0.0

            for line in lines[1:]:

                fields = line.split()

                length = int(fields[2]) - int(fields[1]) + 1
                covbases = int(fields[4])
                depth = float(fields[6])

                genome_size_total += length
                covbases_total += covbases
                depth_weighted += length * depth

            if genome_size_total > 0:

                coverage_pct = (
                    covbases_total / genome_size_total
                    ) * 100

                mean_depth = (
                    depth_weighted / genome_size_total
                )
        identity_stats = calculate_identity_stats(bam)

        resultados[nombre] = {
            "reads": mapped_reads,
            "coverage": coverage_pct,
            "depth": mean_depth,
            "mean_identity":
                identity_stats["mean_identity"],
            "median_identity":
                identity_stats["median_identity"],
            "reads99":
                identity_stats["reads99"],
            "reads97":
                identity_stats["reads97"],
            "reads95":
                identity_stats["reads95"],
            "reads90":
                identity_stats["reads90"],
            "reads85":
                identity_stats["reads85"],
            "reads80":
                identity_stats["reads80"]
        }

    return resultados

# ---------------- INTERPRETACION ----------------

def interpret_pathogens(
    targets_reads,
    output_dir,
    liberibacter_results,
    xylella_results,
    xanthomonas_results,
    pcit_results
):

    resumen = []

    resumen.append("===== PATHOGEN DETECTION SUMMARY =====")
    resumen.append("")

# =====================================================
# LIBERIBACTER
# =====================================================

    liberibacter_total = sum(v["reads"] for v in liberibacter_results.values())

    resumen.append("")
    resumen.append("=== LIBERIBACTER SPECIES DETECTION ===")

    if liberibacter_total >= 1:

        ordenados = sorted(
            liberibacter_results.items(),
            key=lambda x: (
                x[1]["reads95"],
                x[1]["mean_identity"],
                x[1]["reads99"]
            ),
            reverse=True
        )

        dominante = ordenados[0][0]
        dominante_reads = ordenados[0][1]["reads"]
        dominante_coverage = ordenados[0][1]["coverage"]
        dominante_depth = ordenados[0][1]["depth"]
        dominante_identity = ordenados[0][1]["mean_identity"]
        dominante_reads80 = ordenados[0][1]["reads80"]
        dominante_reads85 = ordenados[0][1]["reads85"]
        dominante_reads90 = ordenados[0][1]["reads90"]
        dominante_reads95 = ordenados[0][1]["reads95"]

        segunda_coverage = 0

        if len(ordenados) > 1:
            segunda_coverage = ordenados[1][1]["coverage"]

        diferencia = dominante_coverage - segunda_coverage

        if (
            dominante_reads95 >= 1
            and dominante_identity >= 95
        ):

            resumen.append(
                "Candidatus Liberibacter: POSITIVE"
            )

            resumen.append(
                f"Species assignment: "
                f"{LIBERIBACTER_NAMES[dominante]}"
            )
        else:

            resumen.append(
                "Candidatus Liberibacter: INCONCLUSIVE"
            )

            resumen.append(
                f"Best matching species: "
                f"{LIBERIBACTER_NAMES[dominante]}"
            )

        resumen.append(
            f"Mapped reads: {dominante_reads}"
        )

        resumen.append(
            f"Genome coverage: {dominante_coverage:.2f}%"
        )
        resumen.append(
            f"Mean depth: {dominante_depth:.2f}X"
        )

        resumen.append(
            f"Mean identity: "
            f"{ordenados[0][1]['mean_identity']:.2f}%"
        )

        resumen.append(
            f"Liberibacter total mapped reads: {liberibacter_total}"
        )

        resumen.append(
            f"Dominant species reads95: {dominante_reads95}"
        )

        resumen.append(
            f"Reads >90% identity: "
            f"{ordenados[0][1]['reads90']}"
        )

        resumen.append(
            f"Reads >95% identity: "
            f"{ordenados[0][1]['reads95']}"
        )

        resumen.append(
            f"Reads >97% identity: "
            f"{ordenados[0][1]['reads97']}"
        )

        resumen.append(
            f"Reads >99% identity: "
            f"{ordenados[0][1]['reads99']}"
        )

        resumen.append(
            f"Difference from second-ranked species: "
            f"{diferencia:.2f}%"
        )

        resumen.append("")
        resumen.append("Top matching species:")

        top_hits = ordenados[:4]

        for especie, datos in top_hits:

            ratio95 = (
                datos["reads95"] / datos["reads"] * 100
                if datos["reads"] > 0 else 0
            )

            ratio99 = (
                datos["reads99"] / datos["reads"] * 100
                if datos["reads"] > 0 else 0
            )

            resumen.append(
                f" {especie}: "
                f"{datos['reads']} reads, "
                f"{datos['coverage']:.2f}% genome coverage, "
                f"{datos['depth']:.2f}X, "
                f"Mean identity={datos['mean_identity']:.2f}%, "
                f"Reads >95% identity: {datos['reads95']}, "
                f"Reads >99% identity: {datos['reads99']}  "
                f"95% ratio: {ratio95:.2f}%, "
                f"99% ratio: {ratio99:.2f}%"
            )

    else:

        resumen.append(
            "Candidatus Liberibacter: NEGATIVE"
        )

# =====================================================
# XYLELLA
# =====================================================

    xylella_total = sum(
        v["reads"]
        for v in xylella_results.values()
    )

    resumen.append("")
    resumen.append("=== XYLELLA FASTIDIOSA===")

    if xylella_total >= 5:

        ordenados = sorted(
            xylella_results.items(),
            key=lambda x: x[1]["coverage"],
            reverse=True
        )

        dominante = ordenados[0][0]
        dominante_reads = ordenados[0][1]["reads"]
        dominante_coverage = ordenados[0][1]["coverage"]
        dominante_depth = ordenados[0][1]["depth"]
        dominante_identity = ordenados[0][1]["mean_identity"]
        dominante_reads90 = ordenados[0][1]["reads90"]
        dominante_reads95 = ordenados[0][1]["reads95"]

        segunda_coverage = 0

        if len(ordenados) > 1:
            segunda_coverage = ordenados[1][1]["coverage"]

        diferencia = dominante_coverage - segunda_coverage

        if dominante_reads >= 10 and dominante_coverage >= 0.1:

            resumen.append(
                "Xylella fastidiosa: POSITIVE"
            )
            resumen.append(
                f"Best matching subspecies: {dominante}"
            )

            resumen.append(
                f"Mapped reads: {dominante_reads}"
            )

            resumen.append(
                f"Genome coverage: {dominante_coverage:.2f}%"
            )

            resumen.append(
                f"Mean depth: {dominante_depth:.2f}X"
            )

            resumen.append(
                f"Mean identity: "
                f"{ordenados[0][1]['mean_identity']:.2f}%"
            )

            resumen.append(
                f"Xylella total mapped reads: {xylella_total}"
            )

            resumen.append(
                f"Dominant species reads95: {dominante_reads95}"
            )

            resumen.append(
                f"Reads >90% identity: "
                f"{ordenados[0][1]['reads90']}"
            )

            resumen.append(
                f"Reads >95% identity: "
                f"{ordenados[0][1]['reads95']}"
            )

            resumen.append("")
            resumen.append("Top matching references:")

            top_hits = ordenados[:3]

            for especie, datos in top_hits:

                resumen.append(
                    f" {especie}: "
                    f"{datos['reads']} reads, "
                    f"{datos['coverage']:.2f}% genome coverage, "
                    f"{datos['depth']:.2f}X, "
                    f"Mean identity={datos['mean_identity']:.2f}%, "
                    f"Reads >95% identity: {datos['reads95']}"
                )

        elif dominante_coverage > 0:

            resumen.append(
                "Xylella fastidiosa: WEAK DETECTION (INCONCLUSIVE)"
            )

            resumen.append(
                f"Best matching subspecies: {dominante}"
            )

            resumen.append(
                f"Mapped reads: {dominante_reads}"
            )

            resumen.append(
                f"Genome coverage: {dominante_coverage:.2f}%"
            )

            resumen.append(
                f"Mean depth: {dominante_depth:.2f}X"
            )

            resumen.append(
                f"Mean identity: "
                f"{ordenados[0][1]['mean_identity']:.2f}%"
            )

            resumen.append(
                f"Reads >90% identity: "
                f"{ordenados[0][1]['reads90']}"
            )

            resumen.append(
                f"Reads >95% identity: "
                f"{ordenados[0][1]['reads95']}"
            )

            resumen.append("")
            resumen.append("Top matching references:")

            top_hits = ordenados[:3]

            for especie, datos in top_hits:

                resumen.append(
                    f" {especie}: "
                    f"{datos['reads']} reads, "
                    f"{datos['coverage']:.2f}% genome coverage, "
                    f"{datos['depth']:.2f}X, "
                    f"Mean identity={datos['mean_identity']:.2f}%, "
                    f"Reads >95% identity: {datos['reads95']}"
                )

        else:

            resumen.append(
                "Xylella fastidiosa: NEGATIVE"
            )

    else:

        resumen.append(
            "Xylella fastidiosa: NEGATIVE"
        )

# =====================================================
# XANTHOMONAS
# =====================================================

    xanthomonas_total = sum(
        v["reads"]
        for v in xanthomonas_results.values()
    )

    resumen.append("")
    resumen.append("=== XANTHOMONAS CITRI ===")

    if xanthomonas_total >= 5:

        ordenados = sorted(
            xanthomonas_results.items(),
            key=lambda x: x[1]["coverage"],
            reverse=True
        )

        dominante = ordenados[0][0]
        dominante_reads = ordenados[0][1]["reads"]
        dominante_coverage = ordenados[0][1]["coverage"]
        dominante_depth = ordenados[0][1]["depth"]
        dominante_identity = ordenados[0][1]["mean_identity"]
        dominante_reads90 = ordenados[0][1]["reads90"]
        dominante_reads95 = ordenados[0][1]["reads95"]

        segunda_coverage = 0

        if len(ordenados) > 1:
            segunda_coverage = ordenados[1][1]["coverage"]

        diferencia = dominante_coverage - segunda_coverage

        if dominante_reads >= 10 and dominante_coverage >= 0.1:

            resumen.append(
                "Xanthomonas citri: POSITIVE"
            )
            resumen.append(
                f"Best matching pathovar: {dominante}"
            )
            resumen.append(
                f"Mapped reads: {dominante_reads}"
            )

            resumen.append(
                f"Genome coverage: {dominante_coverage:.2f}%"
            )

            resumen.append(
                f"Mean depth: {dominante_depth:.2f}X"
            )

            resumen.append(
                f"Mean identity: "
                f"{ordenados[0][1]['mean_identity']:.2f}%"
            )

            resumen.append(
                f"xanthomonas total mapped reads: {xanthomonas_total}"
            )

            resumen.append(
                f"Dominant species reads95: {dominante_reads95}"
            )

            resumen.append(
                f"Reads >90% identity: "
                f"{ordenados[0][1]['reads90']}"
            )

            resumen.append(
                f"Reads >95% identity: "
                f"{ordenados[0][1]['reads95']}"
            )

            resumen.append(
                f"Difference from second-ranked pathovar: {diferencia:.2f}%"
            )

            resumen.append("")
            resumen.append("Top matching pathovars:")

            top_hits = ordenados[:3]

            for especie, datos in top_hits:

                resumen.append(
                    f" {especie}: "
                    f"{datos['reads']} reads, "
                    f"{datos['coverage']:.2f}% genome coverage, "
                    f"{datos['depth']:.2f}X, "
                    f"Mean identity={datos['mean_identity']:.2f}%, "
                    f"Reads >95% identity: {datos['reads95']}"
                )

        elif dominante_coverage > 0:

            resumen.append(
                "Xanthomonas citri: WEAK DETECTION (INCONCLUSIVE)"
            )

            resumen.append(
                f"Best matching pathovar: {dominante}"
            )
            resumen.append(
                f"Mapped reads: {dominante_reads}"
            )

            resumen.append(
                f"Genome coverage: {dominante_coverage:.2f}%"
            )

            resumen.append(
                f"Mean depth: {dominante_depth:.2f}X"
            )

            resumen.append(
                f"Mean identity: "
                f"{ordenados[0][1]['mean_identity']:.2f}%"
            )

            resumen.append(
                f"Reads >90% identity: "
                f"{ordenados[0][1]['reads90']}"
            )

            resumen.append(
                f"Reads >95% identity: "
                f"{ordenados[0][1]['reads95']}"
            )

            resumen.append(
                f"Difference from second-ranked pathovar: {diferencia:.2f}%"
            )

            resumen.append("")
            resumen.append("Top matching pathovars:")

            top_hits = ordenados[:3]

            for especie, datos in top_hits:

                resumen.append(
                    f" {especie}: "
                    f"{datos['reads']} reads, "
                    f"{datos['coverage']:.2f}% genome coverage, "
                    f"{datos['depth']:.2f}X, "
                    f"Mean identity={datos['mean_identity']:.2f}%, "
                    f"Reads >95% identity: {datos['reads95']}"
                )

        else:

            resumen.append(
                "Xanthomonas citri: NEGATIVE"
            )

    else:

        resumen.append(
            "Xanthomonas citri: NEGATIVE"
        )
    # =====================================================
    # CITRICARPA
    # =====================================================

    pcit_total = len(
        targets_reads.get(
            "P_citricarpa",
            set()
        )
    )

    resumen.append("")
    resumen.append("=== PHYLLOSTICTA CITRICARPA ===")

    if pcit_total >= 3:

        ordenados = sorted(
            pcit_results.items(),
            key=lambda x: x[1]["reads95"],
            reverse=True
        )

        dominante = ordenados[0][0]
        dominante_reads = ordenados[0][1]["reads"]
        dominante_coverage = ordenados[0][1]["coverage"]
        dominante_depth = ordenados[0][1]["depth"]
        dominante_identity = ordenados[0][1]["mean_identity"]
        dominante_reads90 = ordenados[0][1]["reads90"]
        dominante_reads95 = ordenados[0][1]["reads95"]
        dominante_reads97 = ordenados[0][1]["reads97"]
        dominante_reads99 = ordenados[0][1]["reads99"]
        ratio95 = 0
        ratio97 = 0
        ratio99 = 0

        if dominante_reads > 0:
            ratio95 = (dominante_reads95 / dominante_reads) * 100
            ratio97 = (dominante_reads97 / dominante_reads) * 100
            ratio99 = (dominante_reads99 / dominante_reads) * 100

        if dominante_reads >= 10 and dominante_coverage >= 0.1:

            resumen.append(
                "Phyllosticta citricarpa: POSITIVE"
            )

            resumen.append(
                f"Best matching reference: {dominante}"
            )

            resumen.append(
                f"Mapped reads: {dominante_reads}"
            )

            resumen.append(
                f"Genome coverage: {dominante_coverage:.2f}%"
            )

            resumen.append(
                f"Mean depth: {dominante_depth:.2f}X"
            )

            resumen.append(
                f"Mean identity: {dominante_identity:.2f}%"
            )

            resumen.append(
                f"Pcit total mapped reads: {pcit_total}"
            )

            resumen.append(
                f"Dominant species reads95: {dominante_reads95}"
            )

            resumen.append(
                f"Reads >90% identity: {dominante_reads90}"
            )

            resumen.append(
                f"Reads >95% identity: {dominante_reads95}"
            )

            resumen.append(
                f"Reads >97% identity: {dominante_reads97}"
            )

            resumen.append(
                f"Reads >99% identity: {dominante_reads99}"
            )
            resumen.append(
                f"95% identity ratio: {ratio95:.2f}%"
            )

            resumen.append(
                f"97% identity ratio: {ratio97:.2f}%"
            )

            resumen.append(
                f"99% identity ratio: {ratio99:.2f}%"
            )


            resumen.append("")
            resumen.append("Top matching references:")

            top_hits = ordenados[:3]

            for especie, datos in top_hits:

                resumen.append(
                    f" {especie}: "
                    f"{datos['reads']} reads, "
                    f"{datos['coverage']:.2f}% genome coverage, "
                    f"{datos['depth']:.2f}X, "
                    f"Mean identity={datos['mean_identity']:.2f}%, "
                    f"Reads >95% identity: {datos['reads95']}, "
                    f"Reads >99% identity: {datos['reads99']}"
                )

        elif dominante_coverage > 0:

            resumen.append(
                "Phyllosticta citricarpa: WEAK DETECTION (INCONCLUSIVE)"
            )

            resumen.append(
                f"Best matching reference: {dominante}"
            )

            resumen.append(
                f"Mapped reads: {dominante_reads}"
            )

            resumen.append(
                f"Genome coverage: {dominante_coverage:.2f}%"
            )

            resumen.append(
                f"Mean depth: {dominante_depth:.2f}X"
            )

            resumen.append(
                f"Mean identity: {dominante_identity:.2f}%"
            )

            resumen.append(
                f"Reads >90% identity: {dominante_reads90}"
            )

            resumen.append(
                f"Reads >95% identity: {dominante_reads95}"
            )

            resumen.append(
                f"Reads >97% identity: {dominante_reads97}"
            )

            resumen.append(
                f"Reads >99% identity: {dominante_reads99}"
            )

            resumen.append("")
            resumen.append("Top matching references:")

            top_hits = ordenados[:3]

            for especie, datos in top_hits:

                resumen.append(
                    f" {especie}: "
                    f"{datos['reads']} reads, "
                    f"{datos['coverage']:.2f}% genome coverage, "
                    f"{datos['depth']:.2f}X, "
                    f"Mean identity={datos['mean_identity']:.2f}%, "
                    f"Reads >95% identity: {datos['reads95']}"
                )


        else:

            resumen.append(
                "Phyllosticta citricarpa: NEGATIVE"
            )

    else:

        resumen.append(
            "Phyllosticta citricarpa: NEGATIVE"
        )

# -----------------------------------------
# SAVE SUMMARY
# -----------------------------------------
    resumen_file = output_dir / "Pathogen_summary.txt"

    with open(resumen_file, "w") as out:
        out.write("\n".join(resumen))

    content = "<br>".join(resumen)

    content = content.replace(
        "POSITIVE",
        '<span class="positive">POSITIVE</span>'
    )

    content = content.replace(
        "NEGATIVE",
        '<span class="negative">NEGATIVE</span>'
    )

    content = content.replace(
    "INCONCLUSIVE",
    '<span class="inconclusive">INCONCLUSIVE</span>'
    )

    html_file = output_dir / "Pathogen_summary.html"

    html_content = f"""
    <html>
    <head>
    <title>Pathogen Detection Report</title>
    <style>

    body {{
        font-family: Arial, sans-serif;
        margin: 30px;
        line-height: 1.5;
    }}

    .positive {{
        color: green;
        font-weight: bold;
    }}

    .negative {{
        color: red;
        font-weight: bold;
    }}

    .inconclusive {{
        color: orange;
        font-weight: bold;
    }}

    </style>
    </head>
    <body>

    <h1>Pathogen Detection Summary</h1>

    {content}

    </body>
    </html>
    """

    with open(html_file, "w") as out:
        out.write(html_content)

    log.info(
        f"Summary saved in {resumen_file}"
    )

    log.info(
        f"HTML report saved in {html_file}"
    )

def write_qc_report(
    output_dir,
    raw_reads,
    mean_q_raw,
    mean_q_filtered,
    filtered_reads,
    targets_reads,
    liberibacter_results,
    xylella_results,
    xanthomonas_results,
    pcit_results
):

        report = output_dir / "QC_summary.txt"

        with open(report, "w") as out:

            out.write("===== PIPELINE QC SUMMARY =====\n\n")

            out.write("RAW DATA\n")
            out.write("--------\n")
            out.write(f"Total reads: {raw_reads}\n")
            out.write(f"Mean Q-score: {mean_q_raw}\n\n")

            out.write("AFTER FASTPLONG\n")
            out.write("----------------\n")
            out.write(f"Reads retained: {filtered_reads}\n")
            out.write(f"Reads removed: {raw_reads-filtered_reads}\n")
            out.write(f"Mean Q-score: {mean_q_filtered}\n")

            if raw_reads > 0:
                retention = (
                    filtered_reads/raw_reads
                ) * 100

                out.write(
                    f"Retention: {retention:.2f}%\n"
                )

            out.write("\n")

            out.write("KRAKEN CLASSIFICATION\n")
            out.write("---------------------\n")

            for pathogen, reads in targets_reads.items():
                out.write(
                    f"{pathogen}: {len(reads)} reads\n"
                )

            out.write("\n")

            datasets = [
                ("Liberibacter", liberibacter_results),
                ("Xylella", xylella_results),
                ("Xanthomonas", xanthomonas_results),
                ("P_citricarpa", pcit_results)
            ]

            out.write("VALIDATION RESULTS\n")
            out.write("------------------\n")

            for group_name, results in datasets:

                if not results:
                    continue

                best = max(
                    results.items(),
                    key=lambda x: x[1]["reads"]
                )

                ref = best[0]
                data = best[1]

                out.write(f"\n{group_name}\n")

                out.write(
                    f"Best reference: {ref}\n"
                )

                out.write(
                f"Mapped reads: {data['reads']}\n"
                )

                out.write(
                f"Reads >80% identity: {data['reads80']}\n"
                )

                out.write(
                f"Reads >85% identity: {data['reads85']}\n"
                )

                out.write(
                f"Reads >90% identity: {data['reads90']}\n"
                )

                out.write(
                f"Reads >95% identity: {data['reads95']}\n"
                )
                out.write(
                f"Reads >97% identity: {data['reads97']}\n"
                )

                out.write(
                f"Reads >99% identity: {data['reads99']}\n"
                )

                out.write(
                f"Reads with identity <95%: "
                f"{data['reads'] - data['reads95']}\n"
                )

# ---------------- EXTRAER READS ----------------

def extract_targets(
    kraken_output,
    fastq_file,
    output_dir,
    kraken_db
):

    log.info("Extracting organisms of interest...")

    read_taxid = {}

    with open(kraken_output) as f:

        for line in f:

            parts = line.strip().split("\t")

            if len(parts) < 3:
                continue

            read_id = parts[1]
            taxid = parts[2]

            read_taxid[read_id] = taxid

    targets_reads = defaultdict(set)

    for read, taxid in read_taxid.items():

# Asignaciones específicas a especie
        for name in [
            "CaLas",
            "CaLaf",
            "CaLam",
            "CaLsol",
            "Liberibacter_crescens"
        ]:
            if taxid in TARGETS[name]:
                targets_reads[name].add(read)

# Toda lectura asignada a una especie de Liberibacter
# también pasa al conjunto general de Liberibacter.
                targets_reads["Liberibacter_genus"].add(read)

# Asignación directamente al género Liberibacter
        if taxid in TARGETS["Liberibacter_genus"]:
            targets_reads["Liberibacter_genus"].add(read)

# Asignaciones de los demás patógenos
        for name in [
            "Xanthomonas_citri",
            "Xylella_fastidiosa",
            "P_citricarpa"
        ]:
            if taxid in TARGETS[name]:
                targets_reads[name].add(read)

    with open(fastq_file) as fq:
        lines = fq.readlines()

    for name, reads in targets_reads.items():

        if not reads:
            continue

        out_file = output_dir / f"{name}.fastq"

        with open(out_file, "w") as out:

            for i in range(0, len(lines), 4):

                read_id = lines[i].split()[0][1:]

                if read_id in reads:
                    out.writelines(lines[i:i+4])

        log.info(
            f"{name}: {len(reads)} reads"
        )

    liberibacter_reads = set()

    for especie in [
        "CaLas",
        "CaLaf",
        "CaLam",
        "CaLsol",
        "Liberibacter_crescens",
        "Liberibacter_genus"
    ]:
        liberibacter_reads.update(
            targets_reads.get(especie, set())
        )

    log.info(
        f"Liberibacter candidate reads: {len(liberibacter_reads)}"
    )

    liberibacter_fastq = output_dir / "Liberibacter.fastq"

    with open(fastq_file) as fq:
        lines = fq.readlines()

    with open(liberibacter_fastq, "w") as out:

        for i in range(0, len(lines), 4):

            read_id = lines[i].split()[0][1:]

            if read_id in liberibacter_reads:
                out.writelines(lines[i:i+4])

    liberibacter_results = map_against_references(
       output_dir / "Liberibacter.fastq",
        LIBERIBACTER_REFS,
        output_dir,
        "liberibacter",
        kraken_db
    )

    xylella_results = map_against_references(
        output_dir / "Xylella_fastidiosa.fastq",
        XYLELLA_REFS,
        output_dir,
        "xylella",
        kraken_db
    )

    xanthomonas_results = map_against_references(
        output_dir / "Xanthomonas_citri.fastq",
        XANTHOMONAS_REFS,
        output_dir,
        "xanthomonas",
        kraken_db
    )

    pcit_results = map_against_references(
        output_dir / "P_citricarpa.fastq",
        PCIT_REFS,
        output_dir,
        "pcit",
        kraken_db
    )

    return (
        targets_reads,
        liberibacter_results,
        xylella_results,
        xanthomonas_results,
        pcit_results
    )
# ---------------- PROCESAR MUESTRA ----------------

def process_sample(
    input_fastq,
    output_dir,
    args
):

    filtered_out = output_dir / "filtered.fastq"
    raw_reads = count_fastq_reads(input_fastq)
    mean_q_raw = calculate_mean_qscore(
        input_fastq
    )
    kraken_report = output_dir / "kraken.report"
    kraken_output = output_dir / "kraken.out"

    run_cmd([
        "fastplong",
        "-i", str(input_fastq),
        "-o", str(filtered_out),
        "-q", "10",
        "-u", "60",
        "-h", str(output_dir / "fastplong.html"),
        "-j", str(output_dir / "fastplong.json"),
        "--thread", str(args.threads)
    ])

    if (
        not filtered_out.exists()
        or filtered_out.stat().st_size == 0
    ):
        log.warning(
            f"{output_dir.name}: no readings after filtering, skipping"
        )
        return

    run_cmd([
        "kraken2",
        "--db", args.kraken_db,
        "--threads", str(args.threads),
        "--confidence", "0.05",
        "--report", str(kraken_report),
        "--output", str(kraken_output),
        str(filtered_out)
    ])
    filtered_reads = count_fastq_reads(filtered_out)

    mean_q_filtered = calculate_mean_qscore(filtered_out)

    if not kraken_output.exists():

        log.warning(
            f"{output_dir.name}: Kraken yielded no results, skipping."
        )

        return

    (
        targets_reads,
        liberibacter_results,
        xylella_results,
        xanthomonas_results,
        pcit_results
    ) = extract_targets(
        kraken_output,
        filtered_out,
        output_dir,
        args.kraken_db
    )

    filtered_reads = count_fastq_reads(
        filtered_out
    )

    write_qc_report(
        output_dir,
        raw_reads,
        mean_q_raw,
        mean_q_filtered,
        filtered_reads,
        targets_reads,
        liberibacter_results,
        xylella_results,
        xanthomonas_results,
        pcit_results
    )

    interpret_pathogens(
        targets_reads,
        output_dir,
        liberibacter_results,
        xylella_results,
        xanthomonas_results,
        pcit_results
    )
# ---------------- MAIN ----------------

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--run-dir",
        required=True
    )

    parser.add_argument(
        "--work-dir",
        required=True
    )

    parser.add_argument(
        "--kraken-db",
        required=True
    )

    parser.add_argument(
        "--threads",
        type=int,
        default=4
    )

    args = parser.parse_args()

    check_dependencies()

    run_dir = Path(args.run_dir)
    work_dir = Path(args.work_dir)

    work_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    barcodes = find_barcodes(run_dir)

    if barcodes:

        log.info(
            f"Found {len(barcodes)} barcodes"
        )

        for bc in barcodes:

            fastqs = (
                list(bc.rglob("*.fastq")) +
                list(bc.rglob("*.fastq.gz"))
            )

            if not fastqs:
                continue

            out_dir = work_dir / bc.name

            out_dir.mkdir(
                exist_ok=True
            )

            merged = out_dir / "merged.fastq"

            merge_fastqs(
                fastqs,
                merged
            )

            process_sample(
                merged,
                out_dir,
                args
            )

    else:

        log.info(
            "No barcodes detected. Processing as a single sample."
        )

        fastqs = (
            list(run_dir.rglob("*.fastq")) +
            list(run_dir.rglob("*.fastq.gz"))
        )

        merged = work_dir / "merged.fastq"

        merge_fastqs(
            fastqs,
            merged
        )

        process_sample(
            merged,
            work_dir,
            args
        )

if __name__ == "__main__":
    main()

# GenomeRadar

**GenomeRadar** is a Python tool for exploring the genomic context surrounding specific positions in **GenBank-annotated genomes**. It extracts and annotates flanking genes and other genomic features, enabling systematic analysis of local genomic neighborhoods around positions of biological interest.

---

## 🧭 What does GenomeRadar do?

Given:

* 📁 A folder containing **GenBank files** (`.gb`, `.gbff`, optionally compressed as `.gz`)
* 📄 A **TSV file** listing genomic positions of interest for each file

GenomeRadar will:

* Identify the genomic feature **containing the target position** (if any)
* Extract up to **20 upstream and 20 downstream features** around each position
* Filter features by biologically relevant annotation types
* Compute **relative distances and ordering** with respect to the target position
* Produce a **tab-delimited output table** ready for downstream analysis

---

## 🧬 Supported feature types

By default, GenomeRadar processes the following GenBank feature types:

```
CDS, tRNA, operon, regulatory, mobile_element, ncRNA, mRNA, rRNA,
oriT, protein_bind, repeat_region, tmRNA
```

These can be easily customized in the source code via the `ALLOWED_TYPES` variable.

---

## 📦 Requirements

* Python ≥ 3.8
* Biopython

Install dependencies with:

```bash
pip install biopython
```

---

## 📁 Positions file format

A **TSV file** with one position per line:

```
example_genome.gb	123456
example_genome.gb	789012
```

* Column 1: GenBank filename
* Column 2: Genomic position (integer)

---

## ▶️ Usage

```bash
python genomeradar.py \
  -i genbank_folder/ \
  -p positions.tsv \
  -o output.tsv
```

### Arguments

| Option            | Description                     |
| ----------------- | ------------------------------- |
| `-i, --input`     | Folder containing GenBank files |
| `-p, --positions` | TSV file with genomic positions |
| `-o, --output`    | Output file                     |

---

## 📊 Output format

The output is a TSV file containing:

* Contig ID
* GenBank file name
* Organism
* Taxonomic lineage
* Target position
* Locus tag
* Feature start and end coordinates
* Feature type
* Product / protein description
* Distance to the target position
* Relative order (upstream / downstream)
* Strand (+ / −)

---

## 🧠 Typical use cases

* Analysis of genomic neighborhoods around specific loci
* Context exploration for mutations, insertion sites, or breakpoints
* Comparative genomics and annotation inspection
* Data preparation for downstream statistical or visualization analyses

---

## 🛠️ Customization

* Adjust the number of flanking features using `upstream_count` and `downstream_count`
* Modify the set of allowed feature types in `ALLOWED_TYPES`

---

## 📜 License

This project is distributed under the **MIT License**.

---

## ✨ Project name

**GenomeRadar** reflects the idea of scanning genomic regions around specific coordinates to detect structure, context, and biologically relevant signals.

---

If you would like, I can also:

* Add a workflow or schematic diagram
* Create a `requirements.txt` file
* Adapt the README for a publication or supplementary material
* Provide example input and output files

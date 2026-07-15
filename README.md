# Deduper — Reference-Based PCR Duplicate Removal

A memory-efficient, single-pass tool for removing PCR duplicates from sorted single-end SAM files using known UMIs, without loading the entire file into memory.

## Problem

PCR-based library prep introduces duplicate reads that inflate coverage and bias downstream quantification. Standard duplicate marking (e.g. `samtools markdup`) uses position alone; this tool additionally uses the **UMI embedded in the read name** to distinguish true PCR duplicates from reads that independently mapped to the same position — a common and more rigorous strategy in real variant-calling and single-cell pipelines.

## Approach

Given a coordinate-sorted SAM file of uniquely-mapped single-end reads and a list of 96 known UMIs:

1. Stream the SAM file one read at a time (no full-file load — designed to scale to files with tens of millions of reads).
2. Parse the CIGAR string to compute the read's **true 5' start position**, correcting for soft-clipping on both strands (a duplicate can appear to start at different positions in the SAM file if clipping differs, even though the underlying fragment is identical).
3. Determine strand from the bitwise flag and combine `(chromosome, corrected position, strand, UMI)` into a duplicate-detection key.
4. Discard reads with unrecognized UMIs (sequencing/PCR errors on the UMI itself), and retain only the first occurrence of each key.
5. Write a properly formatted, deduplicated SAM file.

Design and test cases were worked out in pseudocode before implementation — see [`Test_files/Psuedocode.md`](Test_files/Psuedocode.md) and the paired [`Test_files/Input.sam`](Test_files/Input.sam) / [`Test_files/Output.sam`](Test_files/Output.sam) fixtures used to validate correctness.

## Usage

```bash
./Girard_deduper.py -f <sorted_input.sam> -o <deduplicated_output.sam> -u STL96.txt
```

| Flag | Description |
|---|---|
| `-f`, `--file` | Path to a coordinate-sorted SAM file (single-end, uniquely mapped reads) |
| `-o`, `--outfile` | Path to write the deduplicated SAM file |
| `-u`, `--umi` | Text file of known UMIs, one per line ([`STL96.txt`](STL96.txt)) |

**Note:** input must be sorted first, e.g. `samtools sort input.sam -o sorted.sam`.

## Skills demonstrated

Algorithm design from a written spec · CIGAR string parsing · SAM bitwise flag interpretation · streaming/constant-memory file processing · `argparse` CLI design · test-driven validation with hand-crafted fixtures

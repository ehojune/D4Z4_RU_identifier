#!/usr/bin/env python3

import os
import sys
import argparse
import configparser
from dataclasses import dataclass
from operator import attrgetter
from typing import List, Tuple
import matplotlib.pyplot as plt

@dataclass
class BlastResult:
    """Class for storing BLAST results"""
    qseqid: str
    sseqid: str
    pident: float
    length: int
    mismatch: int
    gapopen: int
    qstart: int
    qend: int
    sstart: int
    send: int
    evalue: str
    bitscore: float
    qlen: int
    slen: int

@dataclass
class D4Z4:
    """Class for storing D4Z4 repeat information"""
    pos: Tuple[int, int]
    XapI: str  # Position or 'NA'
    BlnI: str  # Position or 'NA'

@dataclass
class MatchedRead:
    """Class for storing information about reads containing D4Z4 repeats"""
    pos: Tuple[int, int]
    d4z4: List[D4Z4]
    pLAM: Tuple[int, int]
    probe: Tuple[int, int]
    pLAM_seq: str
    read_id: str

def convert_fq2fa(fastq_path: str, fasta_path: str) -> None:
    """Convert FASTQ to FASTA format"""
    command = f"cat {fastq_path} | awk '{{if(NR%4==1) {{printf(\">%s\\n\",substr($0,2));}} else if(NR%4==2) print;}}\' > {fasta_path}"
    os.system(command)

def run_makeblastdb(makeblastdb_path: str, input_fasta: str, output_db: str) -> None:
    """Create BLAST database"""
    os.makedirs(output_db, exist_ok=True)
    command = f"{makeblastdb_path} -in {input_fasta} -dbtype nucl -out {output_db}"
    os.system(command)

def run_blastn(blastn_path: str, query: str, db: str, output: str) -> None:
    """Run BLAST with standard parameters"""
    command = (f"{blastn_path} -query {query} -db {db} -out {output} "
              "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\"")
    os.system(command)

def run_blastn_short(blastn_path: str, query: str, db: str, output: str, word_size: int = None) -> None:
    """Run BLAST with parameters optimized for short sequences"""
    base_command = (f"{blastn_path} -task blastn-short -dust no -evalue 1000 -soft_masking false "
                   f"-query {query} -db {db}")
    if word_size:
        base_command += f" -word_size {word_size}"

    # Run with default output format
    os.system(f"{base_command} -out {output}")
    
    # Run with tabular output format
    os.system(f"{base_command} -out {output}.outfmt6.tsv -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\"")

def parse_blast_results(blast_file: str) -> List[BlastResult]:
    """Parse BLAST results from tabular output"""
    results = []
    with open(blast_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.strip().split('\t')
            results.append(BlastResult(
                fields[0], fields[1], float(fields[2]), int(fields[3]),
                int(fields[4]), int(fields[5]), int(fields[6]), int(fields[7]),
                int(fields[8]), int(fields[9]), fields[10], float(fields[11]),
                int(fields[12]), int(fields[13])
            ))
    return results

def find_overlapping_reads(blast_file1: str, blast_file2: str) -> List[str]:
    """Find reads that appear in both BLAST result files"""
    results1 = parse_blast_results(blast_file1)
    results2 = parse_blast_results(blast_file2)
    
    reads1 = {result.sseqid for result in results1}
    reads2 = {result.sseqid for result in results2}
    
    return list(reads1 & reads2)

def filter_blast_results(results: List[BlastResult], min_coverage: float, min_identity: float) -> List[BlastResult]:
    """Filter BLAST results based on coverage and identity thresholds"""
    return [r for r in results 
            if r.length >= r.qlen * min_coverage and r.pident >= min_identity * 100]

def extract_sequence(fasta_path: str, read_id: str, output_path: str) -> None:
    """Extract a specific sequence from a FASTA file"""
    with open(fasta_path) as fin, open(output_path, 'w') as fout:
        write_seq = False
        for line in fin:
            if line.startswith('>'):
                if line[1:].strip() == read_id:
                    write_seq = True
                    fout.write(line)
                else:
                    write_seq = False
            elif write_seq:
                fout.write(line)
                break

def visualize_read(matched_read: MatchedRead, output_path: str) -> None:
    """Create visualization of D4Z4 repeats, restriction sites, and other features"""
    fig, ax = plt.subplots(figsize=(10, 2))
    direction = 1 if matched_read.pos[0] < matched_read.pos[1] else -1
    
    # Draw main read arrow
    ax.annotate("", xy=(matched_read.pos[1], 0), xytext=(matched_read.pos[0], 0),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))
    
    # Add read ID and positions
    ax.text((matched_read.pos[0] + matched_read.pos[1]) / 2, -0.03, 
            matched_read.read_id, ha='center', va='bottom', color='black', fontsize=8)
    
    # Add pLAM and probe features
    y_pos = 0.03
    for feature, label, color in [(matched_read.pLAM, 'pLAM', 'red'),
                                 (matched_read.probe, 'probe', 'red')]:
        ax.annotate("", xy=(feature[1], y_pos), xytext=(feature[0], y_pos),
                    arrowprops=dict(arrowstyle='->', color=color, lw=2 * direction))
        ax.text((feature[0] + feature[1]) / 2, y_pos - 0.01, 
                label, ha='center', va='top', color=color, fontsize=8)
    
    # Add D4Z4 repeats and restriction sites
    for i, d4z4 in enumerate(matched_read.d4z4):
        # Draw restriction sites
        for site, pos, color in [(d4z4.BlnI, 'BlnI', 'blue'),
                                (d4z4.XapI, 'XapI', 'green')]:
            if site != 'NA':
                ax.plot([float(site), float(site)], 
                       [y_pos + 0.01, y_pos + 0.03], 
                       color=color, lw=2 * direction)
                ax.text(float(site), y_pos + 0.03, pos, 
                       ha='center', va='bottom', color=color, fontsize=6)
        
        # Draw D4Z4 repeat
        ax.annotate("", xy=(d4z4.pos[1], y_pos), xytext=(d4z4.pos[0], y_pos),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2 * direction))
        ax.text((d4z4.pos[0] + d4z4.pos[1]) / 2, y_pos - 0.01,
                f"repeat_{i + 1}", ha='center', va='top', color='red', fontsize=8)
    
    ax.set_ylim(-0.1, 0.1)
    ax.axis('off')
    plt.savefig(output_path, format='png', dpi=1200, bbox_inches='tight')
    plt.close()

def main():
    parser = argparse.ArgumentParser(description='Analyze D4Z4 repeats for FSHD1 diagnosis')
    parser.add_argument('--config', required=True, help='Path to configuration file')
    parser.add_argument('--output', required=True, help='Output directory path')
    parser.add_argument('--min-coverage', type=float, default=0.8,
                        help='Minimum alignment coverage (default: 0.8)')
    parser.add_argument('--min-identity', type=float, default=0.8,
                        help='Minimum alignment identity (default: 0.8)')
    args = parser.parse_args()

    # Create output directories
    os.makedirs(args.output, exist_ok=True)
    temp_dir = os.path.join(args.output, 'tmp')
    os.makedirs(temp_dir, exist_ok=True)

    # Read configuration
    config = configparser.ConfigParser()
    config.read(args.config)
    
    # Process each sample
    for sample_name, fastq_path in config['Sample FASTQ Path'].items():
        print(f"Processing sample: {sample_name}")
        
        # Set up directories for this sample
        sample_dir = os.path.join(args.output, sample_name)
        os.makedirs(sample_dir, exist_ok=True)
        
        # Convert FASTQ to FASTA
        fasta_path = os.path.join(sample_dir, f"{sample_name}.fasta")
        convert_fq2fa(fastq_path, fasta_path)
        
        # Create BLAST database
        run_makeblastdb(config['External Programs']['makeblastdb_path'],
                       fasta_path,
                       os.path.join(temp_dir, sample_name))
        
        # Run BLAST searches
        blast_dir = os.path.join(sample_dir, 'blast')
        os.makedirs(blast_dir, exist_ok=True)
        
        for target in ['pLAM', 'probe', 'D4Z4']:
            run_blastn(config['External Programs']['blastn_path'],
                      config['Reference Sequences Path'][f'{target}_ref'],
                      os.path.join(temp_dir, sample_name),
                      os.path.join(blast_dir, f"{target}.blast.txt"))
        
        # Find reads containing both pLAM and probe
        overlapping_reads = find_overlapping_reads(
            os.path.join(blast_dir, "pLAM.blast.txt"),
            os.path.join(blast_dir, "probe.blast.txt")
        )
        
        if not overlapping_reads:
            print(f"No reads containing both pLAM and probe found in sample {sample_name}")
            continue
        
        # Process each overlapping read
        for read_id in overlapping_reads:
            print(f"Processing read: {read_id}")
            read_dir = os.path.join(sample_dir, read_id)
            os.makedirs(read_dir, exist_ok=True)
            
            # Extract read sequence
            read_fasta = os.path.join(read_dir, f"{read_id}.fasta")
            extract_sequence(fasta_path, read_id, read_fasta)
            
            # Run additional analyses and create visualization
            # (rest of the read-specific analysis steps...)
            
            # Note: Detailed implementation of read-specific analysis steps 
            # (restriction site mapping, D4Z4 repeat counting, etc.) continues here...

if __name__ == "__main__":
    main()

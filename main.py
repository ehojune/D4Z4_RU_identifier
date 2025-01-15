import os
import sys
import argparse
import configparser
from dataclasses import dataclass
from operator import attrgetter
from typing import List, Tuple, Set
import matplotlib.pyplot as plt


# Hardcoded sequences
BLNI_SEQ = "CCTAGG"
XAPI_SEQ = "AAATTC"
PLAM_4QA_SEQ = "CAGAGATATATTAAAATGCCCCCTCCC"
PLAM_10Q_SEQ = "CAGAGATATATCAAAATGCCCCCTCCC"

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
    print(command)
    os.system(command)

def run_makeblastdb(makeblastdb_path: str, input_fasta: str, output_db: str) -> None:
    """Create BLAST database"""
    command = f"{makeblastdb_path} -in {input_fasta} -dbtype nucl -out {output_db}"
    os.makedirs(os.path.dirname(output_db), exist_ok=True)
    os.system(command)

def extract_read_from_fasta_with_readID(readID: str, fasta_path: str, output_path: str) -> int:
    """Extract a specific read from a FASTA file"""
    with open(fasta_path, 'r') as fasta, open(output_path, 'w') as fw:
        while True:
            line = fasta.readline()
            if len(line) < 2:
                return 0
            if line.startswith(f">{readID}"):
                fw.write(line)
                line = fasta.readline()
                fw.write(line)
                print(f"write {output_path}.")
                return 0
            else:
                fasta.readline()

def run_blastn_short(blastn_path: str, query: str, db: str, output: str, word_size: int = None) -> None:
    """Run BLAST with parameters optimized for short sequences"""
    base_command = (f"{blastn_path} -task blastn-short -dust no -evalue 1000"
                   f"-soft_masking false -query {query} -db {db}")
    
    if word_size:
        base_command += f" -word_size {word_size}"
    
    # Run with default output format
    os.system(f"{base_command} -out {output}")
    
    # Run with tabular output format
    os.system(f"{base_command} -out {output}.outfmt6.tsv "
             "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\"")

def run_blastn(blastn_path: str, query: str, db: str, output: str) -> None:
    """Run BLAST with standard parameters"""
    command = (f"{blastn_path} -query {query} -db {db} -out {output} "
              "-outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\"")
    os.system(command)

def browse_blastresult(blast_filename: str) -> List[BlastResult]:
    """Parse BLAST results from file"""
    results = []
    with open(blast_filename, 'r') as fr:
        for line in fr:
            if len(line) < 2 or line.startswith('#') or line.startswith('qseqid'):
                continue
            fields = line.strip().split('\t')
            results.append(BlastResult(
                fields[0], fields[1], float(fields[2]), int(fields[3]),
                int(fields[4]), int(fields[5]), int(fields[6]), int(fields[7]),
                int(fields[8]), int(fields[9]), fields[10], float(fields[11]),
                int(fields[12]), int(fields[13])
            ))
    return results

def filter_blastresult_with_seqid(blastresult_list: List[BlastResult], seqid: str) -> List[BlastResult]:
    """Filter BLAST results for a specific sequence ID"""
    return [result for result in blastresult_list if result.sseqid == seqid]

def write_blastresult_tsv(blastresult_list: List[BlastResult], save_path: str) -> None:
    """Write BLAST results to TSV file"""
    with open(save_path, 'w') as fw:
        fw.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\n")
        for result in blastresult_list:
            fw.write(f"{result.qseqid}\t{result.sseqid}\t{result.pident}\t{result.length}\t"
                    f"{result.mismatch}\t{result.gapopen}\t{result.qstart}\t{result.qend}\t"
                    f"{result.sstart}\t{result.send}\t{result.evalue}\t{result.bitscore}\t"
                    f"{result.qlen}\t{result.slen}\n")

def draw_arrows(matched_read: MatchedRead, save_path: str) -> None:
    """Create visualization of read features"""
    fig, ax = plt.subplots(figsize=(10, 2))
    direction = 1 if matched_read.pos[0] < matched_read.pos[1] else -1
    
    # Draw main read arrow
    ax.annotate("", xy=(matched_read.pos[1], 0), xytext=(matched_read.pos[0], 0),
                arrowprops=dict(arrowstyle='->', color='black', lw=2))
    
    # Add read ID and positions
    ax.text((matched_read.pos[0] + matched_read.pos[1]) / 2, -0.03,
            matched_read.read_id, ha='center', va='bottom', color='black', fontsize=8)
    ax.text(matched_read.pos[0], 0.01, str(matched_read.pos[0]),
            ha='center', va='bottom', color='black', fontsize=6)
    ax.text(matched_read.pos[1], 0.01, str(matched_read.pos[1]),
            ha='center', va='bottom', color='black', fontsize=6)
    
    # Draw pLAM and probe
    index_y = 0.03
    for feature, label in [(matched_read.pLAM, 'pLAM'), (matched_read.probe, 'probe')]:
        ax.annotate("", xy=(feature[1], index_y), xytext=(feature[0], index_y),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2 * direction))
        ax.text((feature[0] + feature[1]) / 2, index_y - 0.01,
                label, ha='center', va='top', color='red', fontsize=8)
    
    # Draw D4Z4 repeats and restriction sites
    for i, d4z4 in enumerate(matched_read.d4z4):
        if d4z4.BlnI != 'NA':
            ax.plot([float(d4z4.BlnI), float(d4z4.BlnI)],
                   [index_y + 0.01, index_y + 0.03],
                   color='blue', lw=2 * direction)
            ax.text(float(d4z4.BlnI), index_y + 0.03,
                   'BlnI', ha='center', va='bottom', color='blue', fontsize=6)
        
        if d4z4.XapI != 'NA':
            ax.plot([float(d4z4.XapI), float(d4z4.XapI)],
                   [index_y + 0.01, index_y + 0.03],
                   color='green', lw=2 * direction)
            ax.text(float(d4z4.XapI), index_y + 0.03,
                   'XapI', ha='center', va='bottom', color='green', fontsize=6)
        
        # Draw D4Z4 repeat
        ax.annotate("", xy=(d4z4.pos[1], index_y), xytext=(d4z4.pos[0], index_y),
                    arrowprops=dict(arrowstyle='->', color='red', lw=2 * direction))
        
        # Add position labels and repeat number
        if direction == -1:
            ax.text(d4z4.pos[0], index_y + 0.02, f"{d4z4.pos[0]}",
                   ha='right', va='bottom', color='red', fontsize=6)
            ax.text(d4z4.pos[1], index_y + 0.02, f"{d4z4.pos[1]}",
                   ha='left', va='bottom', color='red', fontsize=6)
        else:
            ax.text(d4z4.pos[0], index_y + 0.02, f"{d4z4.pos[0]}",
                   ha='left', va='bottom', color='red', fontsize=6)
            ax.text(d4z4.pos[1], index_y + 0.02, f"{d4z4.pos[1]}",
                   ha='right', va='bottom', color='red', fontsize=6)
        
        ax.text((d4z4.pos[0] + d4z4.pos[1]) / 2, index_y - 0.01,
                f"repeat_{i + 1}", ha='center', va='top', color='red', fontsize=8)
    
    # Set plot limits and style
    if direction == -1:
        ax.set_xlim(0, max(matched_read.pos[0],
                          max(d4z4.pos[0] for d4z4 in matched_read.d4z4)) + 1000)
    else:
        ax.set_xlim(0, max(matched_read.pos[1],
                          max(d4z4.pos[1] for d4z4 in matched_read.d4z4)) + 1000)
    
    ax.set_ylim(-0.1, 0.1)
    ax.axis('off')
    plt.savefig(save_path, format='png', dpi=1200, bbox_inches='tight')
    plt.close()

def check_overlap_seqid_in_blastresults(blast_filename1: str, blast_filename2: str) -> List[str]:
    """Find overlapping sequence IDs between two BLAST result files"""
    blastresult1_list = browse_blastresult(blast_filename1)
    blastresult2_list = browse_blastresult(blast_filename2)
    
    seqids1 = {result.sseqid for result in blastresult1_list}
    seqids2 = {result.sseqid for result in blastresult2_list}
    
    return list(seqids1 & seqids2)
    
def extract_tr_from_read(readID: str, read_fasta_path: str, tr_read_blast_path: str, output_path: str):
    """Extract tandem repeats from read"""
    # Get TR positions
    tr_positions = []
    with open(tr_read_blast_path, 'r') as fr:
        for line in fr:
            if line.startswith("qseqid"): continue
            tr_pos = (int(line.split()[8])-1, int(line.split()[9])-1)
            if len(line)>2: tr_positions.append(tr_pos)
            else: break

    # Get sequence
    with open(read_fasta_path, 'r') as fr:
        fr.readline()  # skip header
        seq = fr.readline().strip()

    tr_list = []
    if tr_positions[0][0] > tr_positions[0][1]:  # reverse read
        for tr_pos in tr_positions:
            adj_tr_pos = (tr_pos[1], tr_pos[0]+1)
            tr = revcomp(seq[adj_tr_pos[0]:adj_tr_pos[1]])
            tr_list.append(tr)
        tr_list.reverse()
    else:  # forward read
        for tr_pos in tr_positions:
            adj_tr_pos = (tr_pos[0], tr_pos[1]+1)
            tr = seq[adj_tr_pos[0]:adj_tr_pos[1]]
            tr_list.append(tr)

    # Write TR units
    with open(output_path, 'w') as fw:
        for tr_idx, tr_seq in enumerate(tr_list):
            fw.write(f">{readID}_repeat{tr_idx+1}\n{tr_seq}\n")
            
def revcomp(seq: str) -> str:
    """Return reverse complement of a sequence"""
    comp_map = {"A": "T", "G": "C", "C": "G", "T": "A"}
    return "".join(comp_map[base] for base in seq.upper()[::-1])


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
    
    # Create fasta directory for query sequences
    fasta_dir = os.path.join(args.output, 'fasta')
    os.makedirs(fasta_dir, exist_ok=True)
    
     # Create fasta directory and sequence files
    fasta_dir = os.path.join(args.output, 'fasta')
    os.makedirs(fasta_dir, exist_ok=True)
    
    # Create sequence files with hardcoded sequences
    sequences = [
        ('BlnI_seq', BLNI_SEQ),
        ('XapI_seq', XAPI_SEQ),
        ('pLAM_4qA_seq', PLAM_4QA_SEQ),
        ('pLAM_10q_seq', PLAM_10Q_SEQ)
    ]
    
    for seq_type, seq in sequences:
        with open(os.path.join(fasta_dir, f"{seq_type}.fasta"), 'w') as f:
            f.write(f">{seq_type}\n{seq}\n")
    

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
        overlapping_reads = check_overlap_seqid_in_blastresults(
            os.path.join(blast_dir, "pLAM.blast.txt"),
            os.path.join(blast_dir, "probe.blast.txt")
        )
        
        if not overlapping_reads:
            print(f"No reads containing both pLAM and probe found in sample {sample_name}")
            continue
        
        # Process each overlapping read
        for read_id in overlapping_reads:
            print(f"Processing read: {read_id}")
            read_dir = os.path.join(blast_dir, read_id)
            os.makedirs(read_dir, exist_ok=True)
            
            # Collect all BLAST results for this read
            all_results = []
            for target in ['pLAM', 'probe', 'D4Z4']:
                blast_file = os.path.join(blast_dir, f"{target}.blast.txt")
                read_results = filter_blastresult_with_seqid(
                    browse_blastresult(blast_file), read_id)
                all_results.extend(read_results)
            
            # Write combined results
            write_blastresult_tsv(all_results, 
                                os.path.join(read_dir, f"{read_id}.blastn.tsv"))
            
            # Extract individual BLAST results for D4Z4 repeats
            D4Z4_blastresult = [result for result in all_results if result.qseqid == "D4Z4"]
            write_blastresult_tsv(D4Z4_blastresult, 
                                os.path.join(read_dir, f"{read_id}.D4Z4.blastn.tsv"))
            
            # Process BLAST results for visualization
            read_direction = 1 if all_results[0].sstart < all_results[0].send else -1
            read_length = all_results[0].slen
            
            # Find pLAM and probe positions
            pLAM_pos = (0, 0)
            probe_pos = (0, 0)
            for result in all_results:
                if result.qseqid == "pLAM":
                    pLAM_pos = (result.sstart, result.send)
                elif result.qseqid == "GU550601.1":  # probe ID
                    probe_pos = (result.sstart, result.send)
            
            # Process D4Z4 repeats
            read_d4z4_list = []
            for result in all_results:
                if result.qseqid == "D4Z4":
                    d4z4_pos = D4Z4((result.sstart, result.send), 'NA', 'NA')
                    read_d4z4_list.append(d4z4_pos)
            
            if read_direction == -1:
                read_d4z4_list = read_d4z4_list[::-1]
                
                
            # Create read directory and run BLAST-short
            read_dir = os.path.join(blast_dir, read_id)
            os.makedirs(read_dir, exist_ok=True)

            # Extract read sequence to a separate FASTA file
            read_fasta = os.path.join(read_dir, f"{read_id}.fasta")
            extract_read_from_fasta_with_readID(read_id, fasta_path, read_fasta)
            read_blastdb_path = read_fasta.split(".fasta")[0]
            
            # Make BLAST database for this read
            run_makeblastdb(config['External Programs']['makeblastdb_path'],
                          read_fasta,
                          read_blastdb_path)
            
            # Run BLAST-short for restriction enzyme sites
            for enzyme in ['BlnI_seq', 'XapI_seq']:
                run_blastn_short(
                    config['External Programs']['blastn_path'],
                    os.path.join(fasta_dir, f"{enzyme}.fasta"),
                    read_blastdb_path,
                    os.path.join(read_dir, f"{enzyme}.blast.tsv"),
                    word_size=6
                )


            # Extract D4Z4 repeats from read
            D4Z4_blastresult = [result for result in all_results if result.qseqid == "D4Z4"]
            
            # Write D4Z4 blast results
            D4Z4_blast_output = os.path.join(read_dir, f"{read_id}.D4Z4.blastn.tsv")
            write_blastresult_tsv(D4Z4_blastresult, D4Z4_blast_output)
            
            # Extract and separate D4Z4 repeats
            D4Z4_fasta_output = os.path.join(read_dir, f"{read_id}.D4Z4.fasta")
            extract_tr_from_read(read_id, read_fasta, D4Z4_blast_output, D4Z4_fasta_output)
            
            # Add reference to D4Z4 repeats
            os.system(f"cat {config['Reference Sequences Path']['D4Z4_ref']} {D4Z4_fasta_output} > {D4Z4_fasta_output}.with.reference_tr.fasta")
            
            # Create BLAST database for D4Z4 repeats
            D4Z4_blastdb_path = D4Z4_fasta_output.split(".fasta")[0]
            run_makeblastdb(config['External Programs']['makeblastdb_path'],
                          f"{D4Z4_fasta_output}.with.reference_tr.fasta",
                          D4Z4_blastdb_path)
            
            # Run BLAST-short for restriction sites on D4Z4 repeats
            for enzyme in ['BlnI_seq', 'XapI_seq']:
                run_blastn_short(
                    config['External Programs']['blastn_path'],
                    os.path.join(fasta_dir, f"{enzyme}.fasta"),
                    D4Z4_blastdb_path,
                    os.path.join(read_dir, f"{enzyme}.blast.tsv"),
                    word_size=6
                )
            
            # Clean up temporary BLAST database files
            os.system(f"rm {read_blastdb_path}.n*")
            
            # Process restriction sites
                        # BlnI sites
            BlnI_results = browse_blastresult(os.path.join(read_dir, "BlnI_seq.blast.tsv.outfmt6.tsv"))
            if BlnI_results:
                for blnI_result in BlnI_results:
                    # Only process results for individual repeats
                    if 'repeat' not in blnI_result.sseqid:
                        continue
                    direction = 1 if blnI_result.sstart < blnI_result.send else -1
                    repeat_num = int(blnI_result.sseqid.split('repeat')[-1])
                    if direction == 1 and 2750 <= blnI_result.sstart and blnI_result.sstart < 2850:
                        if read_direction == 1:
                            read_d4z4_list[repeat_num-1].BlnI = read_d4z4_list[repeat_num-1].pos[0] + blnI_result.sstart-1
                        elif read_direction == -1:    
                            read_d4z4_list[repeat_num-1].BlnI = read_d4z4_list[repeat_num-1].pos[1] + blnI_result.sstart-1
            
            # XapI sites
            XapI_results = browse_blastresult(os.path.join(read_dir, "XapI_seq.blast.tsv.outfmt6.tsv"))
            if XapI_results:
                for xapI_result in XapI_results:
                    # Only process results for individual repeats
                    if 'repeat' not in xapI_result.sseqid:
                        continue
                    direction = 1 if xapI_result.sstart < xapI_result.send else -1
                    repeat_num = int(xapI_result.sseqid.split('repeat')[-1])
                    if direction == 1 and 1280 <= xapI_result.sstart and xapI_result.sstart <= 1380:
                        if read_direction == 1:
                            read_d4z4_list[repeat_num-1].XapI = read_d4z4_list[repeat_num-1].pos[0] + xapI_result.sstart-1
                        elif read_direction == -1:
                            read_d4z4_list[repeat_num-1].XapI = read_d4z4_list[repeat_num-1].pos[1] + xapI_result.sstart-1
            
            # Create MatchedRead object and generate visualization
            read_pos = (1, read_length) if read_direction == 1 else (read_length, 1)
            matched_read = MatchedRead(
                pos=read_pos,
                d4z4=read_d4z4_list,
                pLAM=pLAM_pos,
                probe=probe_pos,
                pLAM_seq="NA",
                read_id=read_id
            )
            
            # Generate visualization
            draw_arrows(matched_read, os.path.join(read_dir, f"{read_id}.png"))

if __name__ == "__main__":
    main()

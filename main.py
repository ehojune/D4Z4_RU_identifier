import configparser
import os
import argparse
import sys
from operator import attrgetter


class BlastResult:
    def __init__(self,qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen):
        self.qseqid = qseqid
        self.sseqid = sseqid
        self.pident = pident
        self.length = length
        self.mismatch = mismatch
        self.gapopen = gapopen
        self.qstart = qstart
        self.qend = qend
        self.sstart = sstart
        self.send = send
        self.evalue = evalue
        self.bitscore = bitscore
        self.qlen = qlen
        self.slen = slen

# gunzip function is not used; but let users make sure their input fastq to be unzipped in the manual.
def run_gunzip(gzipped_file_path, output_path):
    command = f"gunzip -c {gzipped_file_path} >{output_path}"
    os.system(command)

def convert_fq2fa(sample_fastq_path, output_fasta_path):
        command = f"cat {sample_fastq_path} "
        command += "awk \'\{if(NR%4==1) {printf(\">%s\\n\",substr($0,2));\} else if(NR%4==2) print;}\' > "
        command += output_fasta_path
        os.system(command)
        print(f"fq file {sample_fastq_path} converted to fa file in {output_fasta_path}")

def run_makeblastdb(makeblastdb_path, fasta_path, output_db_path):
    command = makeblastdb_path
    command += " -in fasta_path -input_type nucl"
    command += f" -out {output_db_path}"
    os.system(command)

def run_blastn(makeblastdb_path, query_path, output_path):
    command = blastn_path
    command += f" -query {query_path} -db {makeblastdb_path} "
    command += f" -out {output_path} "
    command += " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\""

    os.system(command)

def browse_blastresult(blast_filename):
    blastresults = []
    with open(blast_filename, 'r') as fr:
        for line in fr:
            if not line.startswith('#'):
                line_splited = line.strip().split('\t')
                qseqid = line_splited[0]
                sseqid = line_splited[1]
                pident = float(line_splited[2])
                length = int(line_splited[3])
                mismatch = int(line_splited[4])
                gapopen = int(line_splited[5])
                qstart = int(line_splited[6])
                qend = int(line_splited[7])
                sstart = int(line_splited[8])
                send = int(line_splited[9])
                evalue = line_splited[10]
                bitscore = int(line_splited[11])
                qlen = int(line_splited[12])
                slen = int(line_splited[13])
                blastresults.append(BlastResult(qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen))
    return blastresults

def filter_blastresult_with_seqid(blastresult_list, seqid):
    blastresult_with_seqid = []
    for blastresult in blastresult_list:
        if blastresult.sseqid == seqid:
            blastresult_with_seqid.append(blastresult)
    return blastresult_with_seqid

def filter_blastresultlist_with_alignment_cov(blastresult_list, coverage):
    # coverage : 0 to 1
    filtered_list = []
    for blastresult in blastresult_list:
        if blastresult.length >= blastresult.qlen * coverage:
            filtered_list.append(blastresult)
    return filtered_list

def filter_blastresultlist_with_pid(blastresult_list, pident):
    # pid : 0 to 100
    filtered_list = []
    for blastresult in blastresult_list:
        if blastresult.pident >= pident:
            filtered_list.append(blastresult)
    return filtered_list

def check_overlap_seqid_in_blastresults(blast_filename1, blast_filename2):
    #only works for outfmt 7
    matched_seqid_results = []
    blastresult1_list = browse_blastresult(blast_filename1)
    blastresult2_list = browse_blastresult(blast_filename2)
    for blastresult1 in blastresult1_list:
        for blastresult2 in blastresult2_list:
            if blastresult1.sseqid == blastresult2.sseqid:
                matched_seqid_results.append(blastresult1.sseqid)
                #print("Found read ID(s) that both files are containing.")
    return list(set(matched_seqid_results))

def write_blastresult_tsv(blastresult_list, save_path):
    fw = open(save_path, 'w')
    for blastresult in blastresult_list:
        fw.write(f"{blastresult.qseqid}\t{blastresult.sseqid}\t\
            {blastresult.pident}\t{blastresult.length}\t{blastresult.mismatch}\t\
            {blastresult.gapopen}\t{blastresult.qstart}\t{blastresult.qend}\t\
            {blastresult.sstart}\t{blastresult.send}\t{blastresult.evalue}\t\
            {blastresult.bitscore}\t{blastresult.qlen}\t{blastresult.slen}\n")
    fw.close()

def sort_blastresult_list(blastresult_list :list, sortby: str, reverse: bool) -> list:
    #order: if False, becomes reverse.
    #need to have """ from operator import attrgetter """
    return sorted(blastresult_list, key = attrgetter(sortby), reverse = reverse)



def print_without_overlap(print_list):
    print_list = list(set(print_list))
    for i in print_list:
        print(i)

"""
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='search a read ID that exists in both blast result file of [query:plam, db:sample_reads] and [query:probe, db:sample_reads]')
    parser.add_argument('--output', required=True, help='full path of output directory')
    args = parser.parse_args()

    output_path = args.output
    _path_blastresult_plam = args.plam

    matched_ids = find_matched_id(_path_blastresult_probe, _path_blastresult_plam)
    print_without_overlap(matched_ids)
"""

def run_porechop():
    pass


if __name__ == "__main__":

    #Read args
    parser = argparse.ArgumentParser(description='search a read ID that exists in both blast result file of [query:plam, db:sample_reads] and [query:probe, db:sample_reads]')
    parser.add_argument('--output', required=True, help='full path of output directory')
    parser.add_argument('--config', required=True, help='path for config file')
    args = parser.parse_args()
    output_path = args.output
    config_path = args.config


    #make output & temp output directory
    tmpdir_path = output_path + '/tmp'
    os.makedirs(tmpdir_path, exist_ok=True)


    #Read config
    config = configparser.ConfigParser()
    config.read(config_path)
    sample_fastq_dict = config['Sample FASTQ Path']
    D4Z4_ref = config['Reference Sequences Path']['D4Z4_ref']
    probe_ref = config['Reference Sequences Path']['probe_ref']
    pLAM_ref = config['Reference Sequences Path']['pLAM_ref']
    BlnI = config['Enzyme Sequences']['BlnI']
    XapI = config['Enzyme Sequences']['XapI']
    pLAM_4qA_seq = config['Plam Feature Sequences']['pLAM_4qA_seq']
    pLAM_10q_seq = config['Plam Feature Sequences']['pLAM_10q_seq']
    blastn_path = config['External Programs']['blastn_path']
    makeblastdb_path = config['External Programs']['makeblastdb_path']
    porechop_path = config['External Programs']['porechop_path']


    #Convert FASTQ to FASTA
    outdir_fasta_path = output_path + '/fasta'
    os.makedirs(outdir_fasta_path, exist_ok=True)
    sample_fasta_dict = {}

    for sample in sample_fastq_dict.keys():
        # get input file path, set output file path
        sample_fastq_path = sample_fastq_dict[sample]
        output_fasta_path = f"{outdir_fasta_path}/{sample}.fa"
        sample_fasta_dict[sample] = output_fasta_path
        # run function
        convert_fq2fa(sample_fastq_path, output_fasta_path)



    # make blastdb for FASTA files
    outdir_blastdb_path = tmpdir_path + '/blastdb'
    os.makedirs(outdir_fasta_path, exist_ok=True)
    sample_blastdb_dict = {}

    for sample in sample_fasta_dict.keys():
        # get input file path, set output file path
        sample_fasta_path = sample_fasta_dict[sample]
        output_sample_blastdb_path = f"{outdir_blastdb_path}/{sample}"
        sample_blastdb_dict[sample] = output_sample_blastdb_path
        os.makedirs(output_sample_blastdb_path, exist_ok=True)
        # run function
        run_makeblastdb(makeblastdb_path, sample_fasta_path, output_sample_blastdb_path)



    # run blastn for pLAM, D4Z4 and probe for samples, find overlapping readID's in blastn results.
    outdir_blastn_path = output_path + '/blast'
    os.makedirs(outdir_blastn_path, exist_ok=True)
    sample_blastn_dict = {}
    matched_readID = {}

    for sample in sample_blastdb_dict.keys():
        # get input file path, set output file path
        sample_blastdb_path = sample_blastdb_dict[sample]
        output_sample_blastn_path = f"{outdir_blastdb_path}/{sample}"
        os.makedirs(output_sample_blastn_path, exist_ok=True)
        output_sample_blastn_pLAM_path = f"{output_sample_blastn_path}/pLAM_{sample}.txt"
        output_sample_blastn_D4Z4_path = f"{output_sample_blastn_path}/D4Z4_{sample}.txt"
        output_sample_blastn_probe_path = f"{output_sample_blastn_path}/probe_{sample}.txt"
        # run function
        run_blastn(sample_blastdb_path, pLAM_ref, output_sample_blastn_pLAM_path)
        run_blastn(sample_blastdb_path, D4Z4_ref, output_sample_blastn_D4Z4_path)
        run_blastn(sample_blastdb_path, probe_ref, output_sample_blastn_probe_path)
        # browse and find overlapping readIDs in blastn results
        pLAM_probe_overlap_seqID = check_overlap_seqid_in_blastresults(output_sample_blastn_pLAM_path, output_sample_blastn_probe_path)
        if not pLAM_probe_overlap_seqID:
            print(f"sample {sample} does not have pLAM & probe overlapping seq.")
            continue
        matched_readID[sample] = pLAM_probe_overlap_seqID
        sample_blastn_dict[sample] = output_sample_blastn_path


    # Terminate pipeline if no sample has pLAM & probe overlapping seq.
    if not matched_readID:
        sys.exit("No sample has pLAM & probe overlapping seq. Terminated.")


    # collect, filter, sort and write blast results of matched seqID
    collected_blastn_result_per_readID = {}

    for sample in sample_blastn_dict.keys():
        # get matched readID for sample
        sample_readID_list = matched_readID[sample]
        # browse sample's blastn file paths
        sample_blastn_path = sample_blastn_dict[sample]
        sample_pLAM_blastresult = browse_blastresult(f"{sample_blastn_path}/pLAM_{sample}.txt")
        sample_D4Z4_blastresult = browse_blastresult(f"{sample_blastn_path}/D4Z4_{sample}.txt")
        sample_probe_blastresult = browse_blastresult(f"{sample_blastn_path}/probe_{sample}.txt")
        sample_blastresult = sample_pLAM_blastresult + sample_D4Z4_blastresult + sample_probe_blastresult
        # filter, sort and write blastn result for each readID
        for readID in sample_readID_list:
            extracted_blastresult_with_readID = filter_blastresult_with_seqid(sample_blastresult, readID)
            filtered_blastresult = filter_blastresultlist_with_alignment_cov(extracted_blastresult_with_readID, 0.9)
            filtered_blastresult = filter_blastresultlist_with_pid(filtered_blastresult, 0.8)
            filtered_blastresult = sort_blastresult_list(filtered_blastresult, 'sstart', False)
            write_blastresult_tsv(f"{sample_blastn_path}/{readID}.blastn.tsv", filtered_blastresult)


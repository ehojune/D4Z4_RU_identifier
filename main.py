import configparser
import os
import argparse

class blastresult:
    def __init__(self, seq_id, pid, length):
        self.seq_id = seq_id
        self.pid = pid
        self.length = length

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
def run_blastn(output_path):
    os.mkdir(output_path+'/blast', exist_ok = True)






    pass









def browse_blastresult(blast_filename):
    blastresults = []
    with open(blast_filename, 'r') as fr:
        for line in fr:
            if not line.startswith('#'):
                line_splited = line.strip().split('\t')
                seq_id = line_splited[1]
                pid = float(line_splited[2])
                length = int(line_splited[3])
                blastresults.append(blastresult(seq_id, pid, length))
    return blastresults

def find_matched_id(blast_filename1, blast_filename2):
    print_list = []
    blastresult1_list = browse_blastresult(blast_filename1)
    blastresult2_list = browse_blastresult(blast_filename2)
    for blastresult1 in blastresult1_list:
        for blastresult2 in blastresult2_list:
            if blastresult1.seq_id == blastresult2.seq_id:
                print_list.append(blastresult1.seq_id)
                #print("Found read ID(s) that both files are containing.")
    return print_list

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
    config.read('config.ini')
    sample_fastq_dict = config['Sample FASTQ Path']
    D4Z4_ref = config['Reference Sequences Path']['D4Z4_ref']
    probe_ref = config['Reference Sequences Path']['probe_ref']
    pLAM_ref = config['Reference Sequences Path']['pLAM_ref']
    BlnI = config['Enzyme Sequences']['BlnI']
    XapI = config['Enzyme Sequences']['XapI']
    pLAM_4qA_seq = config['Plam Feature Sequences']['pLAM_4qA_seq']
    pLAM_10q_seq = config['Plam Feature Sequences']['pLAM_10q_seq']
    blastn_path = config['External Programs']['blastn']
    makeblastdb_path = config['External Programs']['makeblastdb_path']
    porechop_path = config['External Programs']['porechop_path']

    #Convert FASTQ to FASTA
    outdir_fasta_path = output_path + '/fasta'
    os.makedirs(outdir_fasta_path, exist_ok=True)
    sample_fasta_dict = {}
    for sample in sample_fastq_dict.keys():
        sample_fastq_path = sample_fastq_dict[sample]
        output_fasta_path = f"{outdir_fasta_path}/{sample.fa}"
        convert_fq2fa(sample_fastq_path, output_fasta_path)
        sample_fasta_dict[sample] = output_fasta_path

    # make blastdb for fasta files
    outdir_blastdb_path = tmpdir_path + '/blastdb'
    os.makedirs(outdir_fasta_path, exist_ok=True)
    sample_blastdb_dict = {}
    for sample in sample_fastq_dict:
        sample_fasta_path = sample_fastq_dict[sample]
        output_blastdb_path = f"{outdir_fasta_path}/{sample.fa}"
        convert_fq2fa(sample_fastq_path, output_fasta_path)
        sample_blastdb_dict[sample] = output_fasta_path














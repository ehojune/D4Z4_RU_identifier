import configparser
import os
import argparse
import sys
from operator import attrgetter
import matplotlib.pyplot as plt
from multiprocessing import Pool


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

class D4Z4:
	def __init__(self, pos, XapI, BlnI):
		self.pos = pos
		self.XapI = XapI
		self.BlnI = BlnI

class MatchedReads:
	def __init__(self, pos, d4z4, pLAM, probe, pLAM_seq, read_id):
		self.pos = pos
		self.d4z4 = d4z4
		self.pLAM = pLAM
		self.probe = probe
		self.pLAM_seq = pLAM_seq
		self.read_id = read_id

### funcs ###

# gunzip function is not used; but let users make sure their input fastq to be unzipped in the manual.
def run_gunzip(gzipped_file_path, output_path):
	command = f"gunzip -c {gzipped_file_path} >{output_path}"
	os.system(command)

def convert_fq2fa(sample_fastq_path, output_fasta_path):
		command = f"cat {sample_fastq_path} |"
		command += "awk \'{if(NR%4==1) {printf(\">%s\\n\",substr($0,2));} else if(NR%4==2) print;}\' > "
		command += output_fasta_path
		print(command)
		os.system(command)
		print(f"fq file {sample_fastq_path} converted to fa file in {output_fasta_path}")

def run_makeblastdb(makeblastdb_path, fasta_path, output_db_path):
	command = makeblastdb_path
	command += f" -in {fasta_path} -dbtype nucl"
	command += f" -out {output_db_path}"
	os.makedirs(output_db_path, exist_ok=True)
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
			if len(line) == 0: break
			if not line.startswith('#') and not line.startswith('qseqid'):
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
				bitscore = float(line_splited[11])
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
		if blastresult.pident >= pident*100:
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
	fw.write("qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\n")
	for blastresult in blastresult_list:
		fw.write(f"{blastresult.qseqid}\t{blastresult.sseqid}\t{blastresult.pident}\t{blastresult.length}\t{blastresult.mismatch}\t{blastresult.gapopen}\t{blastresult.qstart}\t{blastresult.qend}\t{blastresult.sstart}\t{blastresult.send}\t{blastresult.evalue}\t{blastresult.bitscore}\t{blastresult.qlen}\t{blastresult.slen}\n")
	fw.close()

def sort_blastresult_list(blastresult_list :list, sortby: str, reverse: bool) -> list:
	#order: if False, becomes reverse.
	#need to have """ from operator import attrgetter """
	return sorted(blastresult_list, key = attrgetter(sortby), reverse = reverse)

def extract_read_from_fasta_with_readID(readID, fasta_path, output_path):
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

def extract_tr_positions_from_blastresult_tsv(read_blast_result):
	tr_positions = []
	with open(read_blast_result, 'r') as fr:
		for line in fr:
			if line.startswith("qseqid"): continue #header line
			tr_pos = (int(line.split()[8])-1, int(line.split()[9])-1)
			if len(line)>2: tr_positions.append(tr_pos)
			else: break
	return tr_positions

def extract_seq_from_read_fasta(raw_read_fa):
	with open(raw_read_fa, 'r') as fr:
		fr.readline()
		seq = fr.readline().strip()
	return seq

def revcomp(seq):
	comp_map = {"A": "T", "G": "C", "C": "G", "T": "A"}
	revcomp_seq = ''
	seq = seq[::-1].upper()
	for base in seq:
		revcomp_seq += comp_map[base]
	return revcomp_seq

def extract_tr_from_read(readID, read_fasta_path, tr_read_blast_path, output_path):
	tr_positions = extract_tr_positions_from_blastresult_tsv(tr_read_blast_path)
	seq = extract_seq_from_read_fasta(read_fasta_path)
	tr_list = []

	#If my read is reverse (we can know it by comparing blast's s.start, s.end info)
	if tr_positions[0][0] > tr_positions[0][1]:
		for tr_idx in range(len(tr_positions)):
			tr_pos = tr_positions[tr_idx]
			adj_tr_pos = (tr_pos[1], tr_pos[0]+1)
			tr = revcomp(seq[adj_tr_pos[0]:adj_tr_pos[1]])
			tr_list.append(tr)
		tr_list.reverse()

	#If my read is not-reverse
	else:
		for tr_idx in range(len(tr_positions)):
			tr_pos = tr_positions[tr_idx]
			adj_tr_pos = (tr_pos[0], tr_pos[1]+1)
			tr = seq[adj_tr_pos[0]:adj_tr_pos[1]]
			tr_list.append(tr)

	# write tr units in a file
	with open(output_path, 'w') as fw:
		for tr_idx, tr_seq in enumerate(tr_list):
			fw.write(f">{readID}_repeat{tr_idx+1}\n{tr_seq}\n")

def add_reference_to_tr_extracted_fasta(ref_tr_fa, tr_seperated_fa):
	os.system(f"cat {ref_tr_fa} {tr_seperated_fa} > {tr_seperated_fa}.with.reference_tr.fasta")

def make_query_fasta(query_name, query_target_seq, output_path):
	with open(f"{output_path}/{query_name}.fasta", 'w') as fw:
		fw.write(f">{query_name}\n")
		fw.write(f"{query_target_seq}\n")

def prepare_query_dir(outdir_fasta_path):
	os.makedirs(outdir_fasta_path, exist_ok=True)
	make_query_fasta("BlnI_seq", BlnI_seq, outdir_fasta_path)
	make_query_fasta("XapI_seq", XapI_seq, outdir_fasta_path)
	make_query_fasta("pLAM_4qA_seq", pLAM_4qA_seq, outdir_fasta_path)
	make_query_fasta("pLAM_10q_seq", pLAM_10q_seq, outdir_fasta_path)

def run_blastn_short(makeblastdb_path, query_path, output_path, word_size=False):
	# run with default output format
	command = blastn_path
	command += f" -task blastn-short -dust no -evalue 1000 -soft_masking false -query {query_path} -db {makeblastdb_path} "
	if word_size:
		command += f" -word_size {word_size} "
	command += f" -out {output_path} "
	os.system(command)

	# run again with output format 6
	command = command.split("-out")[0]
	command += f" -out {output_path}.outfmt6.tsv "
	command += " -outfmt \"6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" "
	os.system(command)

def draw_arrows(matched_reads, save_path):
	# 각 숫자에 적절한 간격을 주기 위한 값
	gap = 700

	# 그림 그리기
	fig, ax = plt.subplots(figsize=(10, 2))

	# 방향 설정
	direction = 1 if matched_reads.pos[0] < matched_reads.pos[1] else -1

	# 검은 화살표 그리기 (일반적인 화살표 모양으로 변경)
	ax.annotate("", xy=(matched_reads.pos[1], 0), xytext=(matched_reads.pos[0], 0),
				arrowprops=dict(arrowstyle='->', color='black', lw=2))
	ax.text((matched_reads.pos[0] + matched_reads.pos[1]) / 2, -0.03, matched_reads.read_id, ha='center', va='bottom', color='black', fontsize=8)
	# 검은 화살표 위에 숫자 표시
	ax.text(matched_reads.pos[0], 0.01, str(matched_reads.pos[0]), ha='center', va='bottom', color='black', fontsize=6)
	ax.text(matched_reads.pos[1], 0.01, str(matched_reads.pos[1]), ha='center', va='bottom', color='black', fontsize=6)

	# 빨간 화살표 그리기 (pLAM)
	index_y = 0.03
	ax.annotate("", xy=(matched_reads.pLAM[1], index_y), xytext=(matched_reads.pLAM[0], index_y),
				arrowprops=dict(arrowstyle='->', color='red', lw=2 * direction))
	ax.text((matched_reads.pLAM[0] + matched_reads.pLAM[1]) / 2, index_y -0.01, 'pLAM', ha='center', va='top', color='red', fontsize=8)

	# 빨간 화살표 그리기 (probe)
	ax.annotate("", xy=(matched_reads.probe[1], index_y), xytext=(matched_reads.probe[0], index_y),
				arrowprops=dict(arrowstyle='->', color='red', lw=2 * direction))
	ax.text((matched_reads.probe[0] + matched_reads.probe[1]) / 2, index_y -0.01, 'probe', ha='center', va='top', color='red', fontsize=8)

	# 빨간 화살표 그리기 (d4z4)
	for i, d4z4 in enumerate(matched_reads.d4z4):
		# 특정 위치에 위치된 선 그리기
		if d4z4.BlnI != 'NA':
			ax.plot([d4z4.BlnI, d4z4.BlnI], [index_y + 0.01, index_y + 0.03], color='blue', lw=2 * direction)
			ax.text(d4z4.BlnI, index_y + 0.03, 'BlnI', ha='center', va='bottom', color='blue', fontsize=6)

		if d4z4.XapI != 'NA':
			ax.plot([d4z4.XapI, d4z4.XapI], [index_y + 0.01, index_y + 0.03], color='green', lw=2 * direction)
			ax.text(d4z4.XapI, index_y + 0.03, 'XapI', ha='center', va='bottom', color='green', fontsize=6)

		# 빨간 화살표 그리기 (일반적인 화살표 모양으로 변경)
		ax.annotate("", xy=(d4z4.pos[1], index_y), xytext=(d4z4.pos[0], index_y),
					arrowprops=dict(arrowstyle='->', color='red', lw=2 * direction))
		# 빨간 화살표 위에 시작점과 끝점 표시
		if direction == -1:
			ax.text(d4z4.pos[0], index_y + 0.02, f"{d4z4.pos[0]}", ha='right', va='bottom', color='red', fontsize=6)
			ax.text(d4z4.pos[1], index_y + 0.02, f"{d4z4.pos[1]}", ha='left', va='bottom', color='red', fontsize=6)
		else:
			ax.text(d4z4.pos[0], index_y + 0.02, f"{d4z4.pos[0]}", ha='left', va='bottom', color='red', fontsize=6)
			ax.text(d4z4.pos[1], index_y + 0.02, f"{d4z4.pos[1]}", ha='right', va='bottom', color='red', fontsize=6)
		# 빨간 화살표 아래에 인덱스 표시
		ax.text((d4z4.pos[0] + d4z4.pos[1]) / 2, index_y - 0.01, f"repeat_{i + 1}", ha='center', va='top', color='red', fontsize=8)

	# 그래프 축 설정
	if direction == -1:
		ax.set_xlim(0, max(matched_reads.pos[0], max(d4z4.pos[0] for d4z4 in matched_reads.d4z4)) + 1000)
	else:
		ax.set_xlim(0, max(matched_reads.pos[1], max(d4z4.pos[1] for d4z4 in matched_reads.d4z4)) + 1000)

	ax.set_ylim(-0.1, 0.1)
	ax.axis('off')  # 축 숨기기

	# 그래프 보여주기
	plt.savefig(save_path, format = 'png', dpi = 1200)



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
	samples = list(sample_fastq_dict.keys())
	D4Z4_ref = config['Reference Sequences Path']['D4Z4_ref']
	probe_ref = config['Reference Sequences Path']['probe_ref']
	pLAM_ref = config['Reference Sequences Path']['pLAM_ref']
	BlnI_seq = config['Enzyme Sequences']['BlnI_seq']
	XapI_seq = config['Enzyme Sequences']['XapI_seq']
	pLAM_4qA_seq = config['Plam Feature Sequences']['pLAM_4qA_seq']
	pLAM_10q_seq = config['Plam Feature Sequences']['pLAM_10q_seq']
	blastn_path = config['External Programs']['blastn_path']
	makeblastdb_path = config['External Programs']['makeblastdb_path']
	porechop_path = config['External Programs']['porechop_path']


	# Make BlnI, XapI and pLAM specific sequences to fasta files
	outdir_fasta_path = output_path + '/fasta'
	prepare_query_dir(outdir_fasta_path)





	def process(sample):
		sample_fastq_path = sample_fastq_dict[sample]
		output_fasta_path = f"{outdir_fasta_path}/{sample}.fa"
		convert_fq2fa(sample_fastq_path, output_fasta_path)

	#with Pool(len(samples)) as p:
		#p.map(process, samples)





	def process1(sample):
		run_makeblastdb(makeblastdb_path, f"{outdir_fasta_path}/{sample}.fa",f"{outdir_blastdb_path}/{sample}")


	# make blastdb for FASTA files
	outdir_blastdb_path = tmpdir_path + '/blastdb'
	#os.makedirs(outdir_fasta_path, exist_ok=True)
	#with Pool(len(samples)) as p:
	#	p.map(process1, samples)




	# run blastn for pLAM, D4Z4 and probe for samples, find overlapping readID's in blastn results.
	outdir_blastn_path = output_path + '/blast'
	#os.makedirs(outdir_blastn_path, exist_ok=True)
	matched_readID_dict = {}

	for sample in samples:
		# set input & output files path
		sample_blastdb_path = f"{outdir_blastdb_path}/{sample}"
		output_sample_blastn_path = f"{outdir_blastn_path}/{sample}"
		#os.makedirs(output_sample_blastn_path, exist_ok=True)
		output_sample_blastn_pLAM_path = f"{output_sample_blastn_path}/pLAM_{sample}.txt"
		output_sample_blastn_D4Z4_path = f"{output_sample_blastn_path}/D4Z4_{sample}.txt"
		output_sample_blastn_probe_path = f"{output_sample_blastn_path}/probe_{sample}.txt"
		# run function
		#run_blastn(sample_blastdb_path, pLAM_ref, output_sample_blastn_pLAM_path)
		#run_blastn(sample_blastdb_path, D4Z4_ref, output_sample_blastn_D4Z4_path)
		#run_blastn(sample_blastdb_path, probe_ref, output_sample_blastn_probe_path)
		# browse and find overlapping readIDs in blastn results
		pLAM_probe_overlap_seqID_list = check_overlap_seqid_in_blastresults(output_sample_blastn_pLAM_path, output_sample_blastn_probe_path)
		if not pLAM_probe_overlap_seqID_list:
			print(f"sample {sample} does not have pLAM & probe overlapping seq.")
			continue
		matched_readID_dict[sample] = pLAM_probe_overlap_seqID_list


	# Terminate pipeline if no sample has pLAM & probe overlapping seq.
	if not matched_readID_dict:
		sys.exit("No sample has pLAM & probe overlapping seq. Terminated.")


	# collect, filter, sort and write blast results of matched seqID
	# make an output of blast result of a matched read that contains pLAM, D4Z4 and probe.

	for sample in matched_readID_dict.keys():
		# get matched readID for sample
		sample_readID_list = matched_readID_dict[sample]
		# browse sample's blastn results
		sample_blastn_path = f"{outdir_blastn_path}/{sample}"
		sample_pLAM_blastresult = browse_blastresult(f"{sample_blastn_path}/pLAM_{sample}.txt")
		sample_D4Z4_blastresult = browse_blastresult(f"{sample_blastn_path}/D4Z4_{sample}.txt")
		sample_probe_blastresult = browse_blastresult(f"{sample_blastn_path}/probe_{sample}.txt")
		sample_blastresult = sample_pLAM_blastresult + sample_D4Z4_blastresult + sample_probe_blastresult
		# filter, sort and write blastn result for each readID
		for readID in sample_readID_list:
			os.makedirs(f"{sample_blastn_path}/{readID}", exist_ok=True)
			# make an output of blastn result for pLAM, probe and D4Z4 for each spanning read
			#extracted_blastresult_with_readID = filter_blastresult_with_seqid(sample_blastresult, readID)
			#filtered_blastresult = filter_blastresultlist_with_alignment_cov(extracted_blastresult_with_readID, 0.8)
			#filtered_blastresult = filter_blastresultlist_with_pid(filtered_blastresult, 0.8)
			#sorted_blastresult = sort_blastresult_list(filtered_blastresult, 'sstart', False)
			#write_blastresult_tsv(sorted_blastresult, f"{sample_blastn_path}/{readID}/{readID}.blastn.tsv")
			# make an output of blastn result for only D4Z4 for each spanning read
			#D4Z4_blastresult = [blastresult for blastresult in sorted_blastresult if blastresult.qseqid == "D4Z4"]
			#write_blastresult_tsv(D4Z4_blastresult, f"{sample_blastn_path}/{readID}/{readID}.D4Z4.blastn.tsv")


	# extract matched reads from FASTQ
	print("start extracting matched reads from FASTQ")

	matched_reads_fasta_path = {}
	os.makedirs(f"{output_path}/matched_reads" , exist_ok=True)
	for sample in matched_readID_dict.keys():
		sample_readID_list = matched_readID_dict[sample]
		os.makedirs(f"{output_path}/matched_reads/{sample}", exist_ok=True)
		for readID in sample_readID_list:
			output_read_fasta_dir = f"{output_path}/matched_reads/{sample}/{readID}"
			os.makedirs(output_read_fasta_dir, exist_ok=True)
			output_read_fasta_path = f"{output_read_fasta_dir}/{readID}.fasta"
			sample_fasta_path = f"{outdir_fasta_path}/{sample}.fa"
		#		#	extract_read_from_fasta_with_readID(readID, sample_fasta_path, output_read_fasta_path)
			matched_reads_fasta_path[readID] = output_read_fasta_path
	print("extract D4Z4 from each matched reads FASTQ")

	# extract D4Z4 from each matched reads FASTQ
	for sample in matched_readID_dict.keys():
		sample_readID_list = matched_readID_dict[sample]
		for readID in sample_readID_list:
			read_fasta_path = matched_reads_fasta_path[readID]
			readID_result_path = f"{outdir_blastn_path}/{sample}/{readID}"
			read_D4Z4_blast_path = readID_result_path + f"/{readID}.D4Z4.blastn.tsv"
			#output_D4Z4_seperated_fasta_path = readID_result_path +f"/{readID}.D4Z4.fasta"
			#extract_tr_from_read(readID, read_fasta_path, read_D4Z4_blast_path, output_D4Z4_seperated_fasta_path)
			#add_reference_to_tr_extracted_fasta(D4Z4_ref, output_D4Z4_seperated_fasta_path)

	# make blastdb for each d4z4 of each read of each sample
	# run blastn-short for enzyme, plam
	for sample in matched_readID_dict.keys():
		sample_readID_list = matched_readID_dict[sample]
		for readID in sample_readID_list:
			read_fasta_path = matched_reads_fasta_path[readID]
			readID_result_path = f"{outdir_blastn_path}/{sample}/{readID}"
			output_D4Z4_seperated_fasta_path = readID_result_path +f"/{readID}.D4Z4.fasta"
			D4Z4_blastdb_path = output_D4Z4_seperated_fasta_path.split(".D4Z4")[0]
			read_blastdb_path = read_fasta_path.split(".fasta")[0]

			#run_makeblastdb(makeblastdb_path, output_D4Z4_seperated_fasta_path, D4Z4_blastdb_path)
			#run_makeblastdb(makeblastdb_path, read_fasta_path, read_blastdb_path)

			#run_blastn_short(D4Z4_blastdb_path, f"{outdir_fasta_path}/BlnI_seq.fasta",f"{readID_result_path}/BlnI_seq.blast.tsv", 6)
			#run_blastn_short(D4Z4_blastdb_path, f"{outdir_fasta_path}/XapI_seq.fasta", f"{readID_result_path}/XapI_seq.blast.tsv", 6)
			#run_blastn_short(read_blastdb_path, f"{outdir_fasta_path}/pLAM_4qA_seq.fasta", f"{readID_result_path}/pLAM_4qA_seq.blast.tsv")
			#run_blastn_short(read_blastdb_path, f"{outdir_fasta_path}/pLAM_10q_seq.fasta", f"{readID_result_path}/pLAM_10q_seq.blast.tsv")
			
			os.system(f"rm {D4Z4_blastdb_path}.n*")
			os.system(f"rm {read_blastdb_path}.n*")



	# draw plot
			
	for sample in matched_readID_dict.keys():
		sample_readID_list = matched_readID_dict[sample]
		for readID in sample_readID_list:
			readID_result_path = f"{outdir_blastn_path}/{sample}/{readID}"
			read_blastresults = browse_blastresult(f"{readID_result_path}/{readID}.blastn.tsv")
			read_BlnI_blastresults = browse_blastresult(f"{readID_result_path}/BlnI_seq.blast.tsv.outfmt6.tsv")
			read_XapI_blastresults = browse_blastresult(f"{readID_result_path}/XapI_seq.blast.tsv.outfmt6.tsv")
			
			read_d4z4_list = []
			for blastresult in read_blastresults:
				read_direction = 1 if blastresult.sstart < blastresult.send else -1
				read_length = blastresult.slen
				read_id = blastresult.sseqid
				if blastresult.qseqid == "pLAM":
					pLAM_pos = (blastresult.sstart, blastresult.send)
					continue
				if blastresult.qseqid == "GU550601.1":
					probe_pos = (blastresult.sstart, blastresult.send)
					continue
				else:
					d4z4_pos = D4Z4((blastresult.sstart, blastresult.send), 'NA', 'NA')
					read_d4z4_list.append(d4z4_pos)
			if read_direction == -1:
				read_d4z4_list = read_d4z4_list[::-1]

			if len(read_BlnI_blastresults):
				for blnI_blastresult in read_BlnI_blastresults:
					direction = 1 if blnI_blastresult.sstart < blnI_blastresult.send else -1
					repeat_num = int(blnI_blastresult.sseqid.split('repeat')[-1])
					if direction == 1 and  2750 <= blnI_blastresult.sstart and blnI_blastresult.sstart < 2850:
						if read_direction == 1:
							read_d4z4_list[repeat_num-1].BlnI = read_d4z4_list[repeat_num-1].pos[0] + blnI_blastresult.sstart-1
							continue
						elif read_direction == -1:	
							read_d4z4_list[repeat_num-1].BlnI = read_d4z4_list[repeat_num-1].pos[1] + blnI_blastresult.sstart-1

			if len(read_XapI_blastresults):
				for xapI_blastresult in read_XapI_blastresults:
					direction = 1 if xapI_blastresult.sstart < xapI_blastresult.send else -1
					repeat_num = int(xapI_blastresult.sseqid.split('repeat')[-1])
					if direction == 1 and 1280 <= xapI_blastresult.sstart and xapI_blastresult.sstart <= 1380:
						if read_direction == 1:
							read_d4z4_list[repeat_num-1].XapI = read_d4z4_list[repeat_num-1].pos[0] + xapI_blastresult.sstart-1
							continue
						elif read_direction == -1:
							read_d4z4_list[repeat_num-1].XapI = read_d4z4_list[repeat_num-1].pos[1] + xapI_blastresult.sstart-1

			read_pos = (1, read_length) if read_direction == 1 else (read_length, 1)
			print(read_d4z4_list, pLAM_pos, probe_pos)
			draw_arrows(MatchedReads(read_pos, read_d4z4_list, pLAM_pos, probe_pos, "NA", read_id), f"{readID_result_path}/{readID}.png")







		
			
	
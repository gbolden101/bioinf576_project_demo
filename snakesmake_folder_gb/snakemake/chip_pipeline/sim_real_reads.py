import os
import yaml

def fq_read_cntr(file):
    read_cnt = 0
    with open(file, 'r') as fastq_file:
        for line_num, line in enumerate(fastq_file, start=1):
            if line_num % 4 == 1 and line.startswith('@'):
                read_cnt += 1
    return read_cnt
#source of help https://stackoverflow.com/questions/47204017/opening-and-reading-all-the-files-in-a-directory-in-python-python-beginner
path = '../data/'
def sim_real_fastq(path):
    fileList = os.listdir(path)
    sim_reads = []
    real_reads = []
    for file in fileList:
        if file.endswith(".fastq"):
            file = os.path.join('../data/'+ file)
            if fq_read_cntr(file) > 100:
                #futuresight for snakemake
                real_reads.append(file.replace(".fastq", ""))
            else:
                sim_reads.append(file.replace(".fastq", ""))
    return sim_reads, real_reads
    
real_reads = sim_real_fastq(path)[1]
sim_reads = sim_real_fastq(path)[0]
#snakemake integration and yaml bridge
config_data = {"real_reads": real_reads, "sim_reads": sim_reads}

with open("config_reads.yaml", "w") as config_file:
    yaml.dump(config_data, config_file, default_flow_style=False)

print(f"Filtering complete: {len(real_reads)} real reads, {len(sim_reads)} simulated reads.")

import glob
import subprocess

j = 0
for i in glob.glob("/data/myan/raw_data/pome/pan.raw.fq/11.SNP.smallInDel/02.YS.NoRef.Genome/00.rnaseq/*/*/*.fq.gz"):
    if i.split("/")[-3] not in ["PRJNA802080","PRJNA773212","PRJNA694423"]:
        j += 1
        cmd = ['cp',i,"clean.fastq/"+i.split("/")[-1]]
        print(' '.join(cmd),j)
        
        subprocess.run(cmd)
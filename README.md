# hydractinia_linkage_map
Re-analysis of https://doi.org/10.1186/s12915-023-01532-2 using the chromosome-scale genome

## 1. Download data from NCBI (includes both the F2s and the female)

```
commandsFile = open("SRA_download.sh", 'w')
sampleNames = open("accession_list.txt", "r")
data = sampleNames.readlines()

for rp in data:
        rp = rp.strip()
        sampleName = str(rp) 

        shellScriptName = 'r%s.sh' % (sampleName)
        shellScript = open(shellScriptName, 'w' )

        commandsFile.write('sbatch %s \n' % (shellScriptName))

        shellScript.write("#!/bin/bash\n" )
        shellScript.write("#SBATCH --job-name=%s\n" % (sampleName))
        shellScript.write("#SBATCH --time=0-01:00:00        \n" )
        shellScript.write("#SBATCH --mem=5gb           \n" )
        shellScript.write("#SBATCH --mail-type=NONE          \n" )
        shellScript.write("#SBATCH --ntasks=1                   \n" )
        shellScript.write("#SBATCH --cpus-per-task=5            \n" )
        shellScript.write("#SBATCH --gres=lscratch:100            \n" )
        shellScript.write("#SBATCH --output=k_%s.log\n\n\n" % sampleName)

        shellScript.write("export TMPDIR=/lscratch/$SLURM_JOB_ID \n" )
        shellScript.write("module load sratoolkit \n" )
        shellScript.write("fasterq-dump -p -t /lscratch/$SLURM_JOB_ID -O /data/katzeram/hydractinia_qtl/F2_samples/ "+sampleName+" \n" )
        shellScript.write("echo '<:3)~~~~' \n" )
        

sampleNames.close()
commandsFile.close()
```

```
sh SRA_download.sh
```

Additonally download the male data from zenodo: https://zenodo.org/records/6368105

## 2. BWA-mem2
   index genome:
```
#!/bin/bash
#SBATCH --job-name=bwa_mem2_index
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amanda.katzer@nih.gov
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --output=/data/katzeram/hydractinia_qtl/genome/out_%j.log


cd /data/katzeram/hydractinia_qtl/genome/
module load bwa-mem2
bwa-mem2 index -p hydractiniaT2T GCF_029227915.1_HSymV2.1_genomic.fna
echo '<:3)~~~~'
```

### Map parental samples

Dad
```
#!/bin/bash
#SBATCH --job-name=bwa_mem2_dad
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amanda.katzer@nih.gov
#SBATCH --time=1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50GB
#SBATCH --output=/data/katzeram/hydractinia_qtl/parental_samples/out_%j.log

module load bwa-mem2

bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /data/katzeram/hydractinia_qtl/genome/hydractiniaT2T \
160429_CRATOS_HKM5MBCXX.lane1.12268798.read1.DOWNSAMPLE_0.2.fq.gz \
160429_CRATOS_HKM5MBCXX.lane1.12268798.read2.DOWNSAMPLE_0.2.fq.gz > dad.sam

module load samtools
samtools view -S -b dad.sam > dad.bam
rm dad.sam
echo '<:3)~~~~'
```
Mom
```
#!/bin/bash
#SBATCH --job-name=bwa_mem2_mom
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amanda.katzer@nih.gov
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=50GB
#SBATCH --output=/data/katzeram/hydractinia_qtl/parental_samples/mom_out_%j.log

module load bwa-mem2

bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /data/katzeram/hydractinia_qtl/genome/hydractiniaT2T \
SRR18332201_1.fastq.gz \
SRR18332201_1.fastq.gz > mom.sam

module load samtools
samtools view -S -b mom.sam > mom.bam
rm mom.sam
echo '<:3)~~~~'
```

### Map F2s
```
commandsFile = open("F2_bwa.sh", 'w')
sampleNames = open("accession_list_minus_mom.txt", "r")
data = sampleNames.readlines()

for rp in data:
        rp = rp.strip()
        sampleName = str(rp) 

        shellScriptName = 'r%s.sh' % (sampleName)
        shellScript = open(shellScriptName, 'w' )

        commandsFile.write('sbatch %s \n' % (shellScriptName))

        shellScript.write("#!/bin/bash\n" )
        shellScript.write("#SBATCH --job-name=%s\n" % (sampleName))
        shellScript.write("#SBATCH --time=0-16:00:00        \n" )
        shellScript.write("#SBATCH --mem=50gb           \n" )
        shellScript.write("#SBATCH --mail-type=NONE          \n" )
        shellScript.write("#SBATCH --ntasks=1                   \n" )
        shellScript.write("#SBATCH --cpus-per-task=6            \n" )
        shellScript.write("#SBATCH --output=k_%s.log\n\n\n" % sampleName)

        shellScript.write("cd /data/katzeram/hydractinia_qtl/F2_samples/"+sampleName+" \n" )
        shellScript.write("module load bwa-mem2 \n" )
        shellScript.write("bwa-mem2 mem -t $SLURM_CPUS_PER_TASK /data/katzeram/hydractinia_qtl/genome/hydractiniaT2T "+sampleName+"_1.fastq.gz "+sampleName+"_2.fastq.gz > "+sampleName+".sam \n" )
        shellScript.write("module load samtools \n" )
        shellScript.write("samtools view -S -b "+sampleName+".sam > "+sampleName+".bam \n"  )
        shellScript.write("rm "+sampleName+".sam \n" )
        shellScript.write("echo '<:3)~~~~' \n" )

sampleNames.close()
commandsFile.close()

```

## 3. Mark Duplicate reads with Picard
Make F2 scripts
```
commandsFile = open("F2_bwa.sh", 'w')
sampleNames = open("accession_list_picard.txt", "r")
data = sampleNames.readlines()

for rp in data:
        rp = rp.strip()
        sampleName = str(rp) 

        shellScriptName = 'r%s.sh' % (sampleName)
        shellScript = open(shellScriptName, 'w' )

        commandsFile.write('sbatch %s \n' % (shellScriptName))

        shellScript.write("#!/bin/bash\n" )
        shellScript.write("#SBATCH --job-name=%s\n" % (sampleName))
        shellScript.write("#SBATCH --time=0-06:00:00        \n" )
        shellScript.write("#SBATCH --mem=180gb           \n" )
        shellScript.write("#SBATCH --mail-type=NONE          \n" )
        shellScript.write("#SBATCH --ntasks=1                   \n" )
        shellScript.write("#SBATCH --cpus-per-task=6            \n" )
        shellScript.write("#SBATCH --output=k_%s.log\n\n\n" % sampleName)

        shellScript.write("cd /data/katzeram/hydractinia_qtl/F2_samples/"+sampleName+" \n" )
        shellScript.write("module load samtools \n" )
        shellScript.write("samtools sort -m 60G -@ 3 "+sampleName+".bam -o "+sampleName+".sorted.bam \n" )
        shellScript.write("\n")

        shellScript.write("module load picard \n" )
        shellScript.write("java -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups I="+sampleName+".sorted.bam O="+sampleName+".sorted.rg.bam RGID=1 RGLB=1 RGPL=illumina RGPU=1 RGSM="+sampleName+".sorted.rg.bam \n" )
        shellScript.write("\n")

        shellScript.write("java -Xmx90g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates I="+sampleName+".sorted.rg.bam O="+sampleName+".sorted.rg.md.bam M="+sampleName+"_marked_dup_metrics.txt \n")

        shellScript.write("echo '<:3)~~~~' \n" )

sampleNames.close()
commandsFile.close()
```
Example script
```
#!/bin/bash
#SBATCH --job-name=SRR18332203
#SBATCH --time=0-06:00:00        
#SBATCH --mem=180gb           
#SBATCH --mail-type=NONE          
#SBATCH --ntasks=1                   
#SBATCH --cpus-per-task=6            
#SBATCH --output=k_SRR18332203.log


cd /data/katzeram/hydractinia_qtl/F2_samples/SRR18332203 
module load samtools 
samtools sort -m 60G -@ 3 SRR18332203.bam -o SRR18332203.sorted.bam 

module load picard 
java -jar $PICARDJARPATH/picard.jar AddOrReplaceReadGroups I=SRR18332203.sorted.bam O=SRR18332203.sorted.rg.bam RGID=1 RGLB=1 RGPL=illumina RGPU=1 RGSM=SRR18332203.sorted.rg.bam 

java -Xmx90g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates I=SRR18332203.sorted.rg.bam O=SRR18332203.sorted.rg.md.bam M=SRR18332203_marked_dup_metrics.txt 
echo '<:3)~~~~'
```

Notes:
- Most do not need 180 GB with 6 hours. Only a subset needed that many resources because of samtools. Most used 60 GB at 2 hours. Change the -m under samtools sort to what each THREAD will be using, not the overall amount.
- Picard is very sensitive to having white space after any of the flags. Check second after resource kill error.


## 4. Call variants with GATK

## 5. Filter variant calls

## 6. Genetic map construction

## 7. QTL mapping of sex

## 8. Heterochiasmy 

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

 a. Create dictionary file for genome and index it with Samtools
```
sinteractive --gres=lscratch:50 --cpus-per-task=2 --mem=6g
module load GATK/4.6.0.0

#dictionary file for genome
gatk --java-options "-Xmx5g -Xms5g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" \
CreateSequenceDictionary R=GCF_029227915.1_HSymV2.1_genomic.fna

#index genome with samtools
module load samtools
samtools faidx GCF_029227915.1_HSymV2.1_genomic.fna

#index BAM file
samtools index dad.sorted.rg.md.bam
```

 b. Index each F2
```
#!/bin/bash
#SBATCH --job-name=index_bams
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amanda.katzer@nih.gov
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5GB
#SBATCH --output=/data/katzeram/hydractinia_qtl/F2_samples/F2_gatk_individ_haplotypes/index_bam_out_%j.log

module load samtools

while read sample; do
cd /data/katzeram/hydractinia_qtl/F2_samples/"$sample"
pwd
samtools index "$sample".sorted.rg.md.bam
done <accession_list_minus_mom.txt
```
- Note: if you cannot get a .bai file to create for a sample, you need to rerun the sample from BWA. The large files screw up sometimes. 

 c. Pre-call variants for all the samples
```
commandsFile = open("F2_gatk.sh", 'w')
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
        shellScript.write("#SBATCH --time=0-12:00:00        \n" )
        shellScript.write("#SBATCH --mem=10gb           \n" )
        shellScript.write("#SBATCH --mail-type=NONE          \n" )
        shellScript.write("#SBATCH --ntasks=1                   \n" )
        shellScript.write("#SBATCH --cpus-per-task=4            \n" )
        shellScript.write("#SBATCH --gres=lscratch:50 \n" )
        shellScript.write("#SBATCH --output=k_%s.log\n\n\n" % sampleName)

        shellScript.write("module load GATK \n" )
        shellScript.write('gatk HaplotypeCaller --java-options "-Xmx10g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" -R /data/katzeram/hydractinia_qtl/genome/GCF_029227915.1_HSymV2.1_genomic.fna -I /data/katzeram/hydractinia_qtl/F2_samples/'+sampleName+'/'+sampleName+'.sorted.rg.md.bam --emit-ref-confidence GVCF -O /data/katzeram/hydractinia_qtl/output_gatk_haplotypecaller/'+sampleName+'.out.rawsnps.indels.g.vcf --sequence-dictionary /data/katzeram/hydractinia_qtl/genome/hydractiniagenome.dict \n')
        shellScript.write("echo '<:3)~~~~' \n" )

sampleNames.close()
commandsFile.close()
```

 d. Combine the variants into a single mega GVCF file 
- Highly recommend copying the vcf names into an excel spreadsheet and using concat in excel to add the --variant and the \
```
#!/bin/bash
#SBATCH --job-name=gatk_genotypeGVCF
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amanda.katzer@nih.gov
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --output=/data/katzeram/hydractinia_qtl/output_gatk_haplotypecaller/gatk_genotypegvcf/out_%j.log
#SBATCH --gres=lscratch:50 

module load GATK

gatk --java-options "-Xmx18g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" CombineGVCFs \
 -R /data/katzeram/hydractinia_qtl/genome/GCF_029227915.1_HSymV2.1_genomic.fna \
 --variant dad.out.rawsnps.indels.g.vcf \
 --variant mom.out.rawsnps.indels.g.vcf \
 --variant SRR18332129.out.rawsnps.indels.g.vcf \
 --variant SRR18332130.out.rawsnps.indels.g.vcf \
 --variant SRR18332131.out.rawsnps.indels.g.vcf \
 --variant SRR18332132.out.rawsnps.indels.g.vcf \
 --variant SRR18332133.out.rawsnps.indels.g.vcf \
 --variant SRR18332134.out.rawsnps.indels.g.vcf \
 --variant SRR18332135.out.rawsnps.indels.g.vcf \
 --variant SRR18332136.out.rawsnps.indels.g.vcf \
 --variant SRR18332137.out.rawsnps.indels.g.vcf \
 --variant SRR18332138.out.rawsnps.indels.g.vcf \
 --variant SRR18332139.out.rawsnps.indels.g.vcf \
 --variant SRR18332140.out.rawsnps.indels.g.vcf \
 --variant SRR18332141.out.rawsnps.indels.g.vcf \
 --variant SRR18332142.out.rawsnps.indels.g.vcf \
 --variant SRR18332143.out.rawsnps.indels.g.vcf \
 --variant SRR18332144.out.rawsnps.indels.g.vcf \
 --variant SRR18332145.out.rawsnps.indels.g.vcf \
 --variant SRR18332146.out.rawsnps.indels.g.vcf \
 --variant SRR18332147.out.rawsnps.indels.g.vcf \
 --variant SRR18332148.out.rawsnps.indels.g.vcf \
 --variant SRR18332149.out.rawsnps.indels.g.vcf \
 --variant SRR18332150.out.rawsnps.indels.g.vcf \
 --variant SRR18332151.out.rawsnps.indels.g.vcf \
 --variant SRR18332152.out.rawsnps.indels.g.vcf \
 --variant SRR18332153.out.rawsnps.indels.g.vcf \
 --variant SRR18332154.out.rawsnps.indels.g.vcf \
 --variant SRR18332155.out.rawsnps.indels.g.vcf \
 --variant SRR18332156.out.rawsnps.indels.g.vcf \
 --variant SRR18332157.out.rawsnps.indels.g.vcf \
 --variant SRR18332158.out.rawsnps.indels.g.vcf \
 --variant SRR18332159.out.rawsnps.indels.g.vcf \
 --variant SRR18332160.out.rawsnps.indels.g.vcf \
 --variant SRR18332161.out.rawsnps.indels.g.vcf \
 --variant SRR18332162.out.rawsnps.indels.g.vcf \
 --variant SRR18332163.out.rawsnps.indels.g.vcf \
 --variant SRR18332164.out.rawsnps.indels.g.vcf \
 --variant SRR18332165.out.rawsnps.indels.g.vcf \
 --variant SRR18332166.out.rawsnps.indels.g.vcf \
 --variant SRR18332167.out.rawsnps.indels.g.vcf \
 --variant SRR18332168.out.rawsnps.indels.g.vcf \
 --variant SRR18332169.out.rawsnps.indels.g.vcf \
 --variant SRR18332170.out.rawsnps.indels.g.vcf \
 --variant SRR18332171.out.rawsnps.indels.g.vcf \
 --variant SRR18332172.out.rawsnps.indels.g.vcf \
 --variant SRR18332173.out.rawsnps.indels.g.vcf \
 --variant SRR18332174.out.rawsnps.indels.g.vcf \
 --variant SRR18332175.out.rawsnps.indels.g.vcf \
 --variant SRR18332176.out.rawsnps.indels.g.vcf \
 --variant SRR18332177.out.rawsnps.indels.g.vcf \
 --variant SRR18332178.out.rawsnps.indels.g.vcf \
 --variant SRR18332179.out.rawsnps.indels.g.vcf \
 --variant SRR18332180.out.rawsnps.indels.g.vcf \
 --variant SRR18332181.out.rawsnps.indels.g.vcf \
 --variant SRR18332182.out.rawsnps.indels.g.vcf \
 --variant SRR18332183.out.rawsnps.indels.g.vcf \
 --variant SRR18332184.out.rawsnps.indels.g.vcf \
 --variant SRR18332185.out.rawsnps.indels.g.vcf \
 --variant SRR18332186.out.rawsnps.indels.g.vcf \
 --variant SRR18332187.out.rawsnps.indels.g.vcf \
 --variant SRR18332188.out.rawsnps.indels.g.vcf \
 --variant SRR18332189.out.rawsnps.indels.g.vcf \
 --variant SRR18332190.out.rawsnps.indels.g.vcf \
 --variant SRR18332191.out.rawsnps.indels.g.vcf \
 --variant SRR18332192.out.rawsnps.indels.g.vcf \
 --variant SRR18332193.out.rawsnps.indels.g.vcf \
 --variant SRR18332194.out.rawsnps.indels.g.vcf \
 --variant SRR18332195.out.rawsnps.indels.g.vcf \
 --variant SRR18332196.out.rawsnps.indels.g.vcf \
 --variant SRR18332197.out.rawsnps.indels.g.vcf \
 --variant SRR18332198.out.rawsnps.indels.g.vcf \
 --variant SRR18332199.out.rawsnps.indels.g.vcf \
 --variant SRR18332200.out.rawsnps.indels.g.vcf \
 --variant SRR18332202.out.rawsnps.indels.g.vcf \
 --variant SRR18332203.out.rawsnps.indels.g.vcf \
 --variant SRR18332204.out.rawsnps.indels.g.vcf \
 --variant SRR18332205.out.rawsnps.indels.g.vcf \
 --variant SRR18332206.out.rawsnps.indels.g.vcf \
 --variant SRR18332207.out.rawsnps.indels.g.vcf \
 --variant SRR18332208.out.rawsnps.indels.g.vcf \
 --variant SRR18332209.out.rawsnps.indels.g.vcf \
 --variant SRR18332210.out.rawsnps.indels.g.vcf \
 --variant SRR18332211.out.rawsnps.indels.g.vcf \
 --variant SRR18332212.out.rawsnps.indels.g.vcf \
 --variant SRR18332213.out.rawsnps.indels.g.vcf \
 --variant SRR18332214.out.rawsnps.indels.g.vcf \
 --variant SRR18332215.out.rawsnps.indels.g.vcf \
 --variant SRR18332216.out.rawsnps.indels.g.vcf \
 --variant SRR18332217.out.rawsnps.indels.g.vcf \
 --variant SRR18332218.out.rawsnps.indels.g.vcf \
 --variant SRR18332219.out.rawsnps.indels.g.vcf \
 -O mega.g.vcf.gz
```

 e. Joint genotyping on all samples
```
#!/bin/bash
#SBATCH --job-name=gatk_genotypeGVCF
#SBATCH --mail-type=ALL
#SBATCH --mail-user=amanda.katzer@nih.gov
#SBATCH --time=00:05:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20GB
#SBATCH --output=/data/katzeram/hydractinia_qtl/output_gatk_haplotypecaller/gatk_genotypegvcf/out_%j.log
#SBATCH --gres=lscratch:50 

module load GATK

gatk --java-options "-Xmx18g -Djava.io.tmpdir=/lscratch/$SLURM_JOB_ID" GenotypeGVCFs \
 -R /data/katzeram/hydractinia_qtl/genome/GCF_029227915.1_HSymV2.1_genomic.fna \
 --variant /data/katzeram/hydractinia_qtl/output_gatk_haplotypecaller/mega.g.vcf.gz \
 -O rawvariants.vcf

echo '<:3)~~~~'
```

## 5. Filter variant calls

## 6. Genetic map construction

## 7. QTL mapping of sex

## 8. Heterochiasmy 

# hydractinia_linkage_map
Re-analysis of https://doi.org/10.1186/s12915-023-01532-2 using the chromosome-scale genome

1. Download data from NCBI (includes both the F2s and the female)

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
        

sampleNames.close()
commandsFile.close()
```

```
sh SRA_download.sh
```

Additonally download the male data from zenodo: https://zenodo.org/records/6368105

2. BWA-mem2
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
```

Map parental samples
```

```

Map F2s
```

```

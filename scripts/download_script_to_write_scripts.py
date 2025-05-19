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
	shellScript.write("#SBATCH --time=0-07:00:00        \n" )
	shellScript.write("#SBATCH --mem=5gb           \n" )
	shellScript.write("#SBATCH --mail-type=NONE          \n" )
	shellScript.write("#SBATCH --ntasks=1                   \n" )
	shellScript.write("#SBATCH --cpus-per-task=6            \n" )
	shellScript.write("#SBATCH --gres=lscratch:100            \n" )
	shellScript.write("#SBATCH --output=k_%s.log\n\n\n" % sampleName)

	shellScript.write("export TMPDIR=/lscratch/$SLURM_JOB_ID \n" )
	shellScript.write("module load sratoolkit \n" )
	shellScript.write("mkdir /data/katzeram/hydractinia_qtl/F2_samples/"+sampleName+" \n" )
	shellScript.write("fasterq-dump -p -e 6 -m 5G --split-files --skip-technical -t /lscratch/$SLURM_JOB_ID -O /data/katzeram/hydractinia_qtl/F2_samples/"+sampleName+"/ "+sampleName+" \n" )
	shellScript.write("cd /data/katzeram/hydractinia_qtl/F2_samples/"+sampleName+" \n" )
	shellScript.write("gzip *.fastq \n" )
	

sampleNames.close()
commandsFile.close()

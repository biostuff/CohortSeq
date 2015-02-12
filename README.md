# CohortSeq
## A Toolbox for panel sequencing

##install requirements (besides 'pip'):

We advice the use of virtualenv to allow multiples environments, who already comes with a pip per environment.. 
	pip install virtualenv

After cloning, change into the root of the project directory and run
        virtualenv env
	source env/bin/activate

The dependencies are defined at the requirements.txt file, that you can simply load with
	pip install -r requirements.txt

##Using

python alignment.py --mapper BWA --reference ~/referenceGenomes/hs37d5.fa --sequences "-/bams/someFile.bam" --bam --tempfolder /tmp/ --outputnames /tmp/someFile.bam





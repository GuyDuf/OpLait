Logiciel nécessaire :
Python 3.8.10
snakemake 5.24.1
R 4.1.3
Java 11.0.11
multiqc 1.11

OS : Kubuntu 21.10 x86_64

Les commandes utilisé pour installer les logicels se retrouvent dans le fichier INSTALLATION.txt

Pour lancer le pipeline:

	° Déposer les reads dans le dossier ./data/rawReads
		Les fichiers doivent être nommés selon la structure suivante:
			{id}_L001_{direction}_001.fastq.gz
			ex : IgG1-1_S35_L001_R1_001.fastq.gz
		 		 IgG1-1_S35_L001_R2_001.fastq.gz

	° Dans la racine du dossier executer la commande "snakemake"
	° Les fichier résutants se retrouve dans le dossier ./ouput
	° Les fichiers fasta et fastqc obtenues lors des différentes étapes du snakemake sont disponible dans les sous-dossiers correspondants du dossier ./data

Pour lancer une analyse, les fichiers obtenu par le snakemake ne doivent pas être présents.
Il est possible de nettoyer les dossiers afin de pouvoir lancer le pipeline une seconde fois en executant le script ./source/clean_directory.sh


Pour modifier le pipeline il faut modifier le fichier snakefile
si tu modifie le snakefile assure toi de donner un nom aux input et aux output
Les path ne nécessite pas le ./ au début, il faut éviter de le mettre. Cela pourrais causer des erreurs si on se fie à la documentation de snakemake
	ex. : "data/mergedQC/" au lieu de "./data/mergedQC/"

Les log se retrouve dans le dossier ./.snakemake/log


Pour que le pipeline fonctionne sans problème il faut maintenir la structure de fichier suivante :

.
	/data
		/mergedQC
		/mergedReads
		/rawQC
		/rawReads
		/trimmedQC
		/trimmedReads
	/internal_data
	/logiciel
		/FastQc
		/NGmerge
		/Trimmomatic
		/igblast
	/output
	/out
	/source
	/snakefile



Attention : ne pas supprimer le fichier snakefile

Pour modifier les base de données de référence pour IGHV, IGHD, IGHJ
	Il faut modifier les fichiers ./database/IGHV.fa, ./database/IGHD.fa et ./database/IGHJ.fa
	et supprimer tout les autres fichiers du dossier ./database, les base de données seront recréé lors de l'execution de snakemake

Le dossier ./internal_data permet de faire fonctionner igblast donc il ne faut pas le déplacer

Pour plus d'info sur les logiciels:

multiqc:     https://multiqc.info/
igblast:     https://ncbi.github.io/igblast/
snakemake:   https://snakemake.readthedocs.io/en/stable/
	     https://snakemake.readthedocs.io/en/stable/executing/cli.html
NGmerge:     https://github.com/jsh58/NGmerge
Trimmomatic: http://www.usadellab.org/cms/?page=trimmomatic
			 http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf


Lien github pour le pipeline : https://github.com/GuyDuf/OpLait

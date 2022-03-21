ID,READ = glob_wildcards("data/rawReads/{id}_L001_{read}_001.fastq.gz")

rule all:
    input:
        expand(["data/rawQC/{id}_L001_{read}_001_fastqc.{extension}",
            "data/mergedQC/{id}_fastqc.{extension}",
            "data/trimmedQC/{id}_{paired}_fastqc.{extension}",
            "output/VDJ_{id}.csv",
            "output/multiqc_report.html",
            #"graph/graph{group}.pdf",
            "output/VDJ_{id}.csv_dropped.csv1",
            "output/VDJ_{id}.csv_dropped.csv2",
            #"graph/graph_3.pdf",
            "graph/heatmap{name}.pdf",
            "graph/graphIgBlastDropped.pdf",
            "logiciel/igblast/database/{segment}_clean.fa.{extension2}"],
            id=ID, read=READ, paired=["1P","2P"], extension=["zip","html"], 
            name=["ALL"],
            segment=["IGHV","IGHD","IGHJ"],
            extension2=["ndb","nhr","nin","nog","nos","not","nsq","ntf","nto"])
            #,
            #group=["G1","G2","G3"])

rule rawfastqc:
    input:
        rawread="data/rawReads/{id}_L001_{read}_001.fastq.gz"
    output:
        zip = "data/rawQC/{id}_L001_{read}_001_fastqc.zip",
        html = "data/rawQC/{id}_L001_{read}_001_fastqc.html"
    threads:
        1
    params:
        path="data/rawQC/"
    shell:
        """
        logiciel/FastQC/fastqc {input.rawread} --threads {threads} -o {params.path}
        """

rule trimmomatic:
    input:
        read1="data/rawReads/{id}_L001_R1_001.fastq.gz",
        read2="data/rawReads/{id}_L001_R2_001.fastq.gz"
    output:
        paired_forward="data/trimmedReads/{id}_1P.fastq",
        paired_reverse="data/trimmedReads/{id}_2P.fastq"
    params:
        basename="data/trimmedReads/{id}.fastq",
        log="data/trimmedReads/{id}.log"
    shell:
        """
        java -jar logiciel/Trimmomatic/trimmomatic-0.39.jar PE {input.read1} {input.read2} -baseout {params.basename} LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:100 2> {params.log}
        """
        #	Trimming occurs in the order which the steps are specified on the command line.

rule trimfastqc:
    input:
        trimread1="data/trimmedReads/{id}_1P.fastq",
        trimread2="data/trimmedReads/{id}_2P.fastq"
    output:
        zip="data/trimmedQC/{id}_1P_fastqc.zip",
        html="data/trimmedQC/{id}_1P_fastqc.html",
        zip2="data/trimmedQC/{id}_2P_fastqc.zip",
        html2="data/trimmedQC/{id}_2P_fastqc.html"
    threads:
        1
    params:
        path="data/trimmedQC/"
    shell:
        """
        logiciel/FastQC/fastqc {input.trimread1} {input.trimread2} --threads {threads} -o {params.path}
        """

rule NGmerge:
    input:
        read1="data/trimmedReads/{id}_1P.fastq",
        read2="data/trimmedReads/{id}_2P.fastq"
    output:
        fastq = "data/mergedReads/{id}.fastq",
        log = "data/mergedReads/{id}.log"
    threads:
        1
    shell:
        """
        logiciel/NGmerge/NGmerge -1 {input.read1} -2 {input.read2} -o {output.fastq} -l {output.log} -n {threads}
        """

rule mergedfastqc:
    input:
        mergedread="data/mergedReads/{id}.fastq"
    output:
        zip = "data/mergedQC/{id}_fastqc.zip",
        html = "data/mergedQC/{id}_fastqc.html"
    threads:
        1
    params:
        path="data/mergedQC/"
    shell:
        """
        logiciel/FastQC/fastqc {input.mergedread} --threads {threads} -o {params.path}
        """

rule igblast:
    input:
        mergedread = "output/{id}.fasta",
        IGHV = "logiciel/igblast/database/IGHV_clean.fa",
        IGHD = "logiciel/igblast/database/IGHD_clean.fa",
        IGHJ = "logiciel/igblast/database/IGHJ_clean.fa",
        check = expand(["logiciel/igblast/database/{segment}_clean.fa.{extension}"],
        segment=["IGHV","IGHD","IGHJ"],extension=["ndb","nhr","nin","nog","nos","not","nsq","ntf","nto"])
    output:
        out="output/VDJ_{id}.csv"
    shell:
        """
        logiciel/igblast/bin/igblastn -germline_db_V {input.IGHV} -germline_db_J {input.IGHJ} -germline_db_D {input.IGHD} -organism bovine -query {input.mergedread} -auxiliary_data logiciel/igblast/optional_file/bovine_gl.aux -outfmt 19 > {output.out}
        """

#rule clean:
#    input:
#        IN = expand(["data/mergedReads/{id}.fastq"],id=ID)
#    output:
#        OUT = expand(["output/{id}.fasta",id=ID)
#    shell:
#        """
#        Rscript /mnt/c/Stage/IGH/source/nettoyage.R
#        """

rule fastq2fasta:
    input:
        IN = "data/mergedReads/{id}.fastq"
    output:
        OUT = "output/{id}.fasta"
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input.IN} > {output.OUT}
        """



#rule graphduplicate:
#    input:
#        IN = expand(["data/rawReads/{id}_L001_R1_001.fastq.gz"],id=ID)
	    #IN = expand(["output/VDJ_{id}.csv"],id=ID)
#    output:
#        OUT = "graph/graph_3.pdf"
#    shell:
#        """
#        Rscript ./source/graph3.R {ID}
#        """

rule multiqc:
    input:
       # IN = "output/graph.pdf"
       IN = expand(["output/VDJ_{id}.csv"],id=ID)
    output:
        OUT = "output/multiqc_report.html"
    shell:
        """
        multiqc ./output ./data -o output
        """

rule database_IGHV:
    input:
        IN = "logiciel/igblast/database/IGHV_clean.fa"
    output:
        expand(["logiciel/igblast/database/IGHV_clean.fa.{extension}"], 
            extension=["ndb","nhr","nin","nog","nos","not","nsq","ntf","nto"])
    shell:
        """
        logiciel/igblast/bin/makeblastdb -parse_seqids -dbtype nucl -in {input.IN}
        """

rule database_IGHD:
    input:
        IN = "logiciel/igblast/database/IGHD_clean.fa"
    output:
        expand(["logiciel/igblast/database/IGHD_clean.fa.{extension}"], 
            extension=["ndb","nhr","nin","nog","nos","not","nsq","ntf","nto"])
    shell:
        """"
        logiciel/igblast/bin/makeblastdb -parse_seqids -dbtype nucl -in {input.IN}
        """


rule database_IGHJ:
    input:
        IN = "logiciel/igblast/database/IGHJ_clean.fa"
    output:
        expand(["logiciel/igblast/database/IGHJ_clean.fa.{extension}"], 
            extension=["ndb","nhr","nin","nog","nos","not","nsq","ntf","nto"])
    shell:
        """
        logiciel/igblast/bin/makeblastdb -parse_seqids -dbtype nucl -in {input.IN}
        """

rule name_cleanup:
    input:
        IN = "logiciel/igblast/database/{segment}.fa"
    output:
        OUT = "logiciel/igblast/database/{segment}_clean.fa"
    shell:
        """
        logiciel/igblast/bin/edit_imgt_file.pl {input.IN} > {output.OUT}
        """

rule heatmap:
    input:
        IN = expand(["output/VDJ_{id}.csv_dropped.csv1"], id = ID)
    output:
        OUT = expand(["graph/heatmap{name}.pdf"],
        name=["ALL"])
    shell:
        """
        Rscript ./source/heatmappe.R {ID}
        """

rule IGblastcsvdrop:
    input: 
        IN = expand(["output/VDJ_{id}.csv"], id = ID)
    output:
        OUT1 = expand(["output/VDJ_{id}.csv_dropped.csv1"], id = ID),
        OUT2 = expand(["output/VDJ_{id}.csv_dropped.csv2"], id = ID)
    shell:
        """
        ./source/rmcol.sh output
        ./source/rmcol2.sh output
        """

rule graph:
    input:
        IN  = expand(["data/rawReads/{id}_L001_R1_001.fastq.gz"],id=ID),
        IN2 = expand(["output/VDJ_{id}.csv_dropped.csv1"], id = ID)
	    #IN = expand(["output/VDJ_{id}.csv"],id=ID)
    output:
        OUT = expand(["graph/graph{group}.pdf"], group = ["G1","G2","G3"])
    shell:
        """
        Rscript ./source/graph.R {ID}
        """

rule graphIgBlastDropped:
    input:
        IN = expand(["output/VDJ_{id}.csv_dropped.csv1"], id = ID)
    output:
        OUT = "graph/graphIgBlastDropped.pdf"
    shell:
        """
        Rscript ./source/IgBlastDropped.R {ID}
        """

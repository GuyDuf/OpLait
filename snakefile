ID,READ = glob_wildcards("data/rawReads/{id}_L001_{read}_001.fastq.gz")

rule all:
    input:
        expand(["data/rawQC/{id}_L001_{read}_001_fastqc.{extension}",
            "data/mergedQC/{id}_fastqc.{extension}",
            "data/trimmedQC/{id}_{paired}_fastqc.{extension}",
            "output/VDJ_{id}.csv",
            "output/multiqc_report.html",
            "graph/moss_{group}.pdf",
            "output/VDJ_{id}.csv_dropped.csv1",
            "output/VDJ_{id}.csv_dropped.csv2",
            "graph/heatmap{name}.pdf",
            "graph/coverage{group}.pdf",
            "graph/graphIgBlastDropped.pdf",
            "csv/done.txt",
            "graph/length.pdf",
            "out/distribution{group}_prop.rds",
            "out/distribution{group}_nbr.rds",
            "logiciel/igblast/database/{segment}_clean.fa.{extension2}"],
            id=ID, read=READ, paired=["1P","2P"], extension=["zip","html"],
            name=["ALL"],
            segment=["IGHV","IGHD","IGHJ"],
            extension2=["ndb","nhr","nin","nog","nos","not","nsq","ntf","nto"],
            group=["G1","G2","G3","GM1","GM2","All"])


# Generate a QC for raw reads
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
# Trimm rawReads
# remove nt with a score lesser than 3 at the head and the tail of the read
# Use a sliding window of 4nt and if the average score dips below 15 it cut the read there
# Only keep reads with a minimum length of 100nt
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

# Generate QC from the trimmed reads
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

# Merge the forward and reverse of kept reads
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

# Generate QC of the merged reads
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

# Find the composition (IGHV, IGHJ, IGHD) of each reads
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

# Translate fastq files to fasta
rule fastq2fasta:
    input:
        IN = "data/mergedReads/{id}.fastq"
    output:
        OUT = "output/{id}.fasta"
    shell:
        """
        sed -n '1~4s/^@/>/p;2~4p' {input.IN} > {output.OUT}
        """

# Generate a overall QC
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

# Generate the database of IGHV to use with igBlast
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

# Generate the database of IGHD to use with igBlast
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

# Generate the database of IGHJ to use with igBlast
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
# Clean up the name to use with makeblastdb
rule name_cleanup:
    input:
        IN = "logiciel/igblast/database/{segment}.fa"
    output:
        OUT = "logiciel/igblast/database/{segment}_clean.fa"
    shell:
        """
        logiciel/igblast/bin/edit_imgt_file.pl {input.IN} > {output.OUT}
        """



# CSV cleanup of igBlast result to save up on ram use
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

# Generate CSV of the length of reads to save up on ram use
rule csvLength:
    input:
        IN  = expand(["data/rawReads/{id}_L001_R1_001.fastq.gz"],id=ID),
        IN2 = expand(["output/VDJ_{id}.csv_dropped.csv2"], id = ID)
    output:
        OUT = "csv/done.txt"
    shell:
        """
        Rscript ./source/length_csv.R {ID}
        """

# Generate graph
rule heatmap_graph:
    input:
        IN = expand(["output/VDJ_{id}.csv_dropped.csv1"], id = ID)
    output:
        OUT = expand(["graph/heatmap{name}.pdf"],name=["ALL"])
    shell:
        """
        Rscript ./source/heatmappe.R {ID}
        """

# Generate graph of the length of reads
rule length_graph:
    input:
        IN  = expand(["data/rawReads/{id}_L001_R1_001.fastq.gz"],id=ID),
        IN2 = "csv/done.txt"
    output:
        OUT = expand(["graph/moss_{group}.pdf"], group = ["G1","G2","G3","GM1","GM2","All"]),
        OUT2 = "graph/length.pdf"
    shell:
        """
        Rscript ./source/graph.R {ID}
        """

#generate graph of the dropped reads
rule droppedReads_graph:
    input:
        IN = expand(["output/VDJ_{id}.csv_dropped.csv1"], id = ID)
    output:
        OUT = "graph/graphIgBlastDropped.pdf"
    shell:
        """
        Rscript ./source/IgBlastDropped.R {ID}
        """

#generate graph of the coverage
rule coverage_graph:
    input:
        IN  = expand(["data/rawReads/{id}_L001_R1_001.fastq.gz"],id=ID),
        IN2 = expand(["output/VDJ_{id}.csv_dropped.csv2"], id = ID),
        IN3 = "graph/graphIgBlastDropped.pdf"
    output:
        OUT = expand(["graph/coverage{group}.pdf"], group = ["G1","G2","G3","GM1","GM2","All"])
    shell:
        """
        Rscript ./source/coverage.R {ID}
        """


#generate rdata file IGH composition
rule rdata:
    input:
        IN = expand(["output/VDJ_{id}.csv_dropped.csv1"], id = ID)
    output:
        OUT = expand(["out/distribution{group}_prop.rds"], group = ["G1","G2","G3","GM1","GM2","All"]),
        OUT2 = expand(["out/distribution{group}_nbr.rds"], group = ["G1","G2","G3","GM1","GM2","All"])
    shell:
        """
        Rscript ./source/matrix.R {ID}
        """

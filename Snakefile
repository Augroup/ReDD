scattergather:
    split=8
# rule Unzip:
#     input:
#         "intermediates/fastq/{sample}.*.gz"
#     output:
#         temp("intermediates/fastq/{sample}.fastq")
#     group: 
#         "simple_fastq_operation"
#     resources:
#         runtime = 5
#     threads: 1
#     shell:
#         "zcat {input} > {output}"
rule merge_fastq:
    input:
        "intermediates/fastq/{sample}/"
    output:
        temp("intermediates/fastq/{sample}.fastq")
    group: 
        "simple_fastq_operation"
    threads: 1
    resources:
        runtime= 5
    shell:
        "cat {input}/* > {output}"
# rule fq2fastq:
#     input:
#         "intermediates/fastq/{sample}.fq"
#     output:
#         temp("intermediates/fastq/{sample}.fastq")
#     group: 
#         "simple_fastq_operation"
#     threads: 1
#     resources:
#         runtime= 5
#     shell:
#         "cp {input} {output}"
rule multiline_fastq_to_four_lines:
    input:
        "intermediates/fastq/{sample}.fastq"
    output:
        temp("intermediates/fastq/{sample}_four_lines.fastq")
    group: 
        "simple_fastq_operation"
    threads: 1
    resources:
        runtime=5
    shell:
        "seqtk seq -l0 {input} > {output}"
rule split_fastq:
    input:
        "intermediates/fastq/{sample}_four_lines.fastq"
    output:
        temp(scatter.split("intermediates/splitted/{{sample}}_{scatteritem}.fastq"))
    group: 
        "simple_fastq_operation"
    threads: 1
    resources:
        runtime=5
    run:
        cmd = """fastqsplitter -i {input}"""
        for f in output:
           cmd += f' -o {f}' 
        shell(cmd)
rule NanoFilt:
    input:
        "intermediates/splitted/{sample}_{scatteritem}.fastq"
    output:
        temp("intermediates/splitted/{sample}_{scatteritem}.filtered.fastq")
    group: 
        "simple_fastq_operation"
    threads: 1
    resources:
        runtime= 5
        
    shell:
        "cat {input} | NanoFilt -q 0 --headcrop 5 --tailcrop 3 --readtype 1D > {output}"
rule U2T:
    input:
        "intermediates/splitted/{sample}_{scatteritem}.filtered.fastq"
    output:
        "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq"
    group: 
        "simple_fastq_operation"
    threads: 1
    resources:
        runtime=5
    shell:
        """awk '{{ if (NR%4 == 2) {{gsub(/U/,"T",$1); print $1}} else print }}' {input} > {output}"""
rule minimap2:
    input:
        ref = "intermediates/reference/genome.fa",
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq"
    output:
        temp("intermediates/mapped_reads/{sample}_{scatteritem}.raw.sam")
    resources:
        runtime=30
    group: 
        "simple_sam_operation"
    threads: 4
    shell:
        "minimap2 -G200k --secondary=no -ax splice -uf -k14 -t {threads} {input.ref} {input.fastq} > {output}"
rule filter_sam:
    input:
        "intermediates/mapped_reads/{sample}_{scatteritem}.raw.sam"
    output:
        temp("intermediates/mapped_reads/{sample}_{scatteritem}.sam")
    group: 
        "simple_sam_operation"
    threads: 1
    group: 
        "simple_sam_operation"
    resources:
        runtime=30
    shell:
        """samtools view -H {input} > {output}; cat {input} | awk '{{if(($2=="0")||($2=="16")) print $0}}' >> {output}"""
rule samtools_sort:
    input:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sam"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam"
    threads: 4
    group: 
        "simple_sam_operation"
    resources:
        runtime=30
    shell:
        "samtools view -bS {input} | samtools sort -@ {threads} -o {output}"
rule samtools_index:
    input:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam.bai"
    threads: 4
    group: 
        "simple_sam_operation"
    resources:
        runtime=10
    shell:
        "samtools index -@ {threads} {input}"
rule nanopolish_index:
    input:
        fast5_dir = "intermediates/fast5/{sample}/",
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq"
    output:
       "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq.index"
    threads: 1
    resources:
        runtime=60*4
    shell:
        "nanopolish index -d {input.fast5_dir} {input.fastq}"
rule nanopolish_eventalign:
    input:
        rules.nanopolish_index.output,
        rules.samtools_index.output,
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq",
        bam = "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam",
        ref = "intermediates/reference/genome.fa"
    output:
        "intermediates/eventalign/{sample}_{scatteritem}_eventalign.txt"
    threads: 2
    resources:
        runtime= 24 * 60
    shell:
        "nanopolish eventalign --reads {input.fastq} --bam {input.bam} --genome {input.ref} --samples -t {threads} --print-read-names --signal-index --scale-events > {output}"
rule extract_signals:
    input:
        ref = "intermediates/reference/genome.fa",
        summary = "intermediates/summary/{sample}_summary.txt",
        bam = "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam",
        event = "intermediates/eventalign/{sample}_{scatteritem}_eventalign.txt",
        # candidate=get_candidate_name
    params:
        center="A",
        nt=4,
        featurenum=5,
        buffersize=1000,
        labeltype="I",
        ref_dump_position = "disk"
    output:
        "intermediates/cache/{sample}_{scatteritem}.hdf5",
    threads: 28
    resources:
        runtime=24 * 60
        
    shell:
        """python scripts/extract_bin_data.py --center {params.center} --nt {params.nt} --labeltype {params.labeltype} --maxbuffer_size {params.buffersize} --featurenum {params.featurenum} \
--Ref {input.ref} --summaryIn {input.summary} --samIn {input.bam} --eventIn {input.event} --cachefile {output} --ref_dump_position {params.ref_dump_position} --THREADS {threads}"""
rule predict:
    input:
         "intermediates/cache/{sample}_{scatteritem}.hdf5",
    output:
         "intermediates/{sample}_{scatteritem}.prediction.txt",
    params:
        device='CPU',     #required, select from 'CPU' 'GPU'
        batch_size=10000, #not required, but can be an optional param for user
        model='general' #let the user select from 'general' 'HEK293T_specific' 'GM12878_specific' 'stemcell_specific'
    threads: 5
    resources:
        runtime=5 * 60
    
    shell:
        """python scripts/predict.py --device {params.device} --threads {threads} --batch_size {params.batch_size} --weights_dir  scripts/models/{params.model}.ckpt --test_cachedata {input} --test_outputfile {output} """


rule merge_predicts:
    input:
        gather.split("intermediates/{{sample}}_{scatteritem}.prediction.txt")
    output:
        "outputs/{sample}.prediction.txt"
    threads: 1
    resources:
        runtime= 5
    shell:
        "cat {input} > {output}"
# rule merge_fast5:
#     input:
#         gather.split("intermediates/cache/{{sample}}_{scatteritem}.hdf5")
#     output:
#         "outputs/{sample}_fast5_list.txt"
#     threads: 1
#     resources:
#         runtime=5
        
#     shell:
#         "ls {input} > {output}"
    
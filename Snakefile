scattergather:
    split=8
num_split = workflow._scatter['split']
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
        "intermediates/fastq/{sample}.fastq"
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
        "intermediates/fastq/{sample}_four_lines.fastq"
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
        scatter.split("intermediates/splitted/{{sample}}_{scatteritem}.fastq")
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
        "intermediates/splitted/{sample}_{scatteritem}.filtered.fastq"
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
        "intermediates/mapped_reads/{sample}_{scatteritem}.raw.sam"
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
        "intermediates/mapped_reads/{sample}_{scatteritem}.sam"
    group: 
        "simple_sam_operation"
    threads: 1
    group: 
        "simple_sam_operation"
    resources:
        runtime=30
    shell:
        """samtools view -H {input} > {output}; cat {input} | awk '{{if(($2=="0")||($2=="16")) print $0}}' >> {output}"""
checkpoint samtools_sort:
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
        # nanopolish_index will cost 24 core hourse per flowcell
        runtime= 24 * 60 // num_split // 1
    shell:
        "nanopolish index -d {input.fast5_dir} {input.fastq}"
checkpoint nanopolish_eventalign:
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
        # event_align will cost 150 core hourse per flowcell
        runtime=  150 * 60 // num_split
    shell:
        "nanopolish eventalign --reads {input.fastq} --bam {input.bam} --genome {input.ref} --samples -t {threads} --print-read-names --signal-index --scale-events > {output}"
checkpoint extract_signals:
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
        # extract_signal will cost 672 core hourse per flowcell
        runtime= 672 * 60 // num_split // 28
    shell:
        """python scripts/extract_bin_data.py --center {params.center} --nt {params.nt} --labeltype {params.labeltype} --maxbuffer_size {params.buffersize} --featurenum {params.featurenum} \
--Ref {input.ref} --summaryIn {input.summary} --samIn {input.bam} --eventIn {input.event} --cachefile {output} --ref_dump_position {params.ref_dump_position} --THREADS {threads}"""
checkpoint predict:
    input:
         "intermediates/cache/{sample}_{scatteritem}.hdf5",
    output:
         "intermediates/cache/{sample}_{scatteritem}.prediction.txt",
    params:
        device='CPU',     #required, select from 'CPU' 'GPU'
        batch_size=10000, #not required, but can be an optional param for user
        model='general' #let the user select from 'general' 'HEK293T_specific' 'GM12878_specific' 'stemcell_specific'
    threads: 5
    resources:
        # extract_signal will cost 240 core hourse per flowcell
        runtime= 240 * 60 // num_split // 5
    
    shell:
        """python scripts/predict.py --device {params.device} --threads {threads} --batch_size {params.batch_size} --weights_dir  scripts/models/{params.model}.ckpt --test_cachedata {input} --test_outputfile {output} """


rule merge_predicts:
    input:
        gather.split("intermediates/cache/{{sample}}_{scatteritem}.prediction.txt")
    output:
        "outputs/{sample}.prediction.txt"
    threads: 1
    group: 
        "prediction_aggregation"
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
rule generate_site_level_results:
    input: 
        "outputs/{sample}.prediction.txt"
    output: 
        "outputs/{sample}.site.bed"
    group: 
        "prediction_aggregation"
    params:
        coverage_cutoff=5, #not required, defualt is 5, can be choosed by user.
        ratio_cutoff=0 #not required,defualt is 0, can be choosed by user.
    threads: 6
    resources:
        runtime= 3 * 60 // 6
     
    shell:
        """python scripts/generate_sitelev_info.py --input {input} --output {output} --threads {threads} --cov_cutoff {params.coverage_cutoff} --ratio_cutoff {params.ratio_cutoff} --mod_cov_cutoff 0 """

##genome_sites filtering pipeline
checkpoint genome_site_filtering: #we need parameters to select from genome and cdna
    input:
        site_bed="outputs/{sample}.site.bed",
        alu_bed="intermediates/reference/Hg38_Alu.merge.bed",
        snp_bed="intermediates/reference/hg38_snp151.bed",
        m6A_motif_bed="intermediates/reference/analysis_extract_m6A_motif.bed",
        BEDI_txt="intermediates/reference/REDIportal_hg38.txt"
    output:
        temp_bed_1="outputs/{sample}-1.bed",
        temp_bed_2="outputs/{sample}-2.bed",
        site_tab="outputs/{sample}.flt.tab"
    threads: 1
    group: 
        "prediction_aggregation"
    resources:
        runtime=30
    shell:
        """sh scripts/genome_site_filtering.sh {input.site_bed} {input.alu_bed} {input.snp_bed} {input.m6A_motif_bed} {input.BEDI_txt} {output.temp_bed_1} {output.temp_bed_2} {output.site_tab}"""
rule pre_compute_visualization:
    input:
        site_tab="outputs/{sample}.flt.tab",
        mole_level_bed="outputs/{sample}.prediction.txt",
        ref = "intermediates/reference/genome.fa",
    output:
        directory("outputs/precomputed_visualization/{sample}/")
    threads: 12
    resources:
        runtime= 120
    shell:
        """python scripts/prep_visualization.py --site_level_bed {input.site_tab} --mole_level_bed {input.mole_level_bed} -ref {input.ref} -o {output} -t {threads}"""
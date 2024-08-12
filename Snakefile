scattergather:
    split=8

num_split = workflow._scatter['split']
configfile: "config.yaml"
rule merge_fastq:
    input:
        "intermediates/fastq/{sample}/"
        # config['input_fastq_folder']
    output:
        "intermediates/fastq/{sample}.fastq"
    # group: 
    #     "simple_fastq_operation"
    threads: 1
    resources:
        runtime= 30,
        mem_mb= 100,
    shell:
        "find {input} -name *.fastq -o -name *.fq -type f | xargs cat > {output}"

rule multiline_fastq_to_four_lines:
    input:
        "intermediates/fastq/{sample}.fastq"
    output:
        "intermediates/fastq/{sample}_four_lines.fastq"
    # group: "simple_fastq_operation"
    threads: 1
    resources:
        runtime=30,
        mem_mb= 100,
    shell:
        "seqtk seq -l0 {input} > {output}"
rule split_fastq:
    input:
        "intermediates/fastq/{sample}_four_lines.fastq"
    output:
        scatter.split("intermediates/splitted/{{sample}}_{scatteritem}.fastq")
    # group: 
    #     "simple_fastq_operation"
    threads: 1
    resources:
        runtime=30,
        mem_mb= 100,
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
    # group: 
    #     "simple_fastq_operation"
    threads: 1
    resources:
        runtime= 60,
        mem_mb= 1024,
        
    shell:
        "cat {input} | NanoFilt -q 0 --headcrop 5 --tailcrop 3 --readtype 1D > {output}"
rule U2T:
    input:
        "intermediates/splitted/{sample}_{scatteritem}.filtered.fastq"
    output:
        "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq"
    # group: 
    #     "simple_fastq_operation"
    threads: 1
    resources:
        runtime=30,
        mem_mb= 1024,
    shell:
        """awk '{{ if (NR%4 == 2) {{gsub(/U/,"T",$1); print $1}} else print }}' {input} > {output}"""
if config['reference'] == 'genome':
    rule reference_index:
        input:
            ref_genome = "intermediates/reference/genome.fa",
        output:
            ref_genome_index = "intermediates/reference/genome.fa.fai",
            ref_transcriptome_index = "intermediates/reference/transcriptome.fa.fai",
            ref_transcriptome =  "intermediates/reference/transcriptome.fa",
        resources:
            runtime=60,
            mem_mb= 1024,
        # group: 
        #     "simple_fastq_operation"
        threads: 1
        run:
            shell("""samtools faidx {input.ref_genome}; touch {output.ref_transcriptome_index}; touch {output.ref_transcriptome}""")
elif config['reference'] == 'transcriptome':
    rule reference_index:
        input:
            ref_genome = "intermediates/reference/genome.fa",
            ref_transcriptome =  "intermediates/reference/transcriptome.fa",
        output:
            ref_genome_index = "intermediates/reference/genome.fa.fai",
            ref_transcriptome_index = "intermediates/reference/transcriptome.fa.fai"
        resources:
            runtime=60,
            mem_mb= 1024,
        # group: 
        #     "simple_fastq_operation"
        threads: 1
        run:
            shell("""samtools faidx {input.ref_genome}""")
            shell("""samtools faidx {input.ref_transcriptome}""")
checkpoint minimap2:
    input:
        ref_genome = "intermediates/reference/genome.fa",
        ref_transcriptome =  "intermediates/reference/transcriptome.fa",
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq",
        ref_genome_index = "intermediates/reference/genome.fa.fai"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.raw.sam"
    resources:
        runtime=6 * 60,
        mem_mb=8 * 3 * 1024,
    # group: 
    #     "simple_sam_operation"
    threads: 8
    run:
        if config['reference'] == 'genome':
            shell("minimap2 -G200k --secondary=no -ax splice -uf -k14 -t {threads} {input.ref_genome} {input.fastq} > {output}")
        elif config['reference'] == 'transcriptome':
            shell("minimap2 -G200k --secondary=no -ax map-ont -t {threads} {input.ref_transcriptome} {input.fastq} > {output}")
checkpoint minimap2_genome:
    input:
        ref_genome = "intermediates/reference/genome.fa",
        fastq = "intermediates/fastq/{sample}_four_lines.fastq",
        ref_genome_index = "intermediates/reference/genome.fa.fai"
    output:
        genome_sam = "intermediates/mapped_reads/{sample}.raw.sam",
        genome_filtered_sam = "intermediates/mapped_reads/{sample}.filtered.sam",
        filtered_fastq = "intermediates/fastq/{sample}.filtered.fastq",
        sorted_bam = "igv/{sample}.genome.sorted.bam",
        sorted_bam_bai = "igv/{sample}.genome.sorted.bam.bai"
    resources:
        runtime=6 * 60,
        mem_mb=8 * 3 * 1024,
    # group: 
    #     "simple_sam_operation"
    threads: 8
    shell:
        """
        cat {input.fastq} | NanoFilt -q 0 --headcrop 5 --tailcrop 3 --readtype 1D > {output.filtered_fastq}
        minimap2 -G200k --secondary=no -ax splice -uf -k14 -t {threads} {input.ref_genome} {output.filtered_fastq} > {output.genome_sam}
        samtools view -H {output.genome_sam} > {output.genome_filtered_sam};cat {output.genome_sam} | awk '{{if(($2=="0")||($2=="16")) print $0}}' >> {output.genome_filtered_sam}
        samtools view -bS {output.genome_filtered_sam} | samtools sort -@ {threads} -o {output.sorted_bam}
        samtools index {output.sorted_bam}"""
rule filter_sam:
    input:
        "intermediates/mapped_reads/{sample}_{scatteritem}.raw.sam"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sam"
    # group: 
    #     "simple_sam_operation"
    threads: 1
    resources:
        runtime= 6* 60,
        mem_mb=3 * 1024,
    shell:
        """samtools view -H {input} > {output}; cat {input} | awk '{{if(($2=="0")||($2=="16")) print $0}}' >> {output}"""
rule samtools_sort:
    input:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sam"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam"
    threads: 4
    # group: 
    #     "simple_sam_operation"
    resources:
        runtime= 3 * 60,
        mem_mb=4 * 3 * 1024,
    shell:
        "samtools view -bS {input} | samtools sort -@ {threads} -o {output}"
rule samtools_index:
    input:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam.bai"
    threads: 4
    # group: 
    #     "simple_sam_operation"
    resources:
        runtime= 30,
        mem_mb= 4 * 3 * 1024,
    shell:
        "samtools index -@ {threads} {input}"
rule nanopolish_index:
    input:
        fast5_dir = "intermediates/fast5/{sample}/",
        # fast5_dir = config['input_fast5_folder'],
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq"
    output:
       "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq.index"
    threads: 1
    resources:
        # nanopolish_index will cost 24 core hourse per flowcell
        runtime= 24 * 60 // 3 // 1,
        mem_mb= 9 * 1024
    benchmark: 
        "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq.index.benchmark.tsv"
    
    shell:
        "nanopolish index -d {input.fast5_dir} {input.fastq}"
checkpoint nanopolish_eventalign:
    input:
        rules.nanopolish_index.output,
        rules.samtools_index.output,
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq",
        bam = "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam",
        ref_genome = "intermediates/reference/genome.fa",
        ref_transcriptome =  "intermediates/reference/transcriptome.fa",
    output:
        "intermediates/eventalign/{sample}_{scatteritem}_eventalign.txt"
    threads: 2
    resources:
        # event_align will cost 150 core hours per flowcell
        #runtime=  300 * 60 // 3
        runtime= 24 * 60,
        mem_mb= 9 * 1024
    benchmark: 
        "intermediates/eventalign/{sample}_{scatteritem}_eventalign.benchmark.tsv"
    run:
        if config['reference'] == 'genome':
            shell("nanopolish eventalign --reads {input.fastq} --bam {input.bam} --genome {input.ref_genome} --samples -t {threads} --print-read-names --signal-index --scale-events > {output}")
        elif config['reference'] == 'transcriptome':
            shell("nanopolish eventalign --reads {input.fastq} --bam {input.bam} --genome {input.ref_transcriptome} --samples -t {threads} --print-read-names --signal-index --scale-events > {output}")
checkpoint extract_signals:
    input:
        ref_genome = "intermediates/reference/genome.fa",
        ref_transcriptome =  "intermediates/reference/transcriptome.fa",
        summary = "intermediates/summary/{sample}_summary.txt",
        # summary = config['input_summary_file'],
        bam = "intermediates/mapped_reads/{sample}_{scatteritem}.sorted.bam",
        event = "intermediates/eventalign/{sample}_{scatteritem}_eventalign.txt",
        # candidate=get_candidate_name
    params:
        center=config['center'],
        nt=config['nt'],
        featurenum=config['featurenum'],
        buffersize=config['buffersize'],
        labeltype=config['labeltype'],
        ref_dump_position=config['ref_dump_position']
    output:
        "intermediates/cache/{sample}_{scatteritem}.hdf5",
    #threads: 28
    threads: 6
    resources:
        # extract_signal will cost 672 core hourse per flowcell
        runtime= 2000 * 60 // 3 // 28,
        mem_mb= 6 * 5 * 1024
        # runtime= 60
    
    benchmark: 
        "intermediates/cache/{sample}_{scatteritem}.hdf5.benchmark.tsv"
    
    run:
        if config['reference'] == 'genome':
            shell("""python scripts/extract_bin_data.py --center {params.center} --nt {params.nt} --labeltype {params.labeltype} --maxbuffer_size {params.buffersize} --featurenum {params.featurenum} \
--Ref {input.ref_genome} --summaryIn {input.summary} --samIn {input.bam} --eventIn {input.event} --cachefile {output} --ref_dump_position {params.ref_dump_position} --THREADS {threads}""")
        elif config['reference'] == 'transcriptome':
            shell("""python scripts/extract_bin_data.py --center {params.center} --nt {params.nt} --labeltype {params.labeltype} --maxbuffer_size {params.buffersize} --featurenum {params.featurenum} \
--Ref {input.ref_transcriptome} --summaryIn {input.summary} --samIn {input.bam} --eventIn {input.event} --cachefile {output} --ref_dump_position {params.ref_dump_position} --THREADS {threads}""")
checkpoint predict:
    input:
         "intermediates/cache/{sample}_{scatteritem}.hdf5",
    output:
        "intermediates/cache/{sample}_{scatteritem}.prediction.raw.txt",
    params:
        device = config['device'],#required, select from 'CPU' 'GPU'   
        batch_size=config['batch_size'], #not required, but can be an optional param for user
        model=config['model'] #let the user select from 'general' 'HEK293T_specific' 'GM12878_specific' 'stemcell_specific'
    threads: 5
    resources:
        # extract_signal will cost 240 core hourse per flowcell
        runtime= 480 * 60 // 3 // 5,
         mem_mb= 5 * 5 * 1024
    
    benchmark: 
        "intermediates/cache/{sample}_{scatteritem}.prediction.raw.benchmark.tsv"
    
    shell:
        """python scripts/predict.py --device {params.device} --threads {threads} --batch_size {params.batch_size} --weights_dir  scripts/models/{params.model}.ckpt --test_cachedata {input} --test_outputfile {output} """

#moleoutput = list()
if config['reference'] == 'genome':
    moleoutput="outputs/{sample}.prediction.genome.txt"
elif config['reference'] == 'transcriptome':
    moleoutput="outputs/{sample}.prediction.txt"

checkpoint merge_predicts:
    input:
        gather.split("intermediates/cache/"+"{{sample}}_{scatteritem}.prediction.raw.txt")
    output:
        moleoutput
    threads: 1
    # group: 
    #     "prediction_aggregation"
    resources:
        runtime= 5,
        mem_mb= 1024
    shell:
        "cat {input} > {output}"


mole_input=[]
if config['reference'] == 'genome':
    mole_input.append("outputs/"+"{sample}"+".prediction.genome.txt")
elif config['reference'] == 'transcriptome':
    mole_input.append("outputs/"+"{sample}"+".prediction.txt")

checkpoint generate_site_level_results:
    input: 
        mole_input
    output: 
        "outputs/{sample}.site.bed"
    # group: 
    #     "prediction_aggregation"
    params:
        coverage_cutoff = config['coverage_cutoff'], #not required, defualt is 5, can be choosed by user.
        ratio_cutoff=0 #config['ratio_cutoff'] #not required,defualt is 0.1, can be choosed by user.
    threads: 12
    resources:
        runtime= 12 * 60,
        mem_mb= 12 * 6 * 1024
    shell:
        #"""sh scripts/generate_sitelev_info.sh {input} {output}"""
        """python scripts/generate_sitelev_info.py --input {input} --output {output} --threads {threads} --cov_cutoff {params.coverage_cutoff} --ratio_cutoff {params.ratio_cutoff} --mod_cov_cutoff 0 """
##genome_sites filtering pipeline
if config['reference'] == 'transcriptome':
    rule create_cdna2genome_bed:
        input:
            annotation_gpd="intermediates/reference/gencode.v31.annotation.gpd",
            # annotation_gpd=config['ref_annotation_file']
        output:
            cdna2genome_tab="intermediates/reference/gencode.v31.annotation.cdna2genome.tab",
        resources:
            runtime= 60,
             mem_mb= 6 * 1024
        threads: 1
        shell:
            "perl scripts/convert_gpd_cdna2genome.pl {input.annotation_gpd} {output.cdna2genome_tab}"
else:
    rule create_cdna2genome_bed:
        output:
            cdna2genome_tab="intermediates/reference/gencode.v31.annotation.cdna2genome.tab",
        resources:
            runtime= 1,
            mem_mb=  1024
        threads: 1
        shell:
            "touch {output.cdna2genome_tab}"
if config['filter_m6A']=='True' or config['add_m6A_reference']=='True':
    rule create_m6A_motif:
        input:
            ref_genome = "intermediates/reference/genome.fa"
        output:
            m6A_motif_bed="intermediates/reference/analysis_extract_m6A_motif.bed"
        threads: 1
        # group: 
        #     "prediction_aggregation"
        resources:
            runtime=60,
            mem_mb=  1024
        shell:
            """sh scripts/create_m6A_motif.sh {input.ref_genome} {output.m6A_motif_bed}"""
else:
    rule create_m6A_motif:
        input:
            ref_genome = "intermediates/reference/genome.fa"
        output:
            m6A_motif_bed="intermediates/reference/analysis_extract_m6A_motif.bed"
        threads: 1
        # group: 
        #     "prediction_aggregation"
        resources:
            runtime=1,
            mem_mb=  1024
        shell:
            """touch {output.m6A_motif_bed}"""

if config['reference'] == 'genome':
    site_tab="outputs/"+"{sample}"+".flt.genome.tab",
elif config['reference'] == 'transcriptome':
    site_tab="outputs/"+"{sample}"+".flt.transcriptome.tab"

checkpoint site_filtering: #we need parameters to select from genome and cdna
    input:
        site_bed="outputs/"+"{sample}"+".site.bed",
        alu_bed="intermediates/reference/Hg38_Alu.merge.bed",
        # alu_bed=config['ref_alu_file'],
        snp_bed="intermediates/reference/hg38_snp151.bed",
        # snp_bed=config['ref_snp_file'],
        m6A_motif_bed="intermediates/reference/analysis_extract_m6A_motif.bed",
        REDI_txt="intermediates/reference/REDIportal_hg38.txt",
        # REDI_txt=config['ref_REDIportal_file'],
        cdna2genome_tab="intermediates/reference/gencode.v31.annotation.cdna2genome.tab",
        annotation_gpd="intermediates/reference/gencode.v31.annotation.gpd"
        # annotation_gpd=config['ref_annotation_file']
    output:
        temp_bed_1="outputs/"+"{sample}"+"-1.bed",
        temp_bed_2="outputs/"+"{sample}"+"-2.bed",
        site_tab=site_tab

    params:
        filter_snp=config['filter_snp'],
        filter_m6A=config['filter_m6A'],
        in_alu_coverage=config['in_alu_coverage'],
        out_alu_coverage=config['out_alu_coverage'],
        in_alu_ratio=config['in_alu_ratio'],
        out_alu_ratio=config['out_alu_ratio'],
        coverage_cutoff=config['coverage_cutoff'],
        ratio_cutoff=config['ratio_cutoff']
        
    threads: 6
    # group: 
    #     "prediction_aggregation"
    resources:
        runtime=30,
        mem_mb=  6 * 1024
    run:
        if config['reference'] == 'genome':
            shell("""sh scripts/genome_site_filtering.sh {input.site_bed} {input.alu_bed} {input.snp_bed} {input.m6A_motif_bed} {input.REDI_txt} {output.temp_bed_1} {output.temp_bed_2} {output.site_tab} {params.filter_snp} {params.filter_m6A} {params.in_alu_coverage} {params.out_alu_coverage} {params.in_alu_ratio} {params.out_alu_ratio} {params.coverage_cutoff} {params.ratio_cutoff}""")
        elif config['reference'] == 'transcriptome':
            shell("""sh scripts/transcriptome_site_filtering.sh {input.site_bed} {input.alu_bed} {input.snp_bed} {input.m6A_motif_bed} {input.REDI_txt} {output.temp_bed_1} {output.temp_bed_2} {output.site_tab} {input.cdna2genome_tab} {input.annotation_gpd} {params.filter_snp} {params.filter_m6A}""")

if config['reference'] == 'transcriptome':
    checkpoint cdna2genome:
        input:
            cdna_molecule_input="outputs/"+"{sample}"+".prediction.txt",
            cdna2genome_tab="intermediates/reference/gencode.v31.annotation.cdna2genome.tab",
            cdna_site_tab="outputs/"+"{sample}"+".flt.transcriptome.tab",
        output:
            mole="outputs/"+"{sample}"+".prediction.transcriptome.txt",
            site_tab="outputs/"+"{sample}"+".flt.genome.tab",
        resources:
            runtime= 12 * 60,
            mem_mb=  6 * 10 * 1024
        threads: 6
        shell:
            """sh scripts/cdna2genome.sh {input.cdna_molecule_input} {output.mole} {output.site_tab} {input.cdna2genome_tab} {input.cdna_site_tab}"""

if config['reference'] == 'genome':
    checkpoint pre_compute_visualization:
        input:
            site_tab="outputs/"+"{sample}"+".flt.genome.tab",
            mole_level_bed="outputs/"+"{sample}"+".prediction.genome.txt",
            ref = "intermediates/reference/genome.fa",
            candidate_sites = "intermediates/reference/{sample}.candidate_sites.tab",
            annotation_gpd="intermediates/reference/gencode.v31.annotation.gpd",
            sorted_bam = "igv/{sample}.genome.sorted.bam",
            sorted_bam_bai = "igv/{sample}.genome.sorted.bam.bai"
        params:
            reference_type = "genome"
        output:
            directory("outputs/precomputed_visualization/"+"{sample}"+"/")
        threads: 12
        resources:
            runtime= 12 * 60,
            mem_mb=  12 * 5 * 1024
        shell:
            """python scripts/prep_visualization.py --annotation {input.annotation_gpd} --candidate_sites {input.candidate_sites} --site_level_bed {input.site_tab} --mole_level_bed {input.mole_level_bed} --reference_type {params.reference_type} -ref {input.ref} -o {output} -t {threads}"""
elif config['reference'] == 'transcriptome':
    checkpoint pre_compute_visualization:
        input:
            site_tab="outputs/"+"{sample}"+".flt.genome.tab",
            mole_level_bed="outputs/"+"{sample}"+".prediction.transcriptome.txt",
            ref = "intermediates/reference/genome.fa",
            annotation_gpd="intermediates/reference/gencode.v31.annotation.gpd",
            candidate_sites = "intermediates/reference/{sample}.candidate_sites.tab",
            cdna_site_tab = "outputs/{sample}.flt.transcriptome.tab",
            sorted_bam = "igv/{sample}.genome.sorted.bam",
            sorted_bam_bai = "igv/{sample}.genome.sorted.bam.bai"
        params:
           reference_type = "transcriptome",
           isoform_read_count = 10,
        output:
            sorted_mole_level_bed = "outputs/"+"{sample}"+".prediction.transcriptome.sorted.txt",
            output_dir = directory("outputs/precomputed_visualization/"+"{sample}"+"/")
        threads: 12
        resources:
            runtime= 12 * 60,
            mem_mb=  12 * 5 * 1024
        shell:
            """ 
            sort -S 80% -k 3,3 -k 2,2 {input.mole_level_bed} > {output.sorted_mole_level_bed}
            python scripts/prep_visualization.py --annotation {input.annotation_gpd} --trans_site_level_bed {input.cdna_site_tab} --candidate_sites {input.candidate_sites} --site_level_bed {input.site_tab} --mole_level_bed {output.sorted_mole_level_bed} --reference_type {params.reference_type} -ref {input.ref} -o {output.output_dir} -t {threads}  --isoform_read_count {params.isoform_read_count}
            """

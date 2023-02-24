scattergather:
    split=8

num_split = workflow._scatter['split']
configfile: "config.yaml"

rule merge_fastq:
    input:
        #"intermediates/fastq/{sample}/"
        config['input_fastq_folder']
    output:
        "intermediates/fastq/"+config['sample']+".fastq"
    group: 
        "simple_fastq_operation"
    threads: 1
    resources:
        runtime= 5
    shell:
        "cat {input}/* > {output}"

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
if config['reference'] == 'genome':
    rule reference_index:
        input:
            ref_genome = "intermediates/reference/{sample}.genome.fa",
        output:
            ref_genome_index = "intermediates/reference/{sample}.genome.fa.fai",
            ref_transcriptome_index = "intermediates/reference/{sample}.transcriptome.fa.fai",
            ref_transcriptome =  "intermediates/reference/{sample}.transcriptome.fa",
        resources:
            runtime=30
        group: 
            "simple_fastq_operation"
        threads: 1
        run:
            shell("""samtools faidx {input.ref_genome}; touch {output.ref_transcriptome_index}; touch {output.ref_transcriptome}""")
elif config['reference'] == 'transcriptome':
    rule reference_index:
        input:
            ref_genome = "intermediates/reference/{sample}.genome.fa",
            ref_transcriptome =  "intermediates/reference/{sample}.transcriptome.fa",
        output:
            ref_genome_index = "intermediates/reference/{sample}.genome.fa.fai",
            ref_transcriptome_index = "intermediates/reference/{sample}.transcriptome.fa.fai"
        resources:
            runtime=30
        group: 
            "simple_fastq_operation"
        threads: 1
        run:
            shell("""samtools faidx {input.ref_genome}""")
            shell("""samtools faidx {input.ref_transcriptome}""")
rule minimap2:
    input:
        ref_genome = "intermediates/reference/{sample}.genome.fa",
        ref_transcriptome =  "intermediates/reference/{sample}.transcriptome.fa",
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq",
        ref_genome_index = "intermediates/reference/{sample}.genome.fa.fai"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.raw.sam"
    resources:
        runtime=30,
        # mem_mb=40960,
    # group: 
    #     "simple_sam_operation"
    threads: 8
    run:
        if config['reference'] == 'genome':
            shell("minimap2 -G200k --secondary=no -ax splice -uf -k14 -t {threads} {input.ref_genome} {input.fastq} > {output}")
        elif config['reference'] == 'transcriptome':
            shell("minimap2 -G200k --secondary=no -ax map-ont -t {threads} {input.ref_transcriptome} {input.fastq} > {output}")

rule filter_sam:
    input:
        "intermediates/mapped_reads/{sample}_{scatteritem}.raw.sam"
    output:
        "intermediates/mapped_reads/{sample}_{scatteritem}.sam"
    group: 
        "simple_sam_operation"
    threads: 1
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
        #fast5_dir = "intermediates/fast5/{sample}/",
        fast5_dir = config['input_fast5_folder'],
        fastq = "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq"
    output:
       "intermediates/splitted/{sample}_{scatteritem}.U2T.fastq.index"
    threads: 1
    resources:
        # nanopolish_index will cost 24 core hourse per flowcell
        runtime= 24 * 60 // 3 // 1
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
        ref_genome = "intermediates/reference/{sample}.genome.fa",
        ref_transcriptome =  "intermediates/reference/{sample}.transcriptome.fa",
    output:
        "intermediates/eventalign/{sample}_{scatteritem}_eventalign.txt"
    threads: 2
    resources:
        # event_align will cost 150 core hours per flowcell
        #runtime=  300 * 60 // 3
        runtime= 24 * 60
    benchmark: 
        "intermediates/eventalign/{sample}_{scatteritem}_eventalign.benchmark.tsv"
    run:
        if config['reference'] == 'genome':
            shell("nanopolish eventalign --reads {input.fastq} --bam {input.bam} --genome {input.ref_genome} --samples -t {threads} --print-read-names --signal-index --scale-events > {output}")
        elif config['reference'] == 'transcriptome':
            shell("nanopolish eventalign --reads {input.fastq} --bam {input.bam} --genome {input.ref_transcriptome} --samples -t {threads} --print-read-names --signal-index --scale-events > {output}")
checkpoint extract_signals:
    input:
        ref_genome = "intermediates/reference/{sample}.genome.fa",
        ref_transcriptome =  "intermediates/reference/{sample}.transcriptome.fa",
        #summary = "intermediates/summary/{sample}_summary.txt",
        summary = config['input_summary_file'],
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
        runtime= 2000 * 60 // 3 // 28
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
        runtime= 480 * 60 // 3 // 5
    
    benchmark: 
        "intermediates/cache/{sample}_{scatteritem}.prediction.raw.benchmark.tsv"
    
    shell:
        """python scripts/predict.py --device {params.device} --threads {threads} --batch_size {params.batch_size} --weights_dir  scripts/models/{params.model}.ckpt --test_cachedata {input} --test_outputfile {output} """

#moleoutput = list()
if config['reference'] == 'genome':
    moleoutput="outputs/"+config['sample']+".prediction.genome.txt"
elif config['reference'] == 'transcriptome':
    moleoutput="outputs/"+config['sample']+".prediction.txt"

checkpoint merge_predicts:
    input:
        gather.split("intermediates/cache/"+config['sample']+"_{scatteritem}.prediction.raw.txt")
    output:
        moleoutput
    threads: 1
    # group: 
    #     "prediction_aggregation"
    resources:
        runtime= 5
    shell:
        "cat {input} > {output}"


mole_input=[]
if config['reference'] == 'genome':
    mole_input.append("outputs/"+config['sample']+".prediction.genome.txt")
elif config['reference'] == 'transcriptome':
    mole_input.append("outputs/"+config['sample']+".prediction.txt")

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
        runtime= 12 * 60
    shell:
        # """sh scripts/generate_sitelev_info.sh {input} {output}"""
        """python scripts/generate_sitelev_info.py --input {input} --output {output} --threads {threads} --cov_cutoff {params.coverage_cutoff} --ratio_cutoff {params.ratio_cutoff} --mod_cov_cutoff 0 """
##genome_sites filtering pipeline
if config['reference'] == 'transcriptome':
    rule create_cdna2genome_bed:
        input:
            #annotation_gpd="intermediates/reference/{sample}.gencode.v31.annotation.gpd",
            annotation_gpd=config['ref_annotation_file']
        output:
            cdna2genome_tab="intermediates/reference/"+config['sample']+".gencode.v31.annotation.cdna2genome.tab",
        resources:
            runtime= 60
        threads: 1
        shell:
            "perl scripts/convert_gpd_cdna2genome.pl {input.annotation_gpd} {output.cdna2genome_tab}"
else:
    rule create_cdna2genome_bed:
        output:
            cdna2genome_tab="intermediates/reference/"+config['sample']+".gencode.v31.annotation.cdna2genome.tab",
        resources:
            runtime= 1
        threads: 1
        shell:
            "touch {output.cdna2genome_tab}"
if config['filter_m6A']=='True' or config['add_m6A_reference']=='True':
    rule create_m6A_motif:
        input:
            ref_genome = "intermediates/reference/{sample}.genome.fa"
        output:
            m6A_motif_bed="intermediates/reference/{sample}.analysis_extract_m6A_motif.bed"
        threads: 1
        # group: 
        #     "prediction_aggregation"
        resources:
            runtime=60
        shell:
            """sh scripts/create_m6A_motif.sh {input.ref_genome} {output.m6A_motif_bed}"""
else:
    rule create_m6A_motif:
        input:
            ref_genome = "intermediates/reference/{sample}.genome.fa"
        output:
            m6A_motif_bed="intermediates/reference/{sample}.analysis_extract_m6A_motif.bed"
        threads: 1
        # group: 
        #     "prediction_aggregation"
        resources:
            runtime=1
        shell:
            """touch {output.m6A_motif_bed}"""

if config['reference'] == 'genome':
    site_tab="outputs/"+config['sample']+".flt.genome.tab",
elif config['reference'] == 'transcriptome':
    site_tab="outputs/"+config['sample']+".flt.transcriptome.tab"

rule site_filtering: #we need parameters to select from genome and cdna
    input:
        site_bed="outputs/"+config['sample']+".site.bed",
        #alu_bed="intermediates/reference/"+config['sample']+".Hg38_Alu.merge.bed",
        alu_bed=config['ref_alu_file'],
        #snp_bed="intermediates/reference/"+config['sample']+".hg38_snp151.bed",
        snp_bed=config['ref_snp_file'],
        m6A_motif_bed="intermediates/reference/"+config['sample']+".analysis_extract_m6A_motif.bed",
        #REDI_txt="intermediates/reference/"+config['sample']+".REDIportal_hg38.txt",
        REDI_txt=config['ref_REDIportal_file'],
        cdna2genome_tab="intermediates/reference/"+config['sample']+".gencode.v31.annotation.cdna2genome.tab",
        #annotation_gpd="intermediates/reference/{sample}.gencode.v31.annotation.gpd"
        annotation_gpd=config['ref_annotation_file']
    output:
        temp_bed_1="outputs/"+config['sample']+"-1.bed",
        temp_bed_2="outputs/"+config['sample']+"-2.bed",
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
        runtime=30
    run:
        if config['reference'] == 'genome':
            shell("""sh scripts/genome_site_filtering.sh {input.site_bed} {input.alu_bed} {input.snp_bed} {input.m6A_motif_bed} {input.REDI_txt} {output.temp_bed_1} {output.temp_bed_2} {output.site_tab} {params.filter_snp} {params.filter_m6A} {params.in_alu_coverage} {params.out_alu_coverage} {params.in_alu_ratio} {params.out_alu_ratio} {params.coverage_cutoff} {params.ratio_cutoff}""")
        elif config['reference'] == 'transcriptome':
            shell("""sh scripts/transcriptome_site_filtering.sh {input.site_bed} {input.alu_bed} {input.snp_bed} {input.m6A_motif_bed} {input.REDI_txt} {output.temp_bed_1} {output.temp_bed_2} {output.site_tab} {input.cdna2genome_tab} {input.annotation_gpd} {params.filter_snp} {params.filter_m6A}""")

if config['reference'] == 'transcriptome':
    rule cdna2genome:
        input:
            cdna_molecule_input="outputs/"+config['sample']+".prediction.txt",
            cdna2genome_tab="intermediates/reference/"+config['sample']+".gencode.v31.annotation.cdna2genome.tab",
            cdna_site_tab="outputs/"+config['sample']+".flt.transcriptome.tab",
        output:
            mole="outputs/"+config['sample']+".prediction.transcriptome.txt",
            site_tab="outputs/"+config['sample']+".flt.genome.tab",
        resources:
            runtime= 12 * 60
        threads: 6
        shell:
            """sh scripts/cdna2genome.sh {input.cdna_molecule_input} {output.mole} {output.site_tab} {input.cdna2genome_tab} {input.cdna_site_tab}"""

if config['reference'] == 'genome':
    rule pre_compute_visualization:
        input:
            site_tab="outputs/"+config['sample']+".flt.genome.tab",
            mole_level_bed="outputs/"+config['sample']+".prediction.genome.txt",
            ref = "intermediates/reference/"+config['sample']+".genome.fa",
            candidate_sites = config['ref_candidate_sites_file']
        params:
            reference_type = "genome"
        output:
            directory("outputs/precomputed_visualization/"+config['sample']+"/")
        threads: 12
        resources:
            runtime= 12 * 60
        shell:
            """python scripts/prep_visualization.py --candidate_sites {input.candidate_sites} --site_level_bed {input.site_tab} --mole_level_bed {input.mole_level_bed} --reference_type {params.reference_type} -ref {input.ref} -o {output} -t {threads}"""
elif config['reference'] == 'transcriptome':
    rule pre_compute_visualization:
        input:
            site_tab="outputs/"+config['sample']+".flt.genome.tab",
            mole_level_bed="outputs/"+config['sample']+".prediction.transcriptome.txt",
            ref = "intermediates/reference/"+config['sample']+".genome.fa",
            candidate_sites = config['ref_candidate_sites_file']
        params:
           reference_type = "transcriptome"
           isoform_read_count = 10
        output:
            sorted_mole_level_bed="outputs/"+config['sample']+".prediction.transcriptome.sorted.txt",
            mole_level_isoform_read_counts="outputs/"+config['sample']+".isoform_read_counts.txt",
            directory("outputs/precomputed_visualization/"+config['sample']+"/")
        threads: 12

        resources:
            runtime= 12 * 60
        shell:
            """ 
            sort -S 80% -k 3,3 -k 2,2 {output.mole_level_bed} > {output.sorted_mole_level_bed}
            # awk -v OFS='\t' '!seen[$2, $3]++ {print $2,$3}' {output.sorted_mole_level_bed} > {output.mole_level_isoform_read_counts}
            # python scripts/subsample_reads.py -i {output.mole_level_isoform_read_counts} > {output.selected_reads}
            # awk -v OFS='\t' 'NR==FNR{check[$0];next} $2 in check' {output.selected_reads} {output.mole_level_bed} | sort -S 80% -k 3,3 -k 2,2 > {output.sorted_mole_level_bed}
            python scripts/prep_visualization.py --candidate_sites {input.candidate_sites} --site_level_bed {input.site_tab} --mole_level_bed {input.sorted_mole_level_bed} --reference_type {params.reference_type} -ref {input.ref} -o {output} -t {threads} --isoform_read_list {output.mole_level_isoform_read_counts} --isoform_read_count {params.isoform_read_count}
            """

#"cat 'transID\ttrans_position\tnum_edit\tread_cov\tedit_ratio\teverage_edit_score\tchr\tposition\tgene_id\tstrand\tREDIportal\tSNP\tm6A_motif\tALU' {input.site_tab} > {input.site_tab}" && "cat 'transID\ttrans_position\tnum_edit\tread_cov\tedit_ratio\teverage_edit_score\tchr\tposition\tgene_id\tstrand\tREDIportal\tSNP\tm6A_motif\tALU' {input.mole_level_bed} > {input.mole_level_bed}

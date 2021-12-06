SAMPLES = ['HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
           'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
           'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD']
RADIUS = ['1']
TMPDIR = "/scratch/tereiter"

rule all:
    input: 
        expand("outputs/sourmash_gather/{sample}_gather_gtdb-rs202-genomic-reps.csv", sample = SAMPLES)

rule fastp:
    input:
        R1="inputs/raw/{sample}_R1.fastq.gz",
        R2="inputs/raw/{sample}_R2.fastq.gz",
    output:
        R1="outputs/fastp/{sample}_R1.trim.fq.gz",
        R2="outputs/fastp/{sample}_R2.trim.fq.gz",
        json = "outputs/fastp/{sample}.json",
        html = "outputs/fastp/{sample}.html",
    conda: "envs/fastp.yml"
    threads: 1
    resources:
        mem_mb=16000,
        tmpdir = TMPDIR
    benchmark: "benchmarks/{sample}_fastp.tsv"
    shell:'''
    fastp --in1 {input.R1} \
      --in2 {input.R2} \
      --out1 {output.R1} \
      --out2 {output.R2} \
      --detect_adapter_for_pe \
      --qualified_quality_phred 4 \
      --length_required 31 --correction \
      --json {output.json} \
      --html {output.html}
    '''

rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
    output:
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz',
        human_r1='outputs/bbduk/{sample}_R1.human.fq.gz',
        human_r2='outputs/bbduk/{sample}_R2.human.fq.gz'
    input: 
        r1="outputs/fastp/{sample}_R1.trim.fq.gz",
        r2="outputs/fastp/{sample}_R2.trim.fq.gz",
        human='inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    conda: 'envs/bbmap.yml'
    threads: 1
    resources:
        mem_mb=64000,
        tmpdir = TMPDIR
    benchmark: "benchmarks/{sample}_rmhost.tsv"
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 
        'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        'outputs/bbduk/{sample}_R2.nohost.fq.gz'
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'envs/khmer.yml'
    threads: 1
    resources:
        mem_mb=61000,
        tmpdir = TMPDIR
    benchmark: "benchmarks/{sample}_kmer_trim.tsv"
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule sourmash_sketch:
    input: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: 'outputs/sourmash_sigs/{sample}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}_sketch.tsv"
    shell:"""
    sourmash sketch dna -p k=31,scaled=2000 -o {output} --name {wildcards.sample} {input}
    """

rule gather:
    input:
        sig="outputs/sourmash_sigs/{sample}.sig",
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic-reps.k31.sbt.zip",
    output: 
        csv="outputs/sourmash_gather/{sample}_gather_gtdb-rs202-genomic-reps.csv",
        matches="outputs/sourmash_gather/{sample}_gather_gtdb-rs202-genomic-reps.matches",
        un="outputs/sourmash_gather/{sample}_gather_gtdb-rs202-genomic-reps.un"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 128000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}_gather.tsv"
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

#rule spacegraphcats:
#    input: 
#        #query = "outputs/arg90_matches/cfxA4_AY769933.fna", 
#        conf = "conf/{sample}_r{radius}_conf.yml",
#        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz"
#    output:
#        "outputs/sgc_arg_queries_r{radius}/{sample}_k31_r{radius}_search_oh0/cfxA4_AY769933.fna.cdbg_ids.reads.gz",
#        "outputs/sgc_arg_queries_r{radius}/{sample}_k31_r{radius}_search_oh0/cfxA4_AY769933.fna.contigs.sig"
#    params: outdir = lambda wildcards: "outputs/sgc_genome_queries" + wildcards.radius
#    conda: "envs/spacegraphcats.yml"
#    resources:
#        mem_mb = 64000
#    threads: 1
#    benchmark: "benchmarks/{sample}_r{radius}_spacegraphcats.tsv"
#    shell:'''
#    python -m spacegraphcats run {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
#    '''

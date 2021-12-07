SAMPLES = ['HSM67VF9', 'HSM67VFD', 'HSM67VFJ', 'HSM6XRQB',
           'HSM6XRQI', 'HSM6XRQK', 'HSM6XRQM', 'HSM6XRQO',
           'HSM7CYY7', 'HSM7CYY9', 'HSM7CYYB', 'HSM7CYYD']
TMPDIR = "/scratch/tereiter"

class Checkpoint_GatherResults:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url: 
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_genome_accs(self):
        gather_csv = f'outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv'
        assert os.path.exists(gather_csv)

        genome_accs = []
        with open(gather_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['name'].split(' ')[0]
               genome_accs.append(acc)
        print(f'loaded {len(genome_accs)} accessions from {gather_csv}.')

        return genome_accs

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'gather_gtdb_rep_to_shared_assemblies'; 
        # this will trigger exception until that rule has been run.
        checkpoints.gather_to_sgc_queries.get(**w)

        # parse accessions in gather output file
        genome_accs = self.get_genome_accs()

        p = expand(self.pattern, acc=genome_accs, **w)
        return p

rule all:
    input: 
        Checkpoint_GatherResults("outputs/nbhd_sketch_tables/{acc}_long.csv")

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
    output: 'outputs/sample_sigs/{sample}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}_sketch.tsv"
    shell:"""
    sourmash sketch dna -p k=31,scaled=2000 -o {output} --name {wildcards.sample} {input}
    """

rule gather_reps:
    input:
        sig="outputs/sample_sigs/{sample}.sig",
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic-reps.k31.sbt.zip",
    output: 
        csv="outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic-reps.csv",
        matches="outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic-reps.matches",
        un="outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic-reps.un"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 128000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}_gather.tsv"
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

rule gather_full:
    input:
        sig="outputs/sample_sigs/{sample}.sig",
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic.k31.sbt.zip",
    output: 
        csv="outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic.csv",
        matches="outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic.matches",
        un="outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic.un"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}_gather_full.tsv"
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

# make this a checkpoint to interact with class Checkpoint_GatherResults
checkpoint gather_to_sgc_queries:
    input:  
        gather=expand("outputs/sample_gather/{sample}_gather_gtdb-rs202-genomic.csv", sample = SAMPLES),
    output: 
        gather_grist = "outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv",
        pdf="figures/common_species_breakdown.pdf"
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/gather_gtdb_to_sgc_queries.R"

rule download_sgc_query_genomes:
    input: 
        gather_grist = "outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv",
        conf = "conf/genome-grist-conf.yml"
    output: "genbank_genomes/{acc}_genomic.fna.gz"
    conda: "envs/genome-grist.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell:'''
    genome-grist run {input.conf} --until make_sgc_conf --nolock
    '''
    
rule make_sgc_genome_query_conf_files:
    input:
        csv = "outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv",
        queries = ancient(Checkpoint_GatherResults("genbank_genomes/{acc}_genomic.fna.gz")),
    output:
        conf = "outputs/sgc_conf/{sample}_k31_r1_conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        query_list = "\n- ".join(input.queries)
        with open(output.conf, 'wt') as fp:
           print(f"""\
catlas_base: {wildcards.sample}
input_sequences:
- outputs/abundtrim/{wildcards.sample}.abundtrim.fq.gz
ksize: 31
radius: 1
paired_reads: true
search:
- {query_list}
""", file=fp)

rule spacegraphcats:
    input: 
        queries = ancient(Checkpoint_GatherResults("genbank_genomes/{acc}_genomic.fna.gz")), 
        conf = ancient("outputs/sgc_conf/{sample}_k31_r1_conf.yml"),
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output:
        "outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/results.csv"
        #"outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz",
    params: outdir = "outputs/sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 300000
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
    '''
    
rule touch_spacegraphcats:
    input: 
        "outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/results.csv",
    output: "outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz"
    resources:
        mem_mb = 500
    threads: 1
    shell:'''
    ls {output}
    '''

rule orpheum_translate_reads:
    input: 
        ref="inputs/orpheum_index/{orpheum_db}.{alphabet}-k{ksize}.nodegraph",
        fasta="outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{genome}.fna.cdbg_ids.reads.gz"
    output:
        pep="outputs/orpheum/{sample}_{acc}.coding.faa",
        nuc="outputs/orpheum/{sample}_{acc}.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/{sample}_{acc}.nuc_noncoding.fna",
        csv="outputs/orpheum/{sample}_{acc}.coding_scores.csv",
        json="outputs/orpheum/{sample}_{acc}.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/orpheum-translate-{sample}-{acc}.txt"
    resources:  
        mem_mb=100000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    orpheum translate --alphabet protein --peptide-ksize 10  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fasta} > {output.pep}
    '''

rule sourmash_sketch_species_genomes:
    input: "outputs/orpheum/{sample}_{acc}.coding.faa"
    output: 'outputs/nbhd_sigs/{sample}_{acc}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}_{acc}_sketch.tsv"
    shell:"""
    sourmash sketch protein -p k=10,scaled=100,protein -o {output} --name {wildcards.acc} {input}
    """

rule convert_signature_to_csv:
    input: 'outputs/nbhd_sigs/{sample}_{acc}.sig'
    output: 'outputs/nbhd_sigs/{sample}_{acc}.csv'
    conda: 'envs/sourmash.yml'
    threads: 1
    benchmark: "benchmarks/{sample}_{acc}_sketch_to_csv.tsv"
    resources:
        mem_mb=2000,
        tmpdir = TMPDIR
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_table_long:
    input: 
        expand("outputs/nbhd_sigs/{sample}_{{acc}}.csv", sample = SAMPLES)
    output: csv = "outputs/nbhd_sketch_tables/{acc}_long.csv"
    conda: 'envs/r.yml'
    threads: 1
    benchmark: "benchmarks/sourmash_sketch_tables_long/{acc}.txt"
    resources:
        mem_mb=64000,
        tmpdir = TMPDIR
    script: "scripts/sketch_csv_to_long.R"
    
# rule pagoo


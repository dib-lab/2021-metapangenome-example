import csv
import sys
import urllib.request
import pdb

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


class Checkpoint_AccToDbs:
    """
    Define a class a la genome-grist to simplify file specification
    from checkpoint (e.g. solve for {acc} wildcard). This approach
    is documented at this url:
    http://ivory.idyll.org/blog/2021-snakemake-checkpoints.html
    """
    def __init__(self, pattern):
        self.pattern = pattern

    def get_acc_dbs(self):
        acc_db_csv = f'outputs/genbank/gather_gtdb-rs202-genomic.x.dbs.csv'
        assert os.path.exists(acc_db_csv)

        acc_dbs = []
        with open(acc_db_csv, 'rt') as fp:
           r = csv.DictReader(fp)
           for row in r:
               acc = row['accession']
               db = row['species']
               acc_db = acc + "-" + db
               acc_dbs.append(acc_db)

        return acc_dbs

    def __call__(self, w):
        global checkpoints

        # wait for the results of 'query_to_species_db';
        # this will trigger exception until that rule has been run.
        checkpoints.query_to_species_db.get(**w)

        # parse accessions in gather output file
        genome_acc_dbs = self.get_acc_dbs()

        p = expand(self.pattern, acc_db=genome_acc_dbs, **w)
        return p


rule all:
    input: 
        #expand("outputs/metabat2_gather/{sample}_done.txt", sample = SAMPLES),
        #Checkpoint_AccToDbs("outputs/metabat2_gather_labelled_bins/{acc_db}.txt"),
        Checkpoint_GatherResults("outputs/nbhd_sketch_tables/{acc}_long.csv"),
        Checkpoint_AccToDbs("outputs/pagoo_species/{acc_db}_binmap.pdf"),
        #Checkpoint_GatherResults(expand("outputs/nbhd_gather/{sample}-{{acc}}_gather_gtdb-rs202-genomic.csv", sample = SAMPLES)),
        expand("outputs/orpheum_compare/{sample}_comp_containment.csv", sample = SAMPLES),
        expand("outputs/orpheum_compare/{sample}_comp.csv", sample = SAMPLES),
        Checkpoint_AccToDbs(expand("outputs/nbhd_species_gather/{sample}-{{acc_db}}_gather_gtdb-rs202-genomic.csv", sample = SAMPLES)),
        Checkpoint_AccToDbs("outputs/compare_all/{acc_db}_containment.csv"),
        Checkpoint_AccToDbs("outputs/compare_core/{acc_db}_containment.csv"),
        # csvs for upset plots
        Checkpoint_AccToDbs("outputs/metabat2_prokka_sigs_all/{acc_db}_all_kmers.csv"),
        Checkpoint_AccToDbs("outputs/roary_sigs_all/{acc_db}_pan_genome_reference_all_genes.csv"),
        Checkpoint_AccToDbs('outputs/nbhd_sigs_species_all/{acc_db}_all_kmers.csv')

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

##############################################################
## de novo assemble and bin to compare against metapangenome
##############################################################


rule split_paired_reads:
    input: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: 
        r1= "outputs/abundtrim_split/{sample}_R1.fq.gz",
        r2= "outputs/abundtrim_split/{sample}_R2.fq.gz",
        singleton = "outputs/abundtrim_split/{sample}_orphan.fq.gz"
    threads: 1
    resources:
        mem_mb=8000,
        tmpdir = TMPDIR
    conda: 'envs/bbmap.yml'
    shell:'''
    repair.sh in1={input} out1={output.r1} out2={output.r2} outs={output.singleton} repair
    '''

rule assemble:
    input: 
        r1 = "outputs/abundtrim_split/{sample}_R1.fq.gz",
        r2 = "outputs/abundtrim_split/{sample}_R2.fq.gz"
    output:
        tmp_dir = directory("outputs/megahit/{sample}_tmp"), 
        fa="outputs/megahit/{sample}.contigs.fa"
    conda: 'envs/megahit.yml'
    threads: 1
    resources:
        mem_mb=32000,
        tmpdir = TMPDIR
    shell:'''
    megahit -1 {input.r1} -2 {input.r2} --out-dir {output.tmp_dir} --out-prefix {wildcards.sample}
    mv {output.tmp_dir}/{wildcards.sample}.contigs.fa {output.fa}
    '''

rule gen_coverages_bowtie_build:
    input: "outputs/megahit/{sample}.contigs.fa"
    output: "outputs/megahit/{sample}.contigs.fa.1.bt2"
    conda: 'envs/bowtie2.yml'
    threads: 1
    resources:
        mem_mb=4000,
        tmpdir = TMPDIR
    shell:'''
    bowtie2-build {input} {input}    
    '''

rule gen_coverages_bowtie:
    input:
        bwt="outputs/megahit/{sample}.contigs.fa.1.bt2",
        fa="outputs/megahit/{sample}.contigs.fa",
        r1 = "outputs/abundtrim_split/{sample}_R1.fq.gz",
        r2 = "outputs/abundtrim_split/{sample}_R2.fq.gz"
    output: "outputs/bowtie2/{sample}_to_sort.bam"
    conda: 'envs/bowtie2.yml'
    threads: 1
    resources:
        mem_mb=6000,
        tmpdir = TMPDIR
    shell:'''
    bowtie2 -x {input.fa} -1 {input.r1} -2 {input.r2} |
        samtools view -bS -o {output}
    '''

rule gen_coverages_sort:
    input: "outputs/bowtie2/{sample}_to_sort.bam"
    output: "outputs/bowtie2/{sample}.bam"
    conda: 'envs/bowtie2.yml'
    threads: 1
    resources:
        mem_mb=6000,
        tmpdir = TMPDIR
    shell:'''
    samtools sort {input} -o {output}
    '''

rule gen_coverages_index:
    input: "outputs/bowtie2/{sample}.bam"
    output: "outputs/bowtie2/{sample}.bam.bai"
    conda: 'envs/bowtie2.yml'
    threads: 1
    resources:
        mem_mb=6000,
        tmpdir = TMPDIR
    shell:'''
    samtools index {input}
    '''

rule get_metabat_depth_file:
    input:
        bam=  "outputs/bowtie2/{sample}.bam",
        indx= "outputs/bowtie2/{sample}.bam.bai"
    output: "outputs/metabat2/{sample}_metabat_depth.txt"
    conda: 'envs/metabat2.yml'
    threads: 1
    resources:
        mem_mb=6000,
        tmpdir = TMPDIR
    shell:'''
    jgi_summarize_bam_contig_depths --outputDepth {output} {input.bam}
    '''

checkpoint metabat:
    input:
        depth="outputs/metabat2/{sample}_metabat_depth.txt",
        fa="outputs/megahit/{sample}.contigs.fa",
    output: directory("outputs/metabat2/{sample}")
    params: outdir = lambda wildcards: "outputs/metabat2/" + wildcards.sample + "/" + "bin" 
    conda: 'envs/metabat2.yml'
    threads: 1
    resources:
        mem_mb=6000,
        tmpdir = TMPDIR
    shell:'''
    metabat2 -m 1500 -i {input.fa} --abdFile {input.depth} -o {params.outdir}
    #touch {output}
    '''

rule metabat_sketch:
    input: "outputs/metabat2/{sample}/bin.{bin}.fa"
    output: 'outputs/metabat2_sigs/{sample}_bin.{bin}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sketch dna -p k=31,scaled=2000 -o {output} --name {wildcards.sample}.bin{wildcards.bin} {input}
    """

rule metabat_gather:
    input:
        sig="outputs/metabat2_sigs/{sample}_bin.{bin}.sig",
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic.k31.sbt.zip",
    output: 
        csv="outputs/metabat2_gather/{sample}_bin.{bin}_gather_gtdb-rs202-genomic.csv",
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

def checkpoint_metabat2_1(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.metabat.get(**wildcards).output[0]    
    file_names = expand("outputs/metabat2_gather/{{sample}}_bin.{bin}_gather_gtdb-rs202-genomic.csv",
                        bin = glob_wildcards(os.path.join(checkpoint_output, "bin.{bin}.fa")).bin)
    return file_names

rule metabat_gather_dummy:
    input: checkpoint_metabat2_1
    output: touch("outputs/metabat2_gather/{sample}_done.txt")
    
     
rule metabat_prokka:
    input: "outputs/metabat2/{sample}/bin.{bin}.fa"
    output: "outputs/metabat2_prokka/{sample}_bin.{bin}.faa"
    conda: 'envs/prokka.yml'
    resources:
        mem_mb = 8000,
        tmpdir = TMPDIR
    threads: 2
    params: 
        outdir = lambda wildcards: 'outputs/metabat2_prokka/'
    shell:'''
    prokka {input} --outdir {params.outdir} --prefix {wildcards.sample}_bin.{wildcards.bin} \
        --metagenome --force --locustag {wildcards.sample}_bin.{wildcards.bin} --cpus {threads} --centre X --compliant
    '''

rule metabat_prokka_sketch:
    input: "outputs/metabat2_prokka/{sample}_bin.{bin}.faa"
    output: 'outputs/metabat2_prokka_sigs/{sample}_bin.{bin}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sketch protein -p k=10,scaled=100,protein -o {output} --name {wildcards.sample}_bin.{wildcards.bin} {input}
    """

def checkpoint_metabat2_2(wildcards):
    # checkpoint_output encodes the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.metabat.get(**wildcards).output[0]    
    file_names = expand("outputs/metabat2_prokka_sigs/{{sample}}_bin.{bin}.sig",
                        bin = glob_wildcards(os.path.join(checkpoint_output, "bin.{bin}.fa")).bin)
    return file_names

rule metabat_prokka_sketch_dummy:
    input: checkpoint_metabat2_2
    output: touch("outputs/metabat2_prokka_sigs/{sample}_done.txt")

rule metabat_generate_lists_of_bins_of_same_species:
    input: 
        lineages = "inputs/gtdb-rs202.taxonomy.v2.csv",
        acc_to_db = "outputs/genbank/gather_gtdb-rs202-genomic.x.dbs.csv",
        dummy1 = expand("outputs/metabat2_gather/{sample}_done.txt", sample = SAMPLES),
        dummy2 = expand("outputs/metabat2_prokka_sigs/{sample}_done.txt", sample = SAMPLES)
    output: txt = "outputs/metabat2_gather_labelled_bins/{acc_db}.txt"
    conda: 'envs/tidy.yml'
    resources:
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    script: "scripts/metabat_generate_lists_of_bins_of_same_species.R"

rule metabat_intersect_sigs_for_core:
    input: "outputs/metabat2_gather_labelled_bins/{acc_db}.txt"
    output: "outputs/metabat2_prokka_sigs_core/{acc_db}_core_kmers_unnamed.sig"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sig intersect -o {output} --from-file {input} 
    """
    
rule metabat_intersect_sigs_for_core_rename:
    input:  "outputs/metabat2_prokka_sigs_core/{acc_db}_core_kmers_unnamed.sig"
    output:  "outputs/metabat2_prokka_sigs_core/{acc_db}_core_kmers.sig"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sig rename -o {output} {input} {wildcards.acc_db}_metabat2_core
    """

rule metabat_merge_species_sketches_for_all:
    input: "outputs/metabat2_gather_labelled_bins/{acc_db}.txt"
    output: "outputs/metabat2_prokka_sigs_all/{acc_db}_all_kmers.sig"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sig merge --name {wildcards.acc_db}_metabat2_all -o {output} --from-file {input} 
    """
    
rule metabat_sig_to_csv_species_for_all:
    input: "outputs/metabat2_prokka_sigs_all/{acc_db}_all_kmers.sig"
    output: "outputs/metabat2_prokka_sigs_all/{acc_db}_all_kmers.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    python scripts/sig_to_csv.py {input} {output}
    """
    
##############################################################
## Pick organisms
###############################################################

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
        lineages = "inputs/gtdb-rs202.taxonomy.v2.csv",
        hmp_metadata = "inputs/hmp2_metadata.csv",
        gather_grist = "outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv",
        pdf="figures/common_species_breakdown.pdf"
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 8000
    threads: 1
    script: "scripts/gather_gtdb_to_sgc_queries.R"

#rule download_sgc_query_genomes:
#    input: 
#        gather_grist = "outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv",
#        conf = "conf/genome-grist-conf.yml"
#    output: "genbank_genomes/genome-grist-done.txt"
#    conda: "envs/genome-grist.yml"
#    resources:
#        mem_mb = 8000
#    threads: 1
#    shell:'''
#    genome-grist run {input.conf} --until make_sgc_conf --nolock
#    touch {output}
#    '''
#    
#rule touch_sgc_query_genomes:
#    input: "genbank_genomes/genome-grist-done.txt"
#    output: "genbank_genomes/{acc}_genomic.fna.gz"
#    resources:
#        mem_mb = 500
#    threads: 1
#    shell:'''
#    ls {output}
#    '''

rule make_genome_info_csv:
    output:
        csvfile = 'outputs/genbank_genomes/{acc}.info.csv'
    conda: "envs/genbank_genomes.yml"
    resources:
        mem_mb = 8000
    threads: 1
    shell: """
        python scripts/genbank_genomes.py {wildcards.acc} \
            --output {output.csvfile}
    """

# combined info.csv
#rule make_combined_info_csv:
#    input:
#        Checkpoint_GatherResults('outputs/genbank_genomes/{acc}.info.csv')
#    output:
#        genomes_info_csv = "outputs/genbank/gather_gtdb-rs202-genomic.genomes.info.csv",
#    conda: "envs/genome-grist.yml"
#    resources:
#        mem_mb = 8000
#    threads: 1
#    shell: """
#        python scripts/combine_csvs.py {input} > {output}
#    """

# download actual genomes!
rule download_matching_genome_wc:
    input:
        csvfile = ancient('outputs/genbank_genomes/{acc}.info.csv')
    output:
        genome = "outputs/genbank_genomes/{acc}_genomic.fna.gz"
    resources:
        mem_mb = 500
    threads: 1
    run:
        with open(input.csvfile, 'rt') as infp:
            r = csv.DictReader(infp)
            rows = list(r)
            assert len(rows) == 1
            row = rows[0]
            acc = row['acc']
            assert wildcards.acc.startswith(acc)
            url = row['genome_url']
            name = row['ncbi_tax_name']

            print(f"downloading genome for acc {acc}/{name} from NCBI...",
                file=sys.stderr)
            with open(output.genome, 'wb') as outfp:
                with urllib.request.urlopen(url) as response:
                    content = response.read()
                    outfp.write(content)
                    print(f"...wrote {len(content)} bytes to {output.genome}",
                        file=sys.stderr)

   
rule make_sgc_genome_query_conf_files:
    input:
        queries = Checkpoint_GatherResults("outputs/genbank_genomes/{acc}_genomic.fna.gz"),
    output:
        conf = "outputs/sgc_conf/{sample}_k31_r1_conf.yml"
    resources:
        mem_mb = 500
    threads: 1
    run:
        # something weird is happening here where snakemake solves the inputs correctly, but 
        # "queries" ends up being the gather csv file. So...parse the CSV file for accession
        # numbers. Sigh.
        gather_csv = input[0]
        genome_accs = []
        with open(gather_csv, 'rt') as fp:
            r = csv.DictReader(fp)
            for row in r:
                acc = row['name'].split(' ')[0]
                acc = "outputs/genbank_genomes/" + acc + "_genomic.fna.gz"
                genome_accs.append(acc)
        
        query_list = "\n- ".join(genome_accs)
        #query_list = "\n- ".join(input.queries)
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
        #queries = ancient(Checkpoint_GatherResults("outputs/genbank_genomes/{acc}_genomic.fna.gz")), 
        conf = ancient("outputs/sgc_conf/{sample}_k31_r1_conf.yml"),
        reads = "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output:
        "outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/results.csv"
        #"outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz",
    params: outdir = "outputs/sgc_genome_queries"
    conda: "envs/spacegraphcats.yml"
    resources:
        mem_mb = 20000
    benchmark: "benchmarks/{sample}_sgc.tsv"
    threads: 1
    shell:'''
    python -m spacegraphcats run {input.conf} extract_reads --nolock --outdir={params.outdir} --rerun-incomplete 
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

rule sourmash_sketch_spacegraphcats_nbhds:
    input: "outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz"
    output: 'outputs/nbhd_sigs/{sample}-{acc}-dna-k31.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}-{acc}_sketch_dna-k31.tsv"
    shell:"""
    sourmash sketch dna -p k=31,scaled=2000 -o {output} --name {wildcards.acc} {input}
    """

rule sourmash_gather_spacegraphcats_nbhds:
    input:
        sig = 'outputs/nbhd_sigs/{sample}-{acc}-dna-k31.sig',
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic.k31.sbt.zip",
    output: 
        csv="outputs/nbhd_gather/{sample}-{acc}_gather_gtdb-rs202-genomic.csv",
        matches="outputs/nbhd_gather/{sample}-{acc}_gather_gtdb-rs202-genomic.matches",
        un="outputs/nbhd_gather/{sample}-{acc}_gather_gtdb-rs202-genomic.un"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}-{acc}_gather_nbhd.tsv"
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db} 
    '''

rule orpheum_translate_reads:
    input: 
        ref="inputs/orpheum_index/gtdb-rs202.protein-k10.nodegraph",
        fasta="outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz"
    output:
        pep="outputs/orpheum/{sample}-{acc}.coding.faa",
        nuc="outputs/orpheum/{sample}-{acc}.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum/{sample}-{acc}.nuc_noncoding.fna",
        csv="outputs/orpheum/{sample}-{acc}.coding_scores.csv",
        json="outputs/orpheum/{sample}-{acc}.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/{sample}-{acc}_orpheum_translate.txt"
    resources:  
        mem_mb=150000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    orpheum translate --alphabet protein --peptide-ksize 10  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fasta} > {output.pep}
    '''

rule sourmash_sketch_species_genomes:
    input: "outputs/orpheum/{sample}-{acc}.coding.faa"
    output: 'outputs/nbhd_sigs/{sample}-{acc}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}-{acc}_sketch.tsv"
    shell:"""
    sourmash sketch protein -p k=10,scaled=100,protein -o {output} --name {wildcards.acc} {input}
    """

rule convert_signature_to_csv:
    input: 'outputs/nbhd_sigs/{sample}-{acc}.sig'
    output: 'outputs/nbhd_sigs/{sample}-{acc}.csv'
    conda: 'envs/sourmash.yml'
    threads: 1
    benchmark: "benchmarks/{sample}-{acc}_sketch_to_csv.tsv"
    resources:
        mem_mb=4000,
        tmpdir = TMPDIR
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_table_long:
    input: 
        expand("outputs/nbhd_sigs/{sample}-{{acc}}.csv", sample = SAMPLES)
    output: csv = "outputs/nbhd_sketch_tables/{acc}_long.csv"
    conda: 'envs/tidy.yml'
    threads: 1
    benchmark: "benchmarks/sourmash_sketch_tables_long/{acc}.txt"
    resources:
        mem_mb=64000,
        tmpdir = TMPDIR
    script: "scripts/sketch_csv_to_long.R"
    

#######################################
## Species-specific databases
#######################################


checkpoint query_to_species_db:
    input:
        gather = "outputs/genbank/gather_gtdb-rs202-genomic.x.genbank.gather.csv",
        lineages = "inputs/gtdb-rs202.taxonomy.v2.csv",
    output: 
        csv = 'outputs/genbank/gather_gtdb-rs202-genomic.x.dbs.csv'
    conda: "envs/tidy.yml"
    threads: 1
    benchmark: "benchmarks/query_to_species_db.tsv"
    resources:
        mem_mb=4000,
        tmpdir = TMPDIR
    script: "scripts/query_to_species_db.R"

rule orpheum_translate_reads_species_db:
    input: 
        ref="inputs/orpheum_index/gtdb-rs202.{species}.protein-k10.nodegraph",
        fasta="outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{acc}_genomic.fna.gz.cdbg_ids.reads.gz"
    output:
        pep="outputs/orpheum_species/{sample}-{acc}-{species}.coding.faa",
        nuc="outputs/orpheum_species/{sample}-{acc}-{species}.nuc_coding.fna",
        nuc_noncoding="outputs/orpheum_species/{sample}-{acc}-{species}.nuc_noncoding.fna",
        csv="outputs/orpheum_species/{sample}-{acc}-{species}.coding_scores.csv",
        json="outputs/orpheum_species/{sample}-{acc}-{species}.summary.json"
    conda: "envs/orpheum.yml"
    benchmark: "benchmarks/{sample}-{acc}-{species}_orpheum_species_translate.txt"
    resources:  
        mem_mb=15000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    orpheum translate --jaccard-threshold 0.39 --alphabet protein --peptide-ksize 10  --peptides-are-bloom-filter --noncoding-nucleotide-fasta {output.nuc_noncoding} --coding-nucleotide-fasta {output.nuc} --csv {output.csv} --json-summary {output.json} {input.ref} {input.fasta} > {output.pep}
    '''

rule sourmash_sketch_species_genomes_species:
    input: "outputs/orpheum_species/{sample}-{acc_db}.coding.faa"
    output: 'outputs/nbhd_sigs_species/{sample}-{acc_db}.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}-{acc_db}_sketch.tsv"
    shell:"""
    sourmash sketch protein -p k=10,scaled=100,protein -o {output} --name {wildcards.acc_db} {input}
    """

rule convert_signature_to_csv_species:
    input: 'outputs/nbhd_sigs_species/{sample}-{acc_db}.sig'
    output: 'outputs/nbhd_sigs_species/{sample}-{acc_db}.csv'
    conda: 'envs/sourmash.yml'
    threads: 1
    benchmark: "benchmarks/{sample}-{acc_db}_sketch_to_csv.tsv"
    resources:
        mem_mb=4000,
        tmpdir = TMPDIR
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_table_long_species:
    input: 
        expand("outputs/nbhd_sigs_species/{sample}-{{acc_db}}.csv", sample = SAMPLES)
    output: csv = "outputs/nbhd_sketch_tables_species/{acc_db}_long.csv"
    conda: 'envs/tidy.yml'
    threads: 1
    benchmark: "benchmarks/sourmash_sketch_tables_long/{acc_db}.txt"
    resources:
        mem_mb=64000,
        tmpdir = TMPDIR
    script: "scripts/sketch_csv_to_long.R"
    

rule install_pagoo:
    output: pagoo = "outputs/pagoo_species/pagoo.txt"
    conda: 'envs/pagoo.yml'
    threads: 1
    resources:
        mem_mb=1000,
        tmpdir = TMPDIR
    script: "scripts/install_pagoo.R"

rule species_pagoo:
    input: 
        csv = "outputs/nbhd_sketch_tables_species/{acc_db}_long.csv",
        pagoo = "outputs/pagoo_species/pagoo.txt"
    output: 
        pca = "outputs/pagoo_species/{acc_db}_pca.pdf",
        binmap = "outputs/pagoo_species/{acc_db}_binmap.pdf"
    conda: 'envs/pagoo.yml'
    threads: 1
    benchmark: "benchmarks/pagoo_{acc_db}.txt"
    resources:
        mem_mb=8000,
        tmpdir = TMPDIR
    script: "scripts/pagoo_plt.R"

# confirm which *strains* are present in kaa-mers
rule sourmash_sketch_species_genomes_species_dna:
    input: "outputs/orpheum_species/{sample}-{acc_db}.nuc_coding.fna",
    output: 'outputs/nbhd_sigs_species/{sample}-{acc_db}-dna-k51.sig' 
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}-{acc_db}_sketch.tsv"
    shell:"""
    sourmash sketch dna -p k=51,scaled=2000 -o {output} --name {wildcards.acc_db} {input}
    """

rule sourmash_gather_orpheum_species:
    input:
        sig = 'outputs/nbhd_sigs_species/{sample}-{acc_db}-dna-k51.sig',
        db="/group/ctbrowngrp/gtdb/databases/ctb/gtdb-rs202.genomic.k51.sbt.zip",
    output: 
        csv="outputs/nbhd_species_gather/{sample}-{acc_db}_gather_gtdb-rs202-genomic.csv",
        matches="outputs/nbhd_species_gather/{sample}-{acc_db}_gather_gtdb-rs202-genomic.matches",
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 16000,
        tmpdir = TMPDIR
    threads: 1
    benchmark: "benchmarks/{sample}-{acc_db}_gather_nbhd_species.tsv"
    shell:'''
    sourmash gather -o {output.csv} --threshold-bp 0 --save-matches {output.matches} --scaled 2000 -k 51 {input.sig} {input.db} 
    '''

# create a signature for the *core* kaa-mer metapangenome
rule generate_list_of_sigs_to_intersect_based_on_kmer_threshold:
    input: csv="outputs/nbhd_sketch_tables_species/{acc_db}_long.csv"
    output: lst="outputs/nbhd_sigs_species_lists/{acc_db}.txt"
    conda: "envs/tidy.yml"
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    script: "scripts/kaa_generate_list_of_sigs_to_intersect.R"

rule intersect_species_sketches_for_core:
    input: "outputs/nbhd_sigs_species_lists/{acc_db}.txt"
    output: 'outputs/nbhd_sigs_species_core/{acc_db}_core_kmers_unnamed.sig'
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sig intersect -o {output} --from-file {input} 
    """

rule rename_species_sketches_for_core:
    input: 'outputs/nbhd_sigs_species_core/{acc_db}_core_kmers_unnamed.sig'
    output: 'outputs/nbhd_sigs_species_core/{acc_db}_core_kmers.sig'
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sig rename -o {output} {input} {wildcards.acc_db}_kaa_core
    """

# create a signature for the *entire* kaa-mer metapangenome
rule merge_species_sketches_for_all:
    input: expand('outputs/nbhd_sigs_species/{sample}-{{acc_db}}.sig', sample = SAMPLES)
    output: 'outputs/nbhd_sigs_species_all/{acc_db}_all_kmers.sig'
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sig merge --name {wildcards.acc_db}_kaa_all -o {output} {input} 
    """

rule species_sketches_sig_to_csv_for_all:
    input: 'outputs/nbhd_sigs_species_all/{acc_db}_all_kmers.sig'
    output: 'outputs/nbhd_sigs_species_all/{acc_db}_all_kmers.csv'
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    python scripts/sig_to_csv.py {input} {output}
    """
    
######################################################
## Compare species-level db and gtdb-level db results
######################################################

rule sourmash_compare_orpheum_outputs:
    input:
        Checkpoint_AccToDbs("outputs/nbhd_sigs_species/{{sample}}-{acc_db}.sig"),
        Checkpoint_GatherResults('outputs/nbhd_sigs/{{sample}}-{acc}.sig') 
    output: 
        comp = "outputs/orpheum_compare/{sample}_comp",
        csv = "outputs/orpheum_compare/{sample}_comp.csv"
    conda: "envs/sourmash.yml"
    benchmark: "benchmarks/{sample}_sourmash_comp_orpheum_outputs.txt"
    resources:  
        mem_mb=15000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    sourmash compare -o {output.comp} --csv {output.csv} {input}
    '''

rule sourmash_compare_orpheum_outputs_max:
    input:
        Checkpoint_AccToDbs("outputs/nbhd_sigs_species/{{sample}}-{acc_db}.sig"),
        Checkpoint_GatherResults('outputs/nbhd_sigs/{{sample}}-{acc}.sig') 
    output: 
        comp = "outputs/orpheum_compare/{sample}_comp_containment",
        csv = "outputs/orpheum_compare/{sample}_comp_containment.csv"
    conda: "envs/sourmash.yml"
    benchmark: "benchmarks/{sample}_sourmash_comp_containment_orpheum_outputs.txt"
    resources:  
        mem_mb=15000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    sourmash compare --max-containment -o {output.comp} --csv {output.csv} {input}
    '''

##############################################################
## Build reference pangenome to compare kmer nbhds against
##############################################################

# Note reference pangenomes were copied over from 2021-panmers,
# to get around wildcard issues solving for species, and all 
# genomes underneath that species.


rule extract_pangenome_fasta_names:
    input: 'outputs/roary/{species}/pan_genome_reference.fa',
    output: 'outputs/roary/{species}/pan_genome_reference_names.txt'
    threads: 1
    resources: 
        mem_mb= 4000,
        tmpdir = TMPDIR
    shell:'''
    grep ">" {input} > {output}
    '''

rule identify_core_gene_sequence_names:
    input: 
        pg="outputs/pagoo_species/pagoo.txt",
        pa='outputs/roary/{acc_db}/gene_presence_absence.csv', 
        names='outputs/roary/{acc_db}/pan_genome_reference_names.txt'
    output: core_names="outputs/roary/{acc_db}/pan_genome_reference_core_gene_names.txt"
    threads: 1
    conda: "envs/pagoo.yml"
    resources: 
        mem_mb= 16000,
        tmpdir = TMPDIR
    script: "scripts/identify_core_gene_sequence_names.R"

rule extract_core_gene_sequences_from_pangenome:
    input:
        fa="outputs/roary/{acc_db}/pan_genome_reference.fa",
        core_names="outputs/roary/{acc_db}/pan_genome_reference_core_gene_names.txt"
    output: "outputs/roary/{acc_db}/pan_genome_reference_core_genes.fa"
    conda: "envs/seqtk.yml"
    resources:
        mem_mb = 2000,
        tmpdir=TMPDIR
    threads: 1
    shell:'''
    seqtk subseq {input.fa} {input.core_names} > {output}
    '''

rule translate_core_gene_sequences_from_pangenome:
    input: "outputs/roary/{acc_db}/pan_genome_reference_core_genes.fa"
    output: "outputs/roary/{acc_db}/pan_genome_reference_core_genes.faa"
    conda: 'envs/emboss.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 2
    shell:'''
    transeq {input} {output}
    '''

rule sketch_core_gene_sequences_from_pangenome:
    input: "outputs/roary/{acc_db}/pan_genome_reference_core_genes.faa"
    output: "outputs/roary_sigs_core/{acc_db}_pan_genome_reference_core_genes.sig"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sketch protein -p k=10,scaled=100,protein -o {output} --name {wildcards.acc_db}_roary_core {input}
    """
 
rule sketch_all_gene_sequences_from_pangenome:
    input: "outputs/roary/{acc_db}/pan_genome_reference_core_genes.faa"
    output: "outputs/roary_sigs_all/{acc_db}_pan_genome_reference_all_genes.sig"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash sketch protein -p k=10,scaled=100,protein -o {output} --name {wildcards.acc_db}_roary_all {input}
    """
 
rule roary_sketches_sig_to_csv_for_all:
    input: "outputs/roary_sigs_all/{acc_db}_pan_genome_reference_all_genes.sig"
    output: "outputs/roary_sigs_all/{acc_db}_pan_genome_reference_all_genes.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    python scripts/sig_to_csv.py {input} {output}
    """
    
#############################################################
## compare containment/similarity across core sequences
#############################################################

rule compare_similarity_pangenomes_core: 
    input:
        kmer='outputs/nbhd_sigs_species_core/{acc_db}_core_kmers.sig',
        gene="outputs/roary_sigs_core/{acc_db}_pan_genome_reference_core_genes.sig",
        bins="outputs/metabat2_prokka_sigs_core/{acc_db}_core_kmers.sig"
    output: "outputs/compare_core/{acc_db}_similarity.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash compare --csv {output} {input}
    """

rule compare_containment_pangenomes_core: 
    input:
        kmer='outputs/nbhd_sigs_species_core/{acc_db}_core_kmers.sig',
        gene="outputs/roary_sigs_core/{acc_db}_pan_genome_reference_core_genes.sig",
        bins="outputs/metabat2_prokka_sigs_core/{acc_db}_core_kmers.sig"
    output: "outputs/compare_core/{acc_db}_containment.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash compare --containment --csv {output} {input}
    """

rule compare_maxcontainment_pangenomes_core: 
    input:
        kmer='outputs/nbhd_sigs_species_core/{acc_db}_core_kmers.sig',
        gene="outputs/roary_sigs_core/{acc_db}_pan_genome_reference_core_genes.sig",
        bins="outputs/metabat2_prokka_sigs_core/{acc_db}_core_kmers.sig"
    output: "outputs/compare_core/{acc_db}_maxcontainment.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash compare --max-containment --csv {output} {input}
    """
    
    
###########################################
## Compare whole pangenomes
###########################################

    
rule compare_similarity_pangenomes_all: 
    input:
        kmer='outputs/nbhd_sigs_species_all/{acc_db}_all_kmers.sig',
        gene="outputs/roary_sigs_all/{acc_db}_pan_genome_reference_all_genes.sig",
        bins="outputs/metabat2_prokka_sigs_all/{acc_db}_all_kmers.sig"
    output: "outputs/compare_all/{acc_db}_similarity.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash compare --csv {output} {input}
    """

rule compare_containment_pangenomes_all: 
    input:
        kmer='outputs/nbhd_sigs_species_all/{acc_db}_all_kmers.sig',
        gene="outputs/roary_sigs_all/{acc_db}_pan_genome_reference_all_genes.sig",
        bins="outputs/metabat2_prokka_sigs_all/{acc_db}_all_kmers.sig"
    output: "outputs/compare_all/{acc_db}_containment.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash compare --containment --csv {output} {input}
    """

rule compare_maxcontainment_pangenomes_all: 
    input:
        kmer='outputs/nbhd_sigs_species_all/{acc_db}_all_kmers.sig',
        gene="outputs/roary_sigs_all/{acc_db}_pan_genome_reference_all_genes.sig",
        bins="outputs/metabat2_prokka_sigs_all/{acc_db}_all_kmers.sig"
    output: "outputs/compare_all/{acc_db}_maxcontainment.csv"
    conda: 'envs/sourmash.yml'
    resources:
        mem_mb = 4000,
        tmpdir = TMPDIR
    threads: 1
    shell:"""
    sourmash compare --max-containment --csv {output} {input}
    """


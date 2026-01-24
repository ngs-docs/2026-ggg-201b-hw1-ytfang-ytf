SAMPLES = ["SRR2584403_1", "SRR2584404_1", "SRR2584405_1", "SRR2584857_1"]
GENOME = ["ecoli-rel606"]

rule make_vcf:
    input:
        expand("outputs/{sample}.x.{genome}.vcf",
               sample=SAMPLES, genome=GENOME),
        expand("outputs/{sample}.x.{genome}.vep.txt",
              sample=SAMPLES, genome=GENOME),
  
rule uncompress_genome:
    input: "{genome}.fa.gz"
    output: "outputs/{genome}.fa"
    shell: """
        gunzip -c {input} > {output}
    """

rule map_reads:
    input:
        reads="{sample}.fastq.gz",
        ref="outputs/{genome}.fa"
    output: "outputs/{sample}.x.{genome}.sam"
    conda: "mapping"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """

rule sam_to_bam:
    input: "outputs/{sample}.x.{genome}.sam"
    output: "outputs/{sample}.x.{genome}.bam",
    conda: "mapping"
    shell: """
        samtools view -b {input} > {output}
     """

rule sort_bam:
    input: "outputs/{sample}.x.{genome}.bam"
    output: "outputs/{sample}.x.{genome}.bam.sorted"
    conda: "mapping"
    shell: """
        samtools sort {input} > {output}
    """

rule index_bam:
    input: "outputs/{sample}.x.{genome}.bam.sorted"
    output: "outputs/{sample}.x.{genome}.bam.sorted.bai"
    conda: "mapping"
    shell: """
        samtools index {input}
    """

rule call_variants:
    input:
        ref="outputs/{genome}.fa",
        bam="outputs/{sample}.x.{genome}.bam.sorted",
        bai="outputs/{sample}.x.{genome}.bam.sorted.bai",
    output:
        pileup="outputs/{sample}.x.{genome}.pileup",
        bcf="outputs/{sample}.x.{genome}.bcf",
        vcf="outputs/{sample}.x.{genome}.vcf",
    conda: "mapping"
    shell: """
        bcftools mpileup -Ou -f {input.ref} {input.bam} > {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view {output.bcf} > {output.vcf}
    """
    
rule tabix:
    input:
        gff="{filename}.gff.gz",
    output:
        tabix_idx='{filename}.gff.gz.tbi',
    shell: """
        tabix {input}
    """

rule predict_effects:
    input:
        fasta="{genome}.fa.gz",
        gff="{genome}.sorted.gff.gz",
        vcf="outputs/{sample}.x.{genome}.vcf",
        tabix_idx='ecoli-rel606.sorted.gff.gz.tbi',
    output:
        txt="outputs/{sample}.x.{genome}.vep.txt",
        html="outputs/{sample}.x.{genome}.vep.txt_summary.html",
        warn="outputs/{sample}.x.{genome}.vep.txt_warnings.txt",
    conda: "vep"
    shell: """
       vep --fasta {input.fasta} --gff {input.gff} -i {input.vcf} -o {output.txt}
    """    
    
    

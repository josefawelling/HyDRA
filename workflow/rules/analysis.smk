from pathlib import Path

'''rule rename_assembly:
    input:
        "results/final_assemblies/{strain}/assembly.fasta"
    output:
        "results/assembly/{strain}_contig.fasta"
    shell:
        "cp {input} {output}"'''

rule prokka:
    input:
        "results/assembly/{strain}/{strain}_assembly.fasta"
        #"results/assembly/{strain}_contig.fasta"
    output:
        gff = "results/analysis/prokka/{strain}/{strain}.gff"
    params:
        outdir = lambda wildcards, output: Path(output.gff).parent
    log:
        "logs/prokka/{strain}.log"
    conda:
        "../envs/prokka.yaml"
    shell:
        "prokka --outdir {params.outdir}/ --force --prefix {wildcards.strain} {input} 2> {log}"

rule abricate:
    input:
        "results/assembly/{strain}/{strain}_assembly.fasta"
        #"results/assembly/{strain}_contig.fasta"
    output:
        all = "results/analysis/abricate/{strain}/{strain}_all.csv",
        info = "results/analysis/abricate/{strain}/{strain}_info.txt"
    log:
        "logs/abricate/{strain}_run.log"
    conda:
        "../envs/abricate.yaml"
    script:
        "../scripts/abricate_run.py"
        
rule abricate_filter:
    input:
        "results/analysis/abricate/{strain}/{strain}_all.csv"
    output:
        args = "results/analysis/abricate/{strain}/{strain}_ARGs.csv",
        filtout = "results/analysis/abricate/{strain}/{strain}_filtoutARGs.csv"
    log:
        "logs/abricate/{strain}_filter.log"
    conda:
        "../envs/abricate.yaml"
    script:
        "../scripts/abricate_filter.py"

## visualize the annotated genome
rule gview:
    input:
        gff = "results/analysis/prokka/{strain}/{strain}.gff",
        gbk = "results/analysis/prokka/{strain}/{strain}.gbk"
    output:
        "results/reports/assembly/{strain}_gview.png"
    params:
        extra = "-l circular",
        jar_file = "resources/gview/gview.jar",
        ## style sheet
        gss = "resources/gview/example_styles/gssExample.gss"
        ## we need a style sheet to make this look cool! (see basic and gss examples)
        #style_sheet = 
    log:
        "logs/gview/{strain}.log"
    conda:
        "../envs/gview.yaml" # env with java
    shell:
        "java -jar {params.jar_file} -i {input.gbk} -o {output} "
        "-g {input.gff} -s {params.gss} {params.extra}"
  
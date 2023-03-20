version 1.0

workflow tobias {
    meta {
        version: "v0.1"
        author: "Mei Knudson; tools from Mette Bensen at Looso Lab https://github.com/loosolab/TOBIAS"
        description: "Workflow to run TOBIAS using outputs from ENCODE bulk ATAC-seq processing pipeline"
    }

    input {
        File bam
        File genome
        File genome_index
        File peaks
        File blacklist
        File motifs
        String prefix
        String docker_image = "mknudson/tobias"
        Int cores = 8
    }

    parameter_meta {
        bam: {
            description: "Unsorted BAM file (.bam)"
        }
        genome: {
            description: "Genome FASTA file (.fasta)"
        }
        genome_index: {
            description: "Genome index file (.fasta.fai)"
        }
        peaks: {
            description: "Narrowpeak peak file (.narrowpeak.gz)"
        }
        blacklist: {
            description: "Genome blacklist file (.bed)"
        }
        motifs: {
            description: "Motif file (PFM, JASPAR, or MEME format)"
        }
    }

    call ATACorrect {
        input:
            bam = bam,
            genome = genome,
            genome_index = genome_index,
            prefix = prefix,
            peaks = peaks,
            blacklist = blacklist,
            docker_image = docker_image,
            cores = cores
    }

    call ScoreBigwig {
        input:
            corrected_bw = ATACorrect.corrected_bw,
            peaks = ATACorrect.bed,
            prefix = prefix,
            docker_image = docker_image,
            cores = cores
    }

    call BINDetect {
        input:
            footprint_bw = ScoreBigwig.footprint_bw,
            motifs = motifs,
            genome = genome,
            genome_index = genome_index,
            peaks = ATACorrect.bed, 
            prefix = prefix,
            docker_image = docker_image,     
            cores = cores
    }

    output {
        File corrected_bw = ATACorrect.corrected_bw
        File footprint_bw = ScoreBigwig.footprint_bw
        File bindetect_results_txt = BINDetect.bindetect_results_txt
        File bindetect_results_xlsx = BINDetect.bindetect_results_xlsx
        File bindetect_figures = BINDetect.bindetect_figures
    }
}

task ATACorrect {
    input {
        File bam
        File genome
        File genome_index
        File peaks
        File blacklist
        String prefix
        String outdir = "ATACorrect"
        String docker_image
        Int cores
    }

    String sorted_bam = "~{prefix}.sorted.bam"
    String bed = "~{prefix}.bed"

    command <<<
        set -e

        # Sort bam prior to indexing
        samtools sort ~{bam} -o ~{sorted_bam}

        # Index bam
        samtools index -@ ~{cores} ~{sorted_bam} ~{sorted_bam}.bai

        # Unzip narrowpeak file and convert to bed
        gunzip -c ~{peaks} | cut -f 1-6 > ~{bed}

        $(which TOBIAS) ATACorrect \
            --bam ~{sorted_bam} \
            --genome ~{genome} \
            --peaks ~{bed} \
            --blacklist ~{blacklist} \
            --outdir ~{outdir} \
            --cores ~{cores} \
            --prefix ~{prefix}
    >>>

    output {
        File bed = "~{bed}"
        File uncorrected_bw = "~{outdir}/~{prefix}_uncorrected.bw"
        File bias_bw = "~{outdir}/~{prefix}_bias.bw"
        File expected_bw = "~{outdir}/~{prefix}_expected.bw"
        File corrected_bw = "~{outdir}/~{prefix}_corrected.bw"
        File atacorrect_pdf = "~{outdir}/~{prefix}_atacorrect.pdf"
    }

    runtime {
        cpu: cores
        docker: docker_image
        memory: "16 GB"
        disks: "local-disk 50 LOCAL"
    }
}

task ScoreBigwig {
    input {
        File corrected_bw
        File peaks
        String prefix
        String docker_image
        Int cores
    }

    String footprints = "~{prefix}.footprints.bw"

    command <<<
        set -e

        $(which TOBIAS) FootprintScores \
            --signal ~{corrected_bw} \
            --regions ~{peaks} \
            --output ~{footprints} \
            --cores ~{cores}
    >>>

    output {
        File footprint_bw = "~{footprints}"
    }

    runtime {
        cpu: cores
        docker: docker_image
        disks: "local-disk 50 LOCAL"
    }
}

task BINDetect {
    input {
        File footprint_bw
        File motifs
        File genome
        File genome_index
        File peaks
        String prefix
        String outdir = "BINDetect"
        String docker_image
        Int cores
    }

    command <<<
        set -e

        $(which TOBIAS) BINDetect \
            --motifs ~{motifs} \
            --signals ~{footprint_bw} \
            --genome ~{genome} \
            --peaks ~{peaks} \
            --prefix ~{prefix} \
            --outdir ~{outdir} \
            --cores ~{cores}
    >>>

    output {
        File bindetect_results_txt = "~{outdir}/~{prefix}_bindetect_results.txt"
        File bindetect_results_xlsx = "~{outdir}/~{prefix}_bindetect_results.xlsx"
        File bindetect_figures = "~{outdir}/~{prefix}_bindetect_figures.pdf"
        File bindetect_distances = "~{outdir}/~{prefix}_bindetect_distances.txt"
    }

    runtime {
        cpu: cores
        docker: docker_image
        disks: "local-disk 50 LOCAL"
    }
}


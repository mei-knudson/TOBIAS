version 1.0

workflow tobias {
    input {
        File bam
        File genome
        File peaks
        File blacklist
        File motifs
        String prefix
        String docker_image = "mknudson/tobias"
        Int cores = 8
    }

    call ATACorrect {
        input:
            bam = bam,
            genome = genome,
            peaks = peaks,
            blacklist = blacklist,
            docker_image = docker_image,
            cores = cores
    }

    call ScoreBigwig {
        input:
            corrected_bw = ATACorrect.corrected_bw,
            peaks = peaks,
            docker_image = docker_image,
            cores = cores
    }

    call BINDetect {
        input:
            footprint_bw = ScoreBigwig.footprint_bw,
            motifs = motifs,
            genome = genome,
            peaks = peaks, 
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
        File peaks
        File blacklist
        String prefix
        String outdir = "ATACorrect"
        String docker_image
        Int cores
    }

    String sorted_bam = "~{prefix}.sorted.bam"

    command <<<
        set -e

        samtools sort -n ~{bam} -o ~{sorted_bam}

        samtools index -@ ~{cores} ~{sorted_bam} ~{sorted_bam}.bai

        $(which TOBIAS) ATACorrect \
            --bam ~{sorted_bam} \
            --genome ~{genome} \
            --peaks ~{peaks} \
            --blacklist ~{blacklist} \
            --outdir ~{outdir} \
            --cores ~{cores} \
            --prefix ~{prefix}
    >>>

    output {
        File sorted_bam = "~{sorted_bam}"
        File uncorrected_bw = "~{outdir}/~{prefix}_uncorrected.bw"
        File bias_bw = "~{outdir}/~{prefix}_bias.bw"
        File expected_bw = "~{outdir}/~{prefix}_expected.bw"
        File corrected_bw = "~{outdir}/~{prefix}_corrected.bw"
        File atacorrect_pdf = "~{outdir}/~{prefix}_atacorrect.pdf"
    }

    runtime {
        cpu: cores
        docker: docker_image
    }
}

task ScoreBigwig {
    input {
        File corrected_bw
        File peaks
        String outfile = "footprints.bw"
        String docker_image
        Int cores
    }

    command <<<
        set -e

        $(which TOBIAS) FootprintScores \
            --signal ~{corrected_bw} \
            --regions ~{peaks} \
            --output ~{outfile} \
            --cores ~{cores}
    >>>

    output {
        File footprint_bw = "~{outfile}"
    }

    runtime {
        cpu: cores
        docker: docker_image
    }
}

task BINDetect {
    input {
        File footprint_bw
        File motifs
        File genome
        File peaks
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
            --outdir ~{outdir} \
            --cores ~{cores}
    >>>

    output {
        File bindetect_results_txt = glob('~{outdir}/*_bindetect_results.txt')[0]
        File bindetect_results_xlsx = glob('~{outdir}/*_bindetect_results.xlsx')[0]
        File bindetect_figures = glob('~{outdir}/bindetect_figures.pdf')[0]
        File distance_matrix = glob('~{outdir}/TF_distance_matrix.txt')[0]
        Array[File] tf_files = glob('~{outdir}/*/*')
    }

    runtime {
        cpu: cores
        docker: docker_image
    }
}

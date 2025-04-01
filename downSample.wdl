version 1.0

struct fastqPair {
    File fastq1
    File fastq2
}

struct bamFile {
  File bam
  File bamIndex
}

workflow downSample {
    input {
        fastqPair? inputFastq
        bamFile? inputBam
        String outputFileNamePrefix
        String downSampleTool
        String downSampleMethod
        Float? downSampleRatio
        Int? downSampleReads
        Int? randomSampleSeed
        Boolean DoSorting = true
        Boolean createIndex = true
        Boolean checkCoverage = false
    }
    parameter_meta {
        inputFastq: "the input fastq file pair"
        inputBam: "the input bam and index"
        outputFileNamePrefix: "Prefix of output file name"
        downSampleTool: "the tool to be used in downsampling, a few options available"
        downSampleMethod: "choose between random/top_reads"
        downSampleRatio: "given a ratio for downsampled reads"
        downSampleReads: "given a number of reads after down sample"
        randomSampleSeed: "the seed for random sampling"
        DoSorting: "whether do sorting after downsample for bam file"
        createIndex: "whether create index for downsampled bam"
        checkCoverage: "whether check coverage for downsampled bam"
    }

    if ( defined (inputFastq) ) {
        call downSampleFastq {
            input:
                fastq1 = select_first([inputFastq]).fastq1,
                fastq2 = select_first([inputFastq]).fastq2,
                outputFileNamePrefix = outputFileNamePrefix,
                downSampleTool = downSampleTool,
                downSampleMethod = downSampleMethod,
                randomSampleSeed = select_first([randomSampleSeed, 0]),
                downSampleRatio = select_first([downSampleRatio, 0]),
                downSampleReads = select_first([downSampleReads, 0])
        }
    }
    if  ( defined (inputBam)) {
        call downSampleBam {
            input:
                bam = select_first([inputBam]).bam,
                bai = select_first([inputBam]).bamIndex,
                downSampleTool = downSampleTool,
                outputFileNamePrefix = outputFileNamePrefix,
                downSampleMethod = downSampleMethod,
                randomSampleSeed = select_first([randomSampleSeed, 0]),
                downSampleRatio = select_first([downSampleRatio, 0]),
                downSampleReads = select_first([downSampleReads, 0]),
                DoSorting = DoSorting,
                createIndex = createIndex,
                checkCoverage = checkCoverage
        }
    }
    meta {
        author: "Gavin Peng"
        email: "gpeng@oicr.on.ca"
        description: "Workflow to downsample fastq or bam files. Can use a combination of differrent method, tools and parameters. Notes: 1) downsample method can choose between random and top_reads, the later can only applied to fastq inputs; 2) fastq downsample tools include seqtk, seqkit; bam downsample tools include samtools, picard; 3) for seqkit, samtools, picard, prefered parameter is downSampleRatio, as use downSampleReads resulting number of reads is not exact, and may include extra compute time"
        dependencies: [
        {
            name: "seqtk/1.3",
            url: "https://github.com/lh3/seqtk"
        },
        {
            name: "seqkit/2.3.1",
            url: "https://github.com/shenwei356/seqkit"
        },
        {
            name: "picard/3.1.0",
            url: "https://broadinstitute.github.io/picard/"
        },
        { 
          name: "samtools/1.16.1",
          url: "https://github.com/samtools/samtools/releases/"
        }
      ]
      output_meta: {
        downSampledFastq1: {
            description: "Output downsampled fastq read1",
            vidarr_label: "downSampledFastq1"
        },
        downSampledFastq2: {
            description: "Output downsampled fastq read2",
            vidarr_label: "downSampledFastq1"
        },
        downSampledBam: {
            description: "Downsampled bam file",
            vidarr_label: "downSampledBam"
        },
        downSampledBai: {
            description: "Downsampled bam file index",
            vidarr_label: "downSampledBai"
        }
      }
    }
    output {
        File? downSampledFastq1 = downSampleFastq.downSampledFastq1
        File? downSampledFastq2 = downSampleFastq.downSampledFastq2
        File? downSampledBam = downSampleBam.downSampledBam
        File? downSampledBai = downSampleBam.downSampledBai 
    }
} 

task downSampleFastq {
    input {
        File fastq1
        File fastq2
        String outputFileNamePrefix
        String downSampleTool
        String downSampleMethod
        Float downSampleRatio
        Int downSampleReads
        Int randomSampleSeed
        Int timeout = 12
        Int memory = 24
        Int threads = 8
        String modules = "seqtk/1.3 seqkit/2.3.1"
    }

    parameter_meta {
        fastq1: "input fastq read1"
        fastq2: "input fastq read2"
        outputFileNamePrefix: "Prefix of output file name"
        downSampleTool: "the tool to be used in downsampling, a few options available"
        downSampleMethod: "choose between random/top_reads"
        downSampleRatio: "given a ratio for downsampled reads"
        downSampleReads: "given a number of reads after down sample"
        randomSampleSeed: "the seed for random sampling"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
    }
    String output_suffix = if (sub(basename(fastq1), ".*\\.gz$", "") != basename(fastq1)) then "fastq.gz" else "fastq"

    command <<<
        valid_downSampleTool=("seqtk/1.3" "seqkit/2.3.1")
        valid_downSampleMethod=("random" "top_reads")

        if [[ ! " ${valid_downSampleTool[@]} " =~  ~{downSampleTool} ]]; then
            echo "valid downSampleTool values are seqtk/1.3 or seqkit/2.3.1" >&2
            exit 1
        fi
        if [[ ! " ${valid_downSampleMethod[@]} " =~  ~{downSampleMethod} ]]; then
            echo "valid downSampleMethod values are random or top_reads" >&2
            exit 1
        fi
        if [ ~{randomSampleSeed} != 0 ]; then
                seed=~{randomSampleSeed}
            else
                seed=42
        fi

        if [ ~{downSampleMethod} = "top_reads" ]; then
            if [ ~{downSampleReads} -ne 0 ]; then
                if [[ ~{output_suffix} == "fastq" ]]; then 
                    head -n $((4 * ~{downSampleReads}))   ~{fastq1}  > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                    head -n $((4 * ~{downSampleReads}))   ~{fastq2}  > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
                else
                    zcat ~{fastq1} | head -n $((4 * ~{downSampleReads})) | gzip  > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                    zcat ~{fastq2} | head -n $((4 * ~{downSampleReads})) | gzip  > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
                fi

            else 
                echo "when downSampleMethod is top_reads, downSampleReads needs provided" >&2
                exit 1
            fi

        elif [ ~{downSampleTool} = "seqtk" ]; then
            #if both downSampleReads and downSampleRatio provided then will use downSampleReads
            if [ ~{downSampleReads} != 0 ]; then
                downSampleFactor=~{downSampleReads}
            elif (( $(echo "~{downSampleRatio} != 0" | bc -l) )); then
                downSampleFactor=~{downSampleRatio}
            else
                echo "downSampleReads or downSampleRatio need provided" >&2
                exit 1
            fi

            if [[ ~{output_suffix} == "fastq" ]]; then 
                seqtk sample -s${seed} ~{fastq1} ${downSampleFactor} > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                seqtk sample -s${seed} ~{fastq2} ${downSampleFactor} > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
            else
                seqtk sample -s${seed} ~{fastq1} ${downSampleFactor} | gzip > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                seqtk sample -s${seed} ~{fastq2} ${downSampleFactor} | gzip > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
            fi

        elif [ ~{downSampleTool} = "seqkit" ]; then
            #if both downSampleReads and downSampleRatio provided then will use downSampleRatio
            if (( $(echo "~{downSampleRatio} != 0" | bc -l) )); then
                if [[ ~{output_suffix} == "fastq" ]]; then 
                    seqkit sample ~{fastq1} -p ~{downSampleRatio} -s ${seed}  > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                    seqkit sample ~{fastq2} -p ~{downSampleRatio} -s ${seed}  > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
                else
                    seqkit sample ~{fastq1} -p ~{downSampleRatio} -s ${seed} | gzip  > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                    seqkit sample ~{fastq2} -p ~{downSampleRatio} -s ${seed} | gzip  > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
                fi
            elif [ ~{downSampleReads} -ne 0 ]; then
                if [[ ~{output_suffix} == "fastq" ]]; then 
                    seqkit sample -p 0.2 ~{fastq1} -s ${seed} | seqkit head -n ~{downSampleReads} > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                    seqkit sample -p 0.2 ~{fastq2} -s ${seed} | seqkit head -n ~{downSampleReads} > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
                else
                    seqkit sample -p 0.2 ~{fastq1} -s ${seed} | seqkit head -n ~{downSampleReads} | gzip > ~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}
                    seqkit sample -p 0.2 ~{fastq2} -s ${seed} | seqkit head -n ~{downSampleReads} | gzip > ~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}
                fi
            else
                echo "downSampleReads or downSampleRatio need provided" >&2
                exit 1
            fi
        fi
    >>>

    output {
        File? downSampledFastq1 = "~{outputFileNamePrefix}.downSampledFastq1.~{output_suffix}"
        File? downSampledFastq2 = "~{outputFileNamePrefix}.downSampledFastq2.~{output_suffix}"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}

task downSampleBam {
    input {
        File bam
        File bai
        String outputFileNamePrefix
        String downSampleTool
        String downSampleMethod
        Float downSampleRatio
        Int downSampleReads
        Int randomSampleSeed
        Boolean DoSorting = true
        Boolean createIndex = true
        Boolean checkCoverage = false
        Int timeout = 12
        Int memory = 24
        Int threads = 8
        String modules = "samtools/1.16.1 picard/3.1.0"
    }

    parameter_meta {
        bam: "the input bam file"
        bai: "the input bam index"
        outputFileNamePrefix: "Prefix of output file name"
        downSampleTool: "the tool to be used in downsampling, a few options available"
        downSampleMethod: "choose between random/top_reads"
        downSampleRatio: "given a ratio for downsampled reads"
        downSampleReads: "given a number of reads after down sample, because of probabilistic sampling, the reads number is not exact"
        randomSampleSeed: "the seed for random sampling"
        timeout: "The hours until the task is killed"
        memory: "The GB of memory provided to the task"
        threads: "The number of threads the task has access to"
        modules: "The modules that will be loaded"
        DoSorting: "whether do sorting after downsample for bam file"
        createIndex: "whether create index for downsampled bam"
        checkCoverage: "whether check coverage for downsampled bam"
    }

    command <<<
        valid_downSampleTool=("samtools", "picard")
        if [[ ! " ${valid_downSampleTool[@]} " =~  ~{downSampleTool} ]]; then
            echo "valid downSampleTool values are samtools or picard for bam dowmsampling" >&2
            exit 1
        fi
        if [[ ~{downSampleMethod} != "random" ]]; then
            echo "valid downSampleMethod for bam downsampling is random" >&2
            exit 1
        fi
        if [ ~{randomSampleSeed} != 0 ]; then
                seed=~{randomSampleSeed}
            else
                seed=42
        fi
        
        if [ ~{downSampleTool} = "samtools" ]; then
            #if both downSampleReads and downSampleRatio provided then will use downSampleRatio
            if (( $(echo "~{downSampleRatio} != 0" | bc -l) )); then
                downSampleFactor=~{downSampleRatio}
            elif [ ~{downSampleReads} != 0 ]; then
                TOTAL_READS=$(samtools view -c ~{bam})
                downSampleFactor=$(echo "scale=6; ~{downSampleReads} / $TOTAL_READS" | bc)
            else
                echo "downSampleReads or downSampleRatio need provided" >&2
                exit 1
            fi
            samtools view -s ${downSampleFactor} -b ~{bam} > ~{outputFileNamePrefix}.downsampled.bam
            samtools index ~{outputFileNamePrefix}.downsampled.bam

        elif [ ~{downSampleTool} = "picard" ]; then
            #if both downSampleReads and downSampleRatio provided then will use downSampleRatio
            if (( $(echo "~{downSampleRatio} != 0" | bc -l) )); then
                downSampleFactor=~{downSampleRatio}
            elif [ ~{downSampleReads} != 0 ]; then
                TOTAL_READS=$(samtools view -c ~{bam})
                downSampleFactor=$(echo "scale=6; ~{downSampleReads} / $TOTAL_READS" | bc)
            else
                echo "downSampleReads or downSampleRatio need provided" >&2
                exit 1
            fi
            export JAVA_OPTS="-Xmx$(echo "scale=0; ~{memory} * 0.8 / 1" | bc)G"
            java -jar ${PICARD_ROOT}/picard.jar DownsampleSam \
            I=~{bam} \
            O=~{outputFileNamePrefix}.downsampled.bam \
            P=${downSampleFactor} \
            RANDOM_SEED=${seed} \
            CREATE_INDEX=true \
            VALIDATION_STRINGENCY=SILENT
            mv ~{outputFileNamePrefix}.downsampled.bai ~{outputFileNamePrefix}.downsampled.bam.bai
        fi
    >>>

    output {
        File? downSampledBam = "~{outputFileNamePrefix}.downsampled.bam"
        File? downSampledBai = "~{outputFileNamePrefix}.downsampled.bam.bai"
    }

    runtime {
        modules: "~{modules}"
        memory:  "~{memory} GB"
        cpu:     "~{threads}"
        timeout: "~{timeout}"
    }
}
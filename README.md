# downSample

Workflow to downsample fastq or bam files. Can use a combination of differrent method, tools and parameters. Notes: 1) downsample method can choose between random and top_reads, the later can only applied to fastq inputs; 2) fastq downsample tools include seqtk, seqkit; bam downsample tools include samtools, picard; 3) for seqkit, samtools, picard, prefered parameter is downSampleRatio, as use downSampleReads resulting number of reads is not exact, and may include extra compute time

## Overview

## Dependencies

* [seqtk 1.3](https://github.com/lh3/seqtk)
* [seqkit 2.3.1](https://github.com/shenwei356/seqkit)
* [picard 3.1.0](https://broadinstitute.github.io/picard/)
* [samtools 1.16.1](https://github.com/samtools/samtools/releases/)


## Usage

### Cromwell
```
java -jar cromwell.jar run downSample.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`outputFileNamePrefix`|String|Prefix of output file name
`downSampleTool`|String|the tool to be used in downsampling, a few options available
`downSampleMethod`|String|choose between random/top_reads


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`inputFastq`|fastqPair?|None|the input fastq file pair
`inputBam`|bamFile?|None|the input bam and index
`downSampleRatio`|Float?|None|given a ratio for downsampled reads
`downSampleReads`|Int?|None|given a number of reads after down sample
`randomSampleSeed`|Int?|None|the seed for random sampling
`DoSorting`|Boolean|true|whether do sorting after downsample for bam file
`createIndex`|Boolean|true|whether create index for downsampled bam
`checkCoverage`|Boolean|false|whether check coverage for downsampled bam


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`downSampleFastq.timeout`|Int|12|The hours until the task is killed
`downSampleFastq.memory`|Int|24|The GB of memory provided to the task
`downSampleFastq.threads`|Int|8|The number of threads the task has access to
`downSampleFastq.modules`|String|"seqtk/1.3 seqkit/2.3.1"|The modules that will be loaded
`downSampleBam.timeout`|Int|12|The hours until the task is killed
`downSampleBam.memory`|Int|24|The GB of memory provided to the task
`downSampleBam.threads`|Int|8|The number of threads the task has access to
`downSampleBam.modules`|String|"samtools/1.16.1 picard/3.1.0"|The modules that will be loaded


### Outputs

Output | Type | Description | Labels
---|---|---|---
`downSampledFastq1`|File?|Output downsampled fastq read1|vidarr_label: downSampledFastq1
`downSampledFastq2`|File?|Output downsampled fastq read2|vidarr_label: downSampledFastq1
`downSampledBam`|File?|Downsampled bam file|vidarr_label: downSampledBam
`downSampledBai`|File?|Downsampled bam file index|vidarr_label: downSampledBai


## Commands

This section lists command(s) run by downSample workflow

* Running downSample

```
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
```
```
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
```

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_

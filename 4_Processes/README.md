# Processes

## Conditional execution

&nbsp;&nbsp;&nbsp;&nbsp;In order to conditionally run a script, add an `if` statement to your code:

```nextflow

params.compress = 'gzip'
params.file2compress = "$baseDir/data/ggal/transcriptome.fa"

process FOO {
    debug true

    input:
    path file

    script:
    if (params.compress == 'gzip')
        """
        echo "gzip -c $file > ${file}.gz"
        """
    else if (params.compress == 'bzip2')
        """
        echo "bzip2 -c $file > ${file}.bz2"
        """
    else
        throw new IllegalArgumentException("Unknown compressor $params.compress")
}

workflow {
    FOO(params.file2compress)
}
```

## When

&nbsp;&nbsp;&nbsp;&nbsp;`when` enables you to conditionally run the whole process. See the example below:

```nextflow
params.dbtype = 'nr'
params.prot = 'data/prots/*.tfa'
proteins = Channel.fromPath(params.prot)

process FIND {
    debug true

    input:
    path fasta
    val type

    when:
    fasta.name =~ /^BB11.*/ && type == 'nr'

    script:
    """
    echo blastp -query $fasta -db nr
    """
}

workflow {
    result = FIND(proteins, params.dbtype)
}
```

## Directives

&nbsp;&nbsp;&nbsp;&nbsp;Directives allow you to optionally define the setting of this process. They are declared at the top of the process body and are used mostly for computing resources:

```nextflow
process FOO {
    cpus 2
    memory 1.GB
    container 'image/name'

    script:
    """
    echo your_command --this --that
    """
}
```

See [whole list of directives](https://www.nextflow.io/docs/latest/process.html#directives). 

## PublishDir directive

&nbsp;&nbsp;&nbsp;&nbsp;If you would like to save important files created during the execution of the process, you can add the `publishDir` directive. It will save all files that match the given pattern in the directory.

Example:

```nextflow
reads_ch = Channel.fromFilePairs('data/ggal/*_{1,2}.fq')

process FOO {
    publishDir "results", pattern: "*.bam"

    input:
    tuple val(sample_id), path(sample_id_paths)

    output:
    tuple val(sample_id), path("*.bam")
    tuple val(sample_id), path("*.bai")

    script:
    """
    echo your_command_here --sample $sample_id_paths > ${sample_id}.bam
    echo your_command_here --sample $sample_id_paths > ${sample_id}.bai
    """
}

workflow {
    FOO(reads_ch)
}
```
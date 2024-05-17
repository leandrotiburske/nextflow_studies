# Welcome to my      <img src="https://upload.wikimedia.org/wikipedia/commons/thumb/e/e1/Logo_Nextflow_%28new%29.png/800px-Logo_Nextflow_%28new%29.png" width="125">     studies repository!


## Introduction

&nbsp;&nbsp;&nbsp;&nbsp;Nextflow is a free open-source programming language designed to create reproducible and scalable workflows! It is an extension of the Groovy programming language and is vastly used to orchestrate [bioinformatics](https://www.youtube.com/watch?v=W-Ov2cUaYQY) pipelines. 

> [!WARNING]
> &nbsp;&nbsp;&nbsp;&nbsp;The materials in this repository are utilized for my personal studies, derived from official [Nextflow tutorials](https://training.nextflow.io/basic_training/), and should be considered as additional resources, not substitutes.

## <img src="https://raw.githubusercontent.com/nf-core/logos/master/nf-core-logos/nf-core-logo.png" width=150>

&nbsp;&nbsp;&nbsp;&nbsp;[nf-core](https://nf-co.re/) is an incredible community project that aims to create a set of pipelines created with Nextflow!

> [!NOTE]
> Don't forget to take a look at the [Hackathons](https://nf-co.re/events/hackathon) page!

## Tips

<img src="https://www.freepnglogos.com/uploads/youtube-play-red-logo-png-transparent-background-6.png" width="20"> &nbsp;&nbsp;Watch the Nextflow & nf-core tutorials on YouTube (available in [English](https://www.youtube.com/watch?v=nzR8DRq13nE&list=PL3xpfTVZLcNhoWxHR0CS-7xzu5eRT8uHo&pp=iAQB), [Portuguese](https://www.youtube.com/watch?v=751E-yOH7H8&list=PL3xpfTVZLcNhi41yDYhyHitUhIcUHIbJg&pp=iAQB), [French](https://www.youtube.com/watch?v=dk7I5SZgAEU&list=PL3xpfTVZLcNhiv9SjhoA1EDOXj9nzIqdS&pp=iAQB), [Spanish](https://www.youtube.com/watch?v=pGUxK2WnwKs&list=PL3xpfTVZLcNhSlCWVoa3GURacuLWeFc8O&pp=iAQB) and [Hindi](https://www.youtube.com/watch?v=84XtbqRkKSk&list=PL3xpfTVZLcNikun1FrSvtXW8ic32TciTJ&pp=iAQB))

<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/8/84/Spotify_icon.svg/1982px-Spotify_icon.svg.png" width ="20"> &nbsp;&nbsp;Listen to  [CHANNELS - The Nextflow Podcast](https://www.nextflow.io/podcasts.html)

## Installation

&nbsp;&nbsp;&nbsp;&nbsp;In order to install Nextflow on your computer, run the following command in your terminal:


```
wget -qO- https://get.nextflow.io | bash
```

&nbsp;&nbsp;&nbsp;&nbsp;Alternative:

```
curl -s https://get.nextflow.io | bash
```

&nbsp;&nbsp;&nbsp;&nbsp;Finally, make sure it is executable:

```
chmod +x nextflow
```

> [!NOTE]
> &nbsp;&nbsp;&nbsp;&nbsp;Bash, Java 11 (or later, up to 21), Git and Docker are required

## Key concepts

&nbsp;&nbsp;&nbsp;&nbsp;Nextflow workflow orchestration is based on two very important concepts: **processes** and **channels**. A Nextflow pipeline is made of multiple processes joined by channels.  

### Processes

&nbsp;&nbsp;&nbsp;&nbsp;A **process** in Nextflow can be written in any programming language executable by Linux, such as Python, Shell Bash and Ruby. A process is isolated from other processes and are executed independently. A process example is the [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) algorithm, used to analyse DNA sequencing quality.

&nbsp;&nbsp;&nbsp;&nbsp;In this documentation, each process will be represented in diagrams by a hexagon node. A process representation is shown below:

```mermaid

flowchart LR;
    Process{{This is a process}};
```

### Channel

&nbsp;&nbsp;&nbsp;&nbsp;The only way different processes can communicate is through **channels**. In this documentation, each channel will be represented by a stadium-shaped node, as shown below:

```mermaid
flowchart LR;
    Channel([This is a Channel]);
```

&nbsp;&nbsp;&nbsp;&nbsp;Channels are asynchronous first-in, first-out ([FIFO](https://www.geeksforgeeks.org/fifo-vs-lifo-approach-in-programming/)) queues. Basically, every input and every output is considered a channel. 

```mermaid
flowchart LR
    Input(["Input\n(Channel)"]) --> Process{{Process}};
    Process{{Process}} --> Output(["Output\n(Channel)"])
```

&nbsp;&nbsp;&nbsp;&nbsp;Finally, a workflow is the execution of the pipeline made of the interaction of processes via channels. Considering the FastQC example, its input channel would be [FASTQ](https://knowledge.illumina.com/software/general/software-general-reference_material-list/000002211) files and its output channel would be output directories:

```mermaid
flowchart LR
subgraph workflow["DNA sequencing quality analysis workflow"]
    direction LR
    FASTQ(["FASTQ files"]) --> FastQC{{"FastQC"}};
    FastQC --> Outdir(["Output directories"])
end
```
# Channels

## Channel types

### Queue channel

&nbsp;&nbsp;&nbsp;&nbsp;Queue channels are:

- **asynchronous** (non-blocking operations),
- **unidirectional** (data flows from producer to consumer) and
- **FIFO** (first in, first out).

### Value channel

&nbsp;&nbsp;&nbsp;&nbsp;Value channels are bound to a single value. They are created through the `Channel.value()` channel factory or with operators that return only one value, like `collect`, `count` and `sum`.

## Channel factories

| Channel factory  | Function | Example |
| ------------- | ------------- | ------------- |
| `Channel.value()` | Create a value channel | `Channel.value('Hello!')` |
| `Channel.of()` | Create a queue channel with the values passed as arguments | `Channel.of(1, 3, 5, 7)` |
| `Channel.fromList()` | Create a channel with the elements of a list object passed as argument |` Channel.fromList(['hi', 'there'])` |
| `Channel.fromPath()` | Create a queue channel with one or more files matching the pattern used as argument | `Channel.fromPath('./data/meta/*.csv')` |
| `Channel.fromFilePairs` | Create a channel emitting file pairs that match the given pattern | `Channel.fromFilePairs('./FASTQ/*_{1,2}.fq')` |
| `Channel.fromSRA()` | Create a channel emitting FASTQ files from the SRA archive that match the given accession number(s) | `Channel.fromSRA(['ERR908507', 'ERR908506', 'ERR908505'])` |
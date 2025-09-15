BBTools are sensitive to filename extensions.  For example, this command:
reformat.sh in=reads.fq out=reads.fa.gz
...will convert reads from fastq format to gzipped fasta.  The recognized sequence file extensions are as follows:

fastq (fq)
fasta (fa, fna, fas, ffn, frn, seq, fsa, faa)
sam
bam [requires samtools]
qual
scarf [input only]
phylip [input only; only supported by phylip2fasta.sh]
header [output only]
oneline [tab delimited 2-column: name and bases]
embl [input only]
gbk [input only]

The recognized compression extensions:

gzip (gz) [can be accelerated by pigz]
zip
bz2 [requires bzip2 or pbzip2 or lbzip2]
fqz [requires fqz_comp]

In order to stream using standard in or standard out, it is recommended to include the format.  For example:
cat data.fq.gz | reformat.sh in=stdin.fq.gz out=stdout.fa > file.fa
This allows the tool to determine the format.  Otherwise it will revert to the default.

BBTools can usually determine the type of sequence data by examining the contents.  To test this, run:
fileformat.sh in=file

...which will print the way the data is detected, e.g. Sanger (ASCII-33) quality, interleaved, etc.  These can normally be overridden with the "qin" and "interleaved" flags.

When BBTools are processing gzipped files, they may, if possible, attempt to spawn a pigz process to accelerate it.  This behavior can be forced with the "pigz=t unpigz=t" flags, or prevented with "pigz=f unpigz=f"; otherwise, the default behavior depends on the tool.  In some cluster configurations, and some Amazon nodes, spawning a process may cause the program to killed with an indication that it used too much virtual memory.  I recommend pigz be enabled unless that scenario occurs.

The most recent extension added is "header".  You can use it like this:
reformat.sh in=reads.fq out=reads.header minlen=100

That will create a file containing headers of reads that pass the "minlen" filter.

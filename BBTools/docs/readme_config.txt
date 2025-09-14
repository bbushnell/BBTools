BBTools Config File Readme
Written by Brian Bushnell
Last updated May 12, 2015

A config file is a text file with a set of parameters that will be added to the command line.
The format is one parameter per line, with the # symbol indicating comments.
To use a config file, use the config=file flag.  For example, take BBDuk:

bbduk.sh in=reads.fq out=trimmed.fq ref=ref.fa k=23 mink=11 hdist=1 tbo tpe

That is equivalent to:

bbduk.sh in=reads.fq out=trimmed.fq ref=ref.fa config=trimadapters.txt
...if trimadapters.txt contained these lines:
k=23
mink=11
hdist=1
tbo
tpe


Any parameter placed AFTER the config file will override the same parameter if it is in the config file.  
For example, in this case k=20 will be used:
bbduk.sh in=reads.fq out=trimmed.fq ref=ref.fa config=trimadapters.txt k=20

But in this case, k=23 will be used, from the config file:
bbduk.sh in=reads.fq out=trimmed.fq ref=ref.fa k=20 config=trimadapters.txt

What are config files for?  Well, mainly, to overcome difficulties like whitespace in file paths, or command lines that are too long.
There are some example config files in bbmap/config/.  They are not used unless you specifically tell a program to use them.

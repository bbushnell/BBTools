BBMap/BBTools readme
Written by Brian Bushnell
Last updated January 2, 2018

The BBTools package is primarily devloped by Brian Bushnell, with some optional JNI and MPI components written by Jonathan Rood.
Some parts have also been written or modified by Shijie Yao, Alex Copeland, and Bryce Foster.


Citation:

Please see citation.txt


License:

The BBTools package is open source and free to use with no restrictions.  For more information, please read Legal.txt and license.txt.


Documentation:

Documentation is in the /bbmap/docs/ directory, and in each tool's shellscript in /bbmap/.
readme.txt: This file.
UsageGuide.txt: Contains basic installation and usage information.  Please read this first!
ToolDescriptions.txt: Contains a list of many BBTools, a description of what they do, and their hardware requirements.
compiling.txt: Information on compiling JNI code.
readme_config.txt: Usage information about config files.
readme_filetypes.txt: More detailed information on file formats supported by BBTools.
changelog.txt: List of changes by version, and current known issues.


Tool-specific Guides:

Some tools have specific guides, like BBDukGuide.txt.  They are in /bbmap/docs/guides/.  For complete documentation of a tool, I recommend that you read UsageGuide.txt first (which covers the shared functionality of all tools), then the tool's specific guide if it has one (such as ReformatGuide.txt), then the tool's shellscript (such as reformat.sh) which lists all of the flags.


Pipelines:

/bbmap/pipelines/ contains shellscripts.  These are different than the ones in /bbmap/, which are wrappers for specific tools.  The pipelines do not print a help message and do not accept any arguments.  They are given to provide examples of the command lines and order of tools used to accomplish specific tasks.


Resources:

/bbmap/resources/ contains various data files.  Most are fasta contaminant sequences.  For more information see /bbmap/resources/contents.txt.

If you have any questions not answered in the documentation, please look at the relevant SeqAnswers thread (linked from here: http://seqanswers.com/forums/showthread.php?t=41057) and post a question there if it is not already answered.  You can also contact JGI's BBTools team at bbtools@lbl.gov, or me at bbushnell@lbl.gov.  But please read the documentation first.

Special thanks for help with shellscripts goes to:
Alex Copeland (JGI), Douglas Jacobsen (JGI/NERSC), Bill Andreopoulos (JGI), sdriscoll (SeqAnswers), Jon Rood (JGI/NERSC), and Elmar Pruesse (UC Denver).

Special thanks for helping to support BBTools goes to Genomax (SeqAnswers).

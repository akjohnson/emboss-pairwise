pairwise -- An EMBOSS application
================================

*pairwise* is an EMBOSS application designed to do many-to-many pairwise alignments of sequences
using *needle*, in a way that can be parceled out for cluster computing.

WARNING!
--------

This code was written for EMBOSS 5.0.0.  It will need to be reviewed and verified before being
used with later versions--the current version of EMBOSS might require some changes.

Additionally, EMBOSS now has an application called "needleall" which also does pairwise needle
comparisons.  However, this implementation has the advantage of being able to specify start and
end sequences, so that the work can be done in chunks and easily spread out on a cluster computer.

Installing Pairwise
-------------------

In order to install pairwise in an EMBOSS installation, download the EMBOSS source code
and copy the "emboss" directory from this distribution into the main EMBOSS directory.

Then you should edit the file:

    emboss/Makefile.am

Add pairwise to the bin_PROGRAMS list:

	octanol oddcomp \
	palindrome pairwise pasteseq patmatdb patmatmotifs \	
	pepcoil pepinfo pepnet pepstats \

And this to the SOURCES section:

	oddcomp_SOURCES = oddcomp.c
	pairwise_SOURCES = pairwise.c
	palindrome_SOURCES = palindrome.c

More instructions on this process are [available on the EMBOSS Sourceforge website](http://emboss.sourceforge.net/developers/program.html).

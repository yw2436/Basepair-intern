DNemulator
==========

DNemulator is a package for simulating DNA sequencing errors,
polymorphisms, cytosine methylation and bisulfite conversion.

Installation
------------

Optionally, you can copy the programs to a standard "bin" directory::

  sudo make install

Or copy them to your personal ~/bin directory::

  make install prefix=~

fasta-polymorph
---------------

This program randomly modifies DNA sequences, using given
polymorphisms and their frequencies.  The polymorphisms should be
supplied in the UCSC Genome "common SNPs" format.  Run it like this::

  fasta-polymorph snp132Common.txt genome.fa > polymorphed.fa

To simulate a diploid genome, give it the same sequences twice::

  fasta-polymorph snp132Common.txt genome.fa genome.fa > diploid.fa

Since the polymorphisms may include insertions and deletions, we might
wish to record the alignments between original and polymorphed
sequences.  The "-a" option writes the alignments to a file::

  fasta-polymorph -a align.seg snp132Common.txt genome.fa > polymorphed.fa

The alignments are in SEG format, which is described in the seg-suite
package.

fasta-random-chunks
-------------------

This program gets random fixed-length chunks of the input sequences.
Run it like this::

  fasta-random-chunks sequences.fa > chunks.fa

Chunks containing uppercase letters other than ACGT are suppressed.
However, lowercase letters are not suppressed.  (This enables us to
abuse lowercase letters to represent things like methylation rates.)

The title line of each chunk receives the following data::

  > serialNumber chunkLength sourceName zeroBasedStartCoordinate strand

Therefore, it is possible to get the alignments of the chunks to the
input sequences, in SEG format, like this::

  awk '/>/ {print $3, $4, $5, $2, $6 == "+" ? 0 : -$3}' chunks.fa > chunks.seg

The program has these options:

-f    Get forward-strand chunks only.
-n N  Get N chunks (default 100000).
-s S  Use a chunk length of S (default 36).
-a A  Simulate unequal sequence abundances, using a Pareto (Type I)
      distribution with Pareto index A.  For typical RNA abundances,
      try "-a 1.25".

fastq-sim
---------

This program simulates DNA sequencing errors.  It joins sequences from
a fasta file to quality symbols from a fastq file, then randomly
mutates the sequences according to the error probabilities represented
by the quality symbols.  Run it like this::

  fastq-sim sequences.fasta qualities.fastq > simulated.fastq

If the number of fasta sequences differs from the number of fastq
sequences, the output will have the minimum of the two.  Likewise, if
the length of a fasta sequence differs from the length of the
corresponding fastq sequence, the output will be the shorter of the
two.

The fastq data is assumed to be in fastq-sanger format, with no
line-wrapping or blank lines.

fasta-paired-chunks
-------------------

This program gets random *pairs* of chunks of the input sequences.
Each pair comes from the ends of a randomly-chosen fragment: the first
chunk is from the start of the fragment's forward strand, and the
second chunk is from the start of the fragment's reverse strand.  Run
it like this::

  fasta-paired-chunks sequences.fa output1.fa output2.fa

The fragment lengths are randomly chosen from a Normal Distribution.
Chunks containing nonstandard letters are suppressed, and chunk titles
are made, in the same way as fasta-random-chunks.

If the output file names are the same, the chunks will be written
interleaved.  The special file name "-" refers to stdin or stdout, as
appropriate.

The program has these options:

-n NUM, --num=NUM    Get NUM pairs of chunks (default 100000).
-l BP, --length=BP   Use a chunk length of BP base-pairs (default 50).
                     You can specify two comma-separated lengths.
-c, --cap            Get fragments from starts of the input sequences
                     only.
-d, --directional    Get fragments from forward strands of the input
                     sequences only.
-f BP, --fraglen=BP  Use a mean fragment length of BP base-pairs
                     (default 250).
-s BP, --sdev=BP     Use a fragment length standard deviation of BP
                     base-pairs (default fraglen/5).

fasta-methyl-sim
----------------

This program reads DNA sequences and randomly assigns a methylation
rate to every cytosine.  Run it like this::

  fasta-methyl-sim sequences.fa > meth.fa

It replaces Cs with symbols representing different methylation rates.
It also replaces Gs with symbols representing methylation rates of the
reverse strand.  The sequences are first converted to uppercase, so
that lowercase letters can be used for methylation rate symbols.

========  ========  ===========  ==============  ==============
C symbol  G symbol  Methylation  Probability in  Probability in
                    rate         cg context      non-cg context
========  ========  ===========  ==============  ==============
C         G         ~0%          0.1             0.96          
c         g         ~10%         0.1             0.01          
d         h         ~20%         0.1             0.01          
v         b         ~50%         0.1             0.01          
t         a         ~100%        0.6             0.01          
========  ========  ===========  ==============  ==============

These symbols abuse the standard code for (ambiguous) bases, in such a
way that complementing them works properly (e.g. v <=> b).

The symbols are interpreted by fasta-bisulf-sim.

fasta-bisulf-sim
----------------

This program simulates bisulfite conversion of cytosines with given
methylation rates.  It replaces the C symbols produced by
fasta-methyl-sim with either "t" (converted) or "C" (not converted).
It replaces the G symbols with "G".  Run it like this::

  fasta-bisulf-sim meth.fa > converted.fa

The conversion rates are:

========  ===============
C symbol  Conversion rate
========  ===============
C         0.99           
c         0.9            
d         0.8            
v         0.5            
t         0              
========  ===============

To simulate conversion of the reverse strand::

  fasta-bisulf-sim -r meth.fa > converted.fa

Example
-------

This example simulates bisulfite-converted DNA reads.  First, we
simulate methylation rates for every cytosine in a genome::

  fasta-methyl-sim genome/chr*.fa > meth.fa

Then, we simulate a polymorphic, diploid genome.  The alignments
between simulated and original chromosomes are written to a file
called align.seg::

  fasta-polymorph -a align.seg snp132Common.txt meth.fa meth.fa > poly.fa

Finally, we get 1 million length-85 chunks, simulate bisulfite
conversion, and simulate sequencing errors.  The error probabilities
are supplied by a file called sample.fastq::

  fasta-random-chunks -n1000000 -s85 poly.fa |
  fasta-bisulf-sim |
  fastq-sim - sample.fastq > reads.fastq

Miscellaneous
-------------

DNemulator is distributed under the GNU General Public License, either
version 3 of the License, or (at your option) any later version.  For
details, see COPYING.txt.

Please cite: "A mostly traditional approach improves alignment of
bisulfite-converted DNA".  MC Frith, R Mori, K Asai.  Nucleic Acids
Research (2012) 40 (13): e100.

Website: http://www.cbrc.jp/dnemulator/

Email (questions, comments, problems):
dnemulator (ATmark) cbrc (dot) jp.

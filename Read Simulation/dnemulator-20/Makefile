# ugh, got to keep these up to date:

progs = fasta-bisulf-sim fasta-methyl-sim fasta-paired-chunks	\
fasta-polymorph fasta-random-chunks fastq-sim

texts = COPYING.txt README.txt ChangeLog.txt

distdir = dnemulator-`hg id -n`

README.html: README.txt
	rst2html README.txt > README.html

log:
	hg log --style changelog > ChangeLog.txt

dist: README.html log
	mkdir ${distdir}
	cp ${progs} Makefile ${texts} README.html ${distdir}
	zip -qrm ${distdir} ${distdir}

prefix = /usr/local
exec_prefix = ${prefix}
bindir = ${exec_prefix}/bin
install:
	mkdir -p ${bindir}
	cp ${progs} ${bindir}

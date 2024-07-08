#!/usr/bin/sh

# to show how host-parasite pair genome sequences files 
# were concatenated to prevent >100% mapping

cat human.fa Pfalciparum.fa > humanPfalciparum.fa
cat human.fa Pberghei.fa > humanPberghei.fa
cat human.fa Pvivax.fa > humanPvivax.fa

cat mouse.fa Pchabaudi.fa > mousePchabaudi.fa
cat mouse.fa Pyoelii.fa > mousePyoelii.fa
cat mouse.fa Pberghei.fa > mousePberghei.fa

cat monkey.fa Pcoatneyi.fa > monkeyPcoatneyi.fa
cat monkey.fa Pcynomolgi.fa > monkeyPcynomolgi.fa
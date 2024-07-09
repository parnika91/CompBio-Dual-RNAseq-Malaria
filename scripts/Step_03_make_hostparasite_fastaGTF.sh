#!/usr/bin/sh

# to show how host-parasite pair genome sequences files 
# were concatenated to prevent >100% mapping

cat human.gtf Pfalciparum.gtf > humanPfalciparum.gtf
cat human.gtf Pberghei.gtf > humanPberghei.gtf
cat human.gtf Pvivax.gtf > humanPvivax.gtf

cat mouse.gtf Pchabaudi.gtf > mousePchabaudi.gtf
cat mouse.gtf Pyoelii.gtf > mousePyoelii.gtf
cat mouse.gtf Pberghei.gtf > mousePberghei.gtf

cat monkey.gtf Pcoatneyi.gtf > monkeyPcoatneyi.gtf
cat monkey.gtf Pcynomolgi.gtf > monkeyPcynomolgi.gtf

# same for gtf files

cat human.gtf Pfalciparum.gtf > humanPfalciparum.gtf
cat human.gtf Pberghei.gtf > humanPberghei.gtf
cat human.gtf Pvivax.gtf > humanPvivax.gtf

cat mouse.gtf Pchabaudi.gtf > mousePchabaudi.gtf
cat mouse.gtf Pyoelii.gtf > mousePyoelii.gtf
cat mouse.gtf Pberghei.gtf > mousePberghei.gtf

cat monkey.gtf Pcoatneyi.gtf > monkeyPcoatneyi.gtf
cat monkey.gtf Pcynomolgi.gtf > monkeyPcynomolgi.gtf

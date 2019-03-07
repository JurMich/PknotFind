


PPPPP KKK  KK                 t   FFFFF              ddd  
 P   P K  K                   t   F    i              d  
 P   P K K                    t   F                   d 
 PPPP  KK      n nnn    ooo  tttt FFFF i  n nnn    dddd
 P     K K     nn   n  o   o  t   F    i  nn   n  d   d
 P     K  K    n    n  o   o  t   F    i  n    n  d   d
 P     K   K   n    n  o   o  t   F    i  n    n  d   d
PPP   KKK  KK nnn  nnn  ooo    tt F    i nnn  nnn  ddddd


PknotFind is a software for sanple pseudoknots non-redundantly.

WARNING!

This software contains some files ov ViennaRNA package:

https://www.tbi.univie.ac.at/RNA/

Citation:  Lorenz, Ronny and Bernhart, Stephan H. and Höner zu Siederdissen, 
Christian and Tafer, Hakim and Flamm, Christoph and Stadler, Peter F. and 
Hofacker, Ivo L.

ViennaRNA Package 2.0

Algorithms for Molecular Biology, 6:1 26, 2011, doi:10.1186/1748-7188-6-26 

The license the applies on these files is the same as the one applying on 
ViennaRNA library (see ViennaRNA_license.txt).

On the rest of files, the following applies:

You can freely use, modify and distribute them as long as the original source 
is mentioned and the software is redistributed without any fee. The software is 
distributed WITHOUT ANY WARRANTY, without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

Version 1.2:

verhaul of adjacent algorithm, now for each (i,j) there is a limited list of
helices up to certain number (default 11, can be set with -l). This reduces
the complexity to n^4 (time), resp. n^2 (memory). Also added option to search
only the pseudoknots above given length (-m).

------------------------------------------------------------------------
|>|>|> INSTALL:

WARNING:
Before the installation, make sure you have the last version of ViennaRNA
library installed!

- clone directory
- cd Pknotfind
- make

|>|>|> Running:
- ./PknotFind -i name_of_input_file -p num_samples

The default number of samples is 100.

Other options:
-j: skew for helices (default 2)
-k: energy penalty (default 900)
-l: length of list+1 (default 11)
-m: minimum pk length (default 6)
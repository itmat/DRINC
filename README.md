# DRINC
Duplication Rate Identifier for NGS Cleanup.


Usage: python DRINK5Cluster.py <stranded info> <radius> <input sam file> <output all dups log file> <output optical dups log file> <output stats file> [filtering mode] [output sam file w/ tags]


For the stranded info,
 - enter 0 for a nonstranded experiment.
 - enter 1 for a stranded experiment.


The radius is the max possible distance between two optical duplicates.
 - If you set this too high you might call legitimate dups as optical dups.
 - If you set it too low you might miss optical dups.


The filtering mode is optional.
 - Enter 0 for just outputing new sam with tags, without filtering.
 - Enter 1 for filtering all dups.
 - Enter 2 for filtering optical dups only.

The output sam file w/ tags field is to be specified if and only if the filtering mode is specified.


#### Written by Peng Song.

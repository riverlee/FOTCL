# Feature Occupancy Transition between Cell Lines (FOTCL)

## Description

Given a list of genome intervals of interest (e.g, bivalent domains which are expressed in ESC but not expressed in HMEC), we would like to see which TFs/DNase/Histone Marks gain/loss the occupancy between ESC and HMEC. 

 
## Example

In the example folder, Esc.expressed_To_Hmec.notexpressed.bed store genome intervals of bivalent domains both in ESC and HMEC, however, they are expressed in ESC but not in HMEC.

```
perl ../FOTCL.pl -i Esc.expressed_To_Hmec.notexpressed.bed -o Esc.expressed_To_Hmec.notexpressed_FOTCL.txt -a H1-HESC -b HMEC
```



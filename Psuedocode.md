# Deduper Psuedocode

# Define the Problem 
I need to remove all PCR duplicates from single-end data using 96 UMIs.


# Psuedocode
1. Read STL96 file to parse UMIs.
2. Place UMIs into set. 
3. Read input sam file by line, skipping headers. 
4. Determine what chromosome read is on. Save as global variable.
5. Create a global set that will hold info for first entries (non-duplicates).
6. Determine if UMI is in UMI set. If it is, add to a tuple. If not, it will not be written to output file.
7. Read bitwise flag (SAM col 2 - bit 16) and determine which strand the read is on. (add to tuple)
8. Determine start position of read (POS - SAM col 4). (add to tuple)
9. Parse CIGAR string (col 6) to figure out if soft clipping has occured. Use info from CIGAR string to adjust position of 5' start of read. (adjust tuple)
10. Evaluate if current tuple is in set. If not, add tuple to set and write to saved output file.
11. Read next line. If chromosome is different than global variable, clear set and continue steps 6-10. 




# High-level functions


parse_UMIs(str: file):
```This function will read through my STL96.txt file and return a set of the UMIs.```
return UMI_set

#parse_UMIs example:

parse_UMIs(STL96.txt)
output: UMI_set

print(UMI_dict)
output: (AACGCCAT,AAGGTACG,AATTCCGG.....TTCGCCTA,TTCGTTCG)

##

bitwise_strand(str: bitwise flag):
```This function will evaluate the bitwise flag and determine which strand the read is on.```
return strand

#bitwise_strand example:

bitwise_strand(83)
output: -

bitwise_strand(203)
output: +

bitwise_strand(159)
output: -


##

cigar_pos(str: cigar string, str: 5prime_position):
```This funuction will evaluate the cigar string and determine a numerical value to adjust the position of the 5' end```
return 5_end_pos

#cigar_pos example

cigar_pos(14M, 14729)
output: 14729 

cigar_pos(2S12M, 14729)
output: 14727

cigar_pos(12M2S, 14729)
output: 14729








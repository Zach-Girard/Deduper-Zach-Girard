#!/usr/bin/env python

# import regex package
import re


import argparse

def get_args():
    parser = argparse.ArgumentParser(description="A program used to determine remove PCR duplicates in a sorted SAM file and return the non-duplicate lines.")
    parser.add_argument("-f", "--file", help="designates absolute file path to sorted sam file.", required=True, type=str)
    parser.add_argument("-o", "--outfile", help="designates absolute file path to sorted sam file.", required=True, type=str)
    parser.add_argument("-u", "--umi", help="designates file containing the list of UMIs.", required=True, type=str)
    parser.print_help(file=None)
    return parser.parse_args()
        
args = get_args()


# Set empty set to hold known UMIs
UMI_set = set()
# Retreive UMIs
with open(args.umi, "r") as file:
    for line in file:
        UMI = line.strip()
        UMI_set.add(UMI)
file.close()


# Function to parse CIGAR string
# This function returns a list of tuples
# Example return: [(5,"S"),(84,"M")]
def CIGAR_smoke(CIGAR):
    parsed_list = []
    for match in re.finditer(r"(\d+)([MIDNSHP=X])", CIGAR):
        number = int(match.group(1))
        letter = match.group(2)
        parsed_list.append((number, letter))
    return parsed_list



# To keep track of unknown UMIs
unknown_UMIs = 0

# Create global value to keep track of non-duplicate info
Non_duplicates = set()

# To keep track of PCR duplicates
Duplicates = 0

# To keep track of header lines
headers = 0

# Total non-duplicates written to output file
lines_written = 0

# Define Chromosome to clean Non-duplicate set and save memory
Chrom = 1

# Go through SAM file and place non-duplicates in output file
with open(args.file, "r") as input, open(args.outfile, "w") as output:
    for line in input:
        if line == "":
            break

        # This writes the header lines to the output file
        if line.startswith("@"):
            headers += 1
            output.write(line)
        
        # Parsing the read lines
        else:
            splitline = line.strip().split()
            # Chromosome of current read
            Current_Chrom = splitline[2]
            # Clear set if new Chromosome
            if Current_Chrom != Chrom:
                Non_duplicates.clear()
                Chrom = Current_Chrom

            # Q Name that contains UMI
            QNAME = splitline[0]
            QNAME = QNAME.split(":")
            # Pulling UMI from Q Name
            UMI = QNAME[7]
            # Pulling left most position of the read
            position = int(splitline[3])
            # pulling the CIGAR string
            CIGAR = splitline[5]
            
            # Checking if there is an unknown UMI
            if UMI not in UMI_set:
                unknown_UMIs += 1

            # If UMI is known, continue
            else:
                # Pulling bitwise flag
                flag = int(splitline[1])
                # Checking for strandedness
                if ((flag & 16) == 16):
                    strand = "minus"
                    # Call CIGAR parsing function
                    parsed_list = CIGAR_smoke(CIGAR)
                    # If there is not soft clipping on the left side, normally adjust position
                    if parsed_list[0][1] != "S":
                        for i in range(len(parsed_list)):
                            if parsed_list[i][1] == "M" or "N" or "D" or "S":
                                position = position + parsed_list[i][0]
                    # If there is soft clipping on left side, skip it
                    else:
                        for i in range(1, len(parsed_list)):
                            if parsed_list[i][1] == "M" or "N" or "D" or "S":
                                position = position + parsed_list[i][0]

    
                # If positive, check for soft clipping
                else:
                    strand = "plus"
                    parsed_list = CIGAR_smoke(CIGAR)
                    if parsed_list[0][1] == "S":
                        position = position - parsed_list[0][0]
                
                # If the unique read info has not been seen, write to file, add to count
                if (position, strand, UMI) not in Non_duplicates:
                    Non_duplicates.add((position, strand, UMI)) # type: ignore
                    output.write(line)
                    lines_written += 1

                # Counts duplicates 
                else:
                    Duplicates += 1



output.close()

# Print stats
print(f"Number of Duplicates: {Duplicates}")
print(f"Number of unknown UMIs: {unknown_UMIs}")
print(f"Number of headers: {headers}")
print(f"Number of Non-duplicate lines: {lines_written}")

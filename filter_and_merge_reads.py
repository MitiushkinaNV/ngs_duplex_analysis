import pysam
import sys
import os
import re
import argparse
import statistics

parser=argparse.ArgumentParser()
parser.add_argument("-f","--file_in",help="Input bam file name", required=True)
parser.add_argument("-o","--file_out",help="Output sam file name", required=True)
parser.add_argument("-c","--overcovered_regions",help="File with overcovered regions")

args=parser.parse_args()

raw_bam_file_name=args.file_in
collapsed_sam_file_name=args.file_out

overcovered_starts_dict={}
overcovered_ends_dict={}
if args.overcovered_regions:
    overcovered_file=open(args.overcovered_regions,'rt')
    for line in overcovered_file:
        overcovered_chr,overcovered_start,overcovered_end=line.split()
        if overcovered_chr in overcovered_starts_dict:
            overcovered_starts_dict[overcovered_chr].append(int(overcovered_start))
            overcovered_ends_dict[overcovered_chr].append(int(overcovered_end))
        else:
            overcovered_starts_dict[overcovered_chr]=[int(overcovered_start),]
            overcovered_ends_dict[overcovered_chr]=[int(overcovered_end),]

samfile=pysam.AlignmentFile(raw_bam_file_name, "rb")
collapsed_samfile=pysam.AlignmentFile(collapsed_sam_file_name, "wb", template=samfile)
current_read_name=""
reads_with_same_name=list()

tlen_vect = []
def process_reads_with_same_name(reads_with_same_name, overcovered_starts_dict, overcovered_ends_dict):
    read1=""
    read2=""
    global tlen_vect
    
    for same_read in reads_with_same_name:
        #check if this is a primary read
        if not same_read.is_secondary:
            if same_read.is_read1:
                if not read1:
                    read1=same_read
                else:
                    #this read is probably split aligned
                    return("Split-aligned R1 or R2")
            elif same_read.is_read2:
                if not read2:
                    read2=same_read
                else:
                    return("Split-aligned R1 or R2")
            else:
                print("Strange: this is not read1 or read2\n")
                print(same_read.tostring())

    if not(read1 and read2):
        return("No R1 or R2")
    elif (read1.mapping_quality<60 or read2.mapping_quality<60):
        if (abs(read1.template_length)==abs(read2.template_length)):
            tlen_vect.append(abs(read1.template_length))
        return("R1 or R2 mapping quality lower than 60")
    elif (re.search('S',read1.cigarstring) or re.search('S',read2.cigarstring)):
        return("Soft-clipped bases in R1 or R2")
    elif (read1.is_reverse and read2.is_reverse):
        return("Both reads have reverse orientation")
    elif ((not read1.is_reverse) and (not read2.is_reverse)):
        return("Both reads have forward orientation")
    elif (read1.template_length>1000 or read1.template_length<-1000):
        return("Template length longer than 1000 bp")
    elif (read1.template_length==0 or read2.template_length==0):
        return("Template length 0")
    elif (not read1.is_proper_pair or not read2.is_proper_pair):
        return("R1 and R2 are not a proper pair")
    else:
        if str(read1.reference_name)==str(read2.reference_name):
            chrom=str(read1.reference_name)
            if (abs(read1.template_length)==abs(read2.template_length)):
                tlen_vect.append(abs(read1.template_length))
        else:
            return("R1 and R2 are not a proper pair")
        if read1.reference_start and read1.reference_end and read2.reference_start and read2.reference_end:
            positions_list=[read1.reference_start,read1.reference_end,read2.reference_start,read2.reference_end]
        else:
            return("R1 and R2 are not a proper pair")
        pos_start=min(positions_list)
        pos_end=max(positions_list)
        in_overcovered=0
        if chrom in overcovered_starts_dict:
            for i in range(0,len(overcovered_starts_dict[chrom])):
                if pos_start>overcovered_ends_dict[chrom][i]:
                    i+=1
                    continue
                elif pos_end<overcovered_starts_dict[chrom][i]:
                    break
                elif pos_start<overcovered_starts_dict[chrom][i] and \
                    pos_end>=overcovered_starts_dict[chrom][i]:
                    in_overcovered=1
                    break
                elif pos_start>=overcovered_starts_dict[chrom][i] and \
                    pos_start<=overcovered_ends_dict[chrom][i]:
                    in_overcovered=1
                    break
        if in_overcovered:
            return("Within overcovered region")
        elif not read1.is_reverse:
            read1.tags += [('rr', '"'+read2.tostring()+'"')]
            read1.tags += [('re', read2.reference_end)]
            collapsed_samfile.write(read1)
            return("Reads passed")
        else:
            read2.tags += [('rr', '"'+read1.tostring()+'"')]
            read2.tags += [('re', read1.reference_end)]
            collapsed_samfile.write(read2)
            return("Reads passed")

number_of_reads=0

read_dict={"Split-aligned R1 or R2":0,
           "No R1 or R2":0,
           "R1 or R2 mapping quality lower than 60":0,
           "Soft-clipped bases in R1 or R2":0,
           "Both reads have reverse orientation":0,
           "Both reads have forward orientation":0,
           "Template length longer than 1000 bp":0,
           "Template length 0":0,
           "R1 and R2 are not a proper pair":0,
           "Within overcovered region":0,
           "Reads passed":0}

for read in samfile.fetch(until_eof=True):
    if not current_read_name:
        current_read_name=read.query_name
        continue
    else:
        if read.query_name==current_read_name:
            reads_with_same_name.append(read)
            continue
        else:
            number_of_reads+=1
            res=process_reads_with_same_name(reads_with_same_name, overcovered_starts_dict, overcovered_ends_dict)
            read_dict[res]+=1
            current_read_name=read.query_name
            reads_with_same_name=[read,]
else:
    number_of_reads+=1
    res=process_reads_with_same_name(reads_with_same_name, overcovered_starts_dict, overcovered_ends_dict)
    read_dict[res]+=1
            
print("Overall number of reads or clusters in bam file")
print(number_of_reads)

for each_reason in read_dict.keys():
    print(each_reason,": ",read_dict[each_reason])

print("Median template length: ",statistics.median(tlen_vect))
samfile.close()
collapsed_samfile.close()

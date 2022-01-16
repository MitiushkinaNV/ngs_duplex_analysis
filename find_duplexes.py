import pysam
import sys
import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-f","--file_in",help="Input file name", required=True)
parser.add_argument("-o","--file_out_prefix",help="Prefix for R1 and R2 sam files", required=True)
parser.add_argument("-n","--read_num",help="Min number of reads need to build a consensus", type=int, default=5)
parser.add_argument("-l","--length_precision",\
        help="Asseptable difference between read lengths in a duplex (plus,minus), bp",type=int,default=2)

args=parser.parse_args()

samfile=pysam.AlignmentFile(args.file_in, "rb")
text_header=samfile.header

read1_file=open(args.file_out_prefix+'R1_temp_file.sam','wt')
print(text_header, file=read1_file, end='')

read2_file=open(args.file_out_prefix+'R2_temp_file.sam','wt')
print(text_header, file=read2_file, end='')

if not os.path.exists('duplex_statistics.txt'):
    duplex_statistics_file=open('duplex_statistics.txt','at')
    print("sample","chrom","pos_start","tlen","reads_f","reads_r", file=duplex_statistics_file)
else:
    duplex_statistics_file=open('duplex_statistics.txt','at')

n_reads_dropped=0
n_read_groups=0

def process_same_reads(same_reads,tlen_dict):
    #define tlen with max number of reads
    global n_reads_dropped
    global n_read_groups
    n_read_groups+=1
    max_reads_with_tlen=0
    current_tlen=0
    for tlen_key in tlen_dict.keys():
        if tlen_dict[tlen_key]>max_reads_with_tlen:
            max_reads_with_tlen=tlen_dict[tlen_key]
            current_tlen=tlen_key
    
    reads_f=[]
    reads_r=[]
    for same_read in same_reads:
        if same_read.is_read1 and abs(abs(same_read.template_length)-current_tlen)<=args.length_precision:
            reads_f.append(same_read)
        elif same_read.is_read2 and abs(abs(same_read.template_length)-current_tlen)<=args.length_precision:
            reads_r.append(same_read)
        else:
            n_reads_dropped+=1
            continue
    print(args.file_out_prefix,str(same_reads[0].reference_name),str(same_reads[0].reference_start),\
          str(current_tlen),str(len(reads_f)),str(len(reads_r)),file=duplex_statistics_file)
    if len(reads_f)>=args.read_num and len(reads_r)>=args.read_num:
        #we have a duplex
        for read_f in reads_f:
            rf,rr=read_f.tostring().split('\trr:Z:"')
            rr=(rr.split('"\tre:i:'))[0]
            read1_file.write(rf+"\n")
            read1_file.write(rr+"\n")
        for read_r in reads_r:
            rf,rr=read_r.tostring().split('\trr:Z:"')
            rr=(rr.split('"\tre:i:'))[0]
            read2_file.write(rf+"\n")
            read2_file.write(rr+"\n")
        return(1) #this means that we have found a duplex
    else:
        return(0) #no duplex

current_chr=""
current_pos=0
current_ref_end=0
tlen_dict={} #number of reads per tlen value
same_reads=list()

duplex_chr=[]
duplex_start=[]
duplex_end=[]
total_reads=0

for read in samfile.fetch(until_eof=True):
    total_reads+=1
    if not current_chr:
        #start read collection
        #print("start read collection")
        current_chr=str(read.reference_name)
        current_pos=read.reference_start
        current_ref_end=read.get_tag("re")
        tlen_dict[abs(read.template_length)]=1
        same_reads=[read,]
        continue
    elif str(read.reference_name)==current_chr and \
        (read.reference_start-current_pos)<=args.length_precision:
        #print("collect read")
        same_reads.append(read)
        if abs(read.template_length) in tlen_dict:
            tlen_dict[abs(read.template_length)]+=1
        else:
            tlen_dict[abs(read.template_length)]=1
        if read.get_tag("re")>current_ref_end:
            current_ref_end=read.get_tag("re")
    elif str(read.reference_name)==current_chr and \
        (read.reference_start-current_pos)>args.length_precision:
        if len(same_reads)>0:
            if process_same_reads(same_reads, tlen_dict):    
                duplex_chr.append(current_chr)
                duplex_start.append(current_pos)
                duplex_end.append(current_ref_end)
        if len(duplex_chr)==0 or \
           (str(read.reference_name)==duplex_chr[-1] and read.reference_start>duplex_end[-1]) or \
            str(read.reference_name)!=duplex_chr[-1]:
            #initiate collection of a new group of reads
            current_chr=str(read.reference_name)
            current_pos=read.reference_start
            current_ref_end=read.get_tag("re")
            tlen_dict={}
            tlen_dict[abs(read.template_length)]=1
            same_reads=[read,]
            #print("initiate new collection")
        else:
            #drop read
            n_reads_dropped+=1
            tlen_dict={}
            same_reads=[]
    else:
        #new chr
        if len(same_reads)>0:
            if process_same_reads(same_reads, tlen_dict):    
                duplex_chr.append(current_chr)
                duplex_start.append(current_pos)
                duplex_end.append(current_ref_end)
        current_chr=str(read.reference_name)
        current_pos=read.reference_start
        current_ref_end=read.get_tag("re")
        tlen_dict={}
        tlen_dict[abs(read.template_length)]=1
        same_reads=[read,]
        #print("new chr")
else:
    if len(same_reads)>0:
        if process_same_reads(same_reads, tlen_dict):    
            duplex_chr.append(current_chr)
            duplex_start.append(current_pos)
            duplex_end.append(current_ref_end)

duplex_positions_file=open(args.file_out_prefix+'_dup_positions_temp_file.txt','wt')
print("We found ",len(duplex_chr)," number of duplexes")
print("among", n_read_groups," read groups with same coordinates")
print("Total number of reads was", total_reads)
print(n_reads_dropped," reads were dropped because they intersected with already found duplexes")
for chrom,pos_start,pos_end in zip(duplex_chr,duplex_start,duplex_end):
    print(chrom,pos_start,pos_end, file=duplex_positions_file)
duplex_positions_file.close()

read1_file.close()
read2_file.close()
samfile.close()
duplex_statistics_file.close()
    

        
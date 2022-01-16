import pysam
import re
import sys
import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument("-f1","--file_in_r1",help="Input bam file for R1", required=True)
parser.add_argument("-f2","--file_in_r2",help="Input bam file for R2", required=True)
parser.add_argument("-d","--dup_positions",help="File with the positions of duplexes", required=True)
parser.add_argument("-r","--reference_file",help="Reference fasta file", required=True)
parser.add_argument("-v","--vcf_file",help="Vcf file with SNP data", required=True)
parser.add_argument("-l","--sample_label",help="Sample label", required=True)
parser.add_argument("-c","--chr_dict",\
                    help="File with chromosome names, which has to be processed, with their names in BAM and VCF files", \
                    required=True)
parser.add_argument("-cn","--cycle_number",help="Number of cycles per single read in paired-end sequencing",type=int,\
            required=True)
parser.add_argument("-n","--read_num",help="Min number of reads to build a consensus", type=int, default=5)
parser.add_argument("-f","--fraction",\
            help="Min fraction of reads with the same nucleotide required to include this nucleotide in a consensus",\
            type=float, default=0.8)
parser.add_argument("-s","--skip_from_end",help="Number of nucleotides to skip from each end of the read",\
            type=int, default=5)
parser.add_argument("-p","--position_statistics_file",help="Name for the output file with main results",\
            default="position_statistics.txt")
parser.add_argument("-af","--max_af",help="Polymorphisms at this population frequency will be excluded", type=float,\
            default=0.01)
parser.add_argument("-m","--mutation_limit",help="Reads with this or higher number of mutations will be exluded",\
                    type=int, default=3)


args=parser.parse_args()

read1_samfile=pysam.AlignmentFile(args.file_in_r1,'rb')
read2_samfile=pysam.AlignmentFile(args.file_in_r2,'rb')
duplex_positions_file=open(args.dup_positions,'rt')
ref_seq=pysam.FastaFile(filename=args.reference_file)
vcf_file = pysam.VariantFile(args.vcf_file)

#print(vcf_file.header)
#vcf_file.close()
#exit()

if not os.path.exists(args.position_statistics_file):
    position_statistics_file=open(args.position_statistics_file,'at')
    print("\t".join(["sample","chrom","mut_position","ref","alt","dup_um_snp","read_start",\
          "read_end","rs","AF","Ndup_read","Nsnp_read","Num_read","r1_count","r2_count"]),file=position_statistics_file)
else:
    position_statistics_file=open(args.position_statistics_file,'at')

vcf_chr_dict={}
chr_dict_file=open(args.chr_dict,'rt')
for line in chr_dict_file:
    bam_chr,vcf_chr=line.strip('\n').split(' ')
    vcf_chr_dict[bam_chr]=vcf_chr
chr_dict_file.close()

def pileup_reads(duplex_pileup):
    final_read={}
    n_reads=0
    for pileupcolumn in duplex_pileup:
        col_bases=pileupcolumn.get_query_sequences(mark_matches=True, add_indels=True)
        if len(col_bases)<args.read_num:
            continue
        n_reads=max(n_reads,len(col_bases))
        base_dict={}
        for each_base in col_bases:
            each_base=each_base.upper()
            if each_base in base_dict:
                base_dict[each_base]+=1
            else:
                base_dict[each_base]=1
        highest_freq=max(base_dict.values())
        if highest_freq/len(col_bases)<args.fraction:
            final_read[pileupcolumn.pos]='N'
        else:
            for each_key in base_dict.keys():
                if base_dict[each_key]==highest_freq:
                    final_read[pileupcolumn.pos]=each_key
    #remove preditermined number of nucleotides from each end
    if len(list(final_read.keys()))>(args.skip_from_end*2):
        read_start_pos=min(list(final_read.keys()))
        read_end_pos=max(list(final_read.keys()))
        for i in range(read_start_pos,read_start_pos+args.skip_from_end):
            if i in final_read:
                del final_read[i]
        for i in range(read_end_pos-args.skip_from_end+1, read_end_pos+1):
            if i in final_read:
                del final_read[i]
        if read_end_pos-read_start_pos+1>=args.cycle_number*2-args.skip_from_end*2:
            for i in range(read_start_pos+args.cycle_number-args.skip_from_end,read_start_pos+args.cycle_number):
                if i in final_read:
                    del final_read[i]
            for i in range(read_end_pos-args.cycle_number+1,read_end_pos-args.cycle_number+args.skip_from_end+1):
                if i in final_read:
                    del final_read[i]
    else:
        final_read={}
    return(final_read,n_reads)

sample_dict={
        'umi_change': 0,
        'umi_n': 0,
        'dup_mut': 0,
        'dup_snp': 0,
        'seq_len': 0,
    }

n_reads=0
n_reads_with_mutations=0

for duplex_coord in duplex_positions_file:
    chrom,pos_start,pos_end=duplex_coord.strip('\n').split(' ')
    read_dict={
        'umi_change': 0,
        'umi_n': 0,
        'dup_mut': 0,
        'dup_snp': 0,
        'seq_len': 0,
        'r1_reads': 0,
        'r2_reads': 0,
    }
    mutations_to_file=[]
    if not chrom in vcf_chr_dict:
        continue
    #print(chrom,pos_start,pos_end)
    final_read1,read_dict['r1_reads']=pileup_reads(read1_samfile.pileup(contig=chrom, start=int(pos_start), \
                                                            end=int(pos_end), stepper='nofilter'))
    final_read2,read_dict['r2_reads']=pileup_reads(read2_samfile.pileup(contig=chrom, start=int(pos_start), \
                                                            end=int(pos_end), stepper='nofilter'))
    reference_seq=ref_seq.fetch(reference=chrom,start=int(pos_start),end=int(pos_end))
    reference_seq=reference_seq.upper()
    for i in range(0,len(reference_seq)):
        pos=i+int(pos_start)
        read_dict['seq_len']+=1
        if not (pos in final_read1) or not (pos in final_read2) or not (reference_seq[i]=="A" \
                or reference_seq[i]=="T" or reference_seq[i]=="G" or reference_seq[i]=="C"):
            read_dict['seq_len']-=1
            continue
        elif (((final_read1[pos]==reference_seq[i]) and (final_read2[pos]==reference_seq[i])) or \
            ((final_read1[pos]=='*') or (final_read2[pos]=='*'))):
            continue
        elif (final_read1[pos]=='N') or (final_read2[pos]=='N') or\
        (final_read1[pos]!=reference_seq[i] and final_read2[pos]!=reference_seq[i] and\
        final_read1[pos]!=final_read2[pos]):
            read_dict['umi_n']+=1
            continue
        else:
            if final_read1[pos]!=reference_seq[i]:
                alt_base=final_read1[pos]
            elif final_read2[pos]!=reference_seq[i]:
                alt_base=final_read2[pos]
            #format alleles for fetching in vcf database
            if re.search('-\d+N+',alt_base):
                num_del=(re.search('\d+',alt_base)).group()
                alt_base=re.sub('-\d+N+','',alt_base)
                ref_base=reference_seq[i:i+int(num_del)+1]
            else:
                alt_base=re.sub('\+\d+','',alt_base)
                ref_base=reference_seq[i]
            alleles_tuple=(ref_base,alt_base)
            #fetch database
            vcf_records=vcf_file.fetch(contig=vcf_chr_dict[chrom],start=pos, stop=pos+1)
            mut_annotated=0
            vcf_record_info_dict={}
            vcf_record_info_list=[]
            for vcf_record in vcf_records:
                if vcf_record.alleles==alleles_tuple:
                    mut_annotated=1
                    #snp found -- find rs and AF
                    vcf_record_info_list=str(vcf_record).strip('\n').split('\t')
                    for info_piece in vcf_record_info_list[7].split(';'):
                        info=info_piece.split('=')
                        if len(info)>1:
                            vcf_record_info_dict[info[0]]=info[1]
                            
            if final_read1[pos]==final_read2[pos]:
                if mut_annotated==1:
                    if float(vcf_record_info_dict.get('AF',0))<args.max_af:
                        read_dict['dup_mut']+=1
                        mutations_to_file.append('\t'.join([args.sample_label,chrom,str(pos),alleles_tuple[0],\
                            alleles_tuple[1],"dup",str(pos_start),str(pos_end),vcf_record_info_list[2],\
                            vcf_record_info_dict.get('AF','0')]))
                    else:
                        read_dict['dup_snp']+=1
                        mutations_to_file.append('\t'.join([args.sample_label,chrom,str(pos),alleles_tuple[0],alleles_tuple[1],\
                          "snp",str(pos_start),str(pos_end),vcf_record_info_list[2],vcf_record_info_dict.get('AF','0')]))
                else:
                    mutations_to_file.append('\t'.join([args.sample_label,chrom,str(pos),alleles_tuple[0],alleles_tuple[1],\
                          "dup",str(pos_start),str(pos_end),'.','0']))
                    read_dict['dup_mut']+=1
            else:
                read_dict['umi_change']+=1
                if mut_annotated==1:
                    mutations_to_file.append('\t'.join([args.sample_label,chrom,str(pos),alleles_tuple[0],alleles_tuple[1],\
                          "umi",str(pos_start),str(pos_end),vcf_record_info_list[2],vcf_record_info_dict.get('AF','0')]))
                else:
                    mutations_to_file.append('\t'.join([args.sample_label,chrom,str(pos),alleles_tuple[0],alleles_tuple[1],\
                          "umi",str(pos_start),str(pos_end),'.','0']))
    if not read_dict['dup_mut']>=args.mutation_limit:       
        n_reads+=1
        if read_dict['dup_mut']>0:
            n_reads_with_mutations+=1
        
        for each_record in mutations_to_file:
            print('\t'.join([each_record,str(read_dict['dup_mut']),str(read_dict['dup_snp']),\
                  str(read_dict['umi_change']),str(read_dict['r1_reads']),str(read_dict['r2_reads'])]),\
                  file=position_statistics_file)
    
        for feature in sample_dict.keys():
            sample_dict[feature]+=read_dict[feature]

print(sample_dict['dup_mut'], "mutations were found in",n_reads_with_mutations,"templates out of",n_reads)
print(sample_dict['seq_len'], "bp were successfully analysed by duplex sequencing")
print("Among these",sample_dict['umi_change'],"changes affected only one strand")
print(sample_dict['dup_snp'],"were polymorphisms with population frequencies higher than",args.max_af)
print(sample_dict['umi_n'],"positions were excluded because it was impossible to build the consensus sequence")
print("Results were added to",args.position_statistics_file)

ref_seq.close()
vcf_file.close()
read1_samfile.close()
read2_samfile.close()
duplex_positions_file.close()
position_statistics_file.close()
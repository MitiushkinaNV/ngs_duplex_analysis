import argparse
parser=argparse.ArgumentParser()

parser.add_argument("-f","--file_in",help="File with depths of coverage", required=True)
parser.add_argument("-o","--file_out",help="File with overcovered regions", default="overcovered.txt")
parser.add_argument("-n","--num_samples",\
    help="Number of samples in which the region has to be covered to filter all reads from this region", \
    required=True, type=int)

args=parser.parse_args()

depth_file=open(args.file_in,'rt')
overcovered_regions=open(args.file_out,'wt')

region_start=0
region_stop=0
current_chr=""

for line in depth_file:
    line_values=line.split("\t")
    chrom=line_values.pop(0)
    position=int(line_values.pop(0))
    samples_covered=0
    for line_value in line_values:
        if int(line_value)>0:
            samples_covered+=1
    if samples_covered>=args.num_samples:
        if current_chr:
            if (current_chr==chrom) and ((position-region_stop)<50):
                region_stop=position
                continue
            else:
                print(current_chr,region_start,region_stop,file=overcovered_regions)
                current_chr=chrom
                #move from 1-based to 0-based right-exclusive intervals
                region_start=position-1
                region_stop=position
        else:
            current_chr=chrom
            region_start=position-1
            region_stop=position
            
print(current_chr,region_start,region_stop,file=overcovered_regions)
    
depth_file.close()
overcovered_regions.close()

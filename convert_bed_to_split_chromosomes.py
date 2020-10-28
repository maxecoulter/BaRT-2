import optparse
from optparse import OptionParser
parser = OptionParser()
(options, args) = parser.parse_args()

infile=args[0]
outfile1=args[1]

#infile="/mnt/shared/scratch/mc42302/201903_RTD2/Pacbio_20_samples/primers_for_RT_PCR_Barke.bed"
conversiondict={"chr1H":205415724,"chr2H":295953730,"chr3H":267034671,"chr4H":273224149,"chr5H":209770093,"chr6H":246027230,"chr7H":316571039}

outfile=open(outfile1,"w")
for line in open(infile):
    line1=line.split("\t")
    midline="\t".join(line1[3:6])
    endline="\t".join(line1[8:])
    chromosome=line1[0]
    start=int(line1[1])
    end=int(line1[2])
    try:
        splitpos=int(conversiondict[chromosome])
        if start>=splitpos:
            newchrom=chromosome+"_part_2"
            newstart=str(int(start)-int(splitpos))
        elif start<splitpos:
            newchrom=chromosome+"_part_1"
            newstart=str(start)
        else:
            print("Error!")
        if end>=splitpos:
            newchrom=chromosome+"_part_2"
            newend=str(int(end)-int(splitpos))
        elif end<splitpos:
            newchrom=chromosome+"_part_1"
            newend=str(end)
        else:
            print("Error!")
        if len(line1)==4 or len(line1)==3:
            try:
                identity=line1[3]
            except IndexError:
                identity="\n"
            outfile.write(newchrom+"\t"+newstart+"\t"+newend+"\t"+identity)
        elif len(line1)==12:
            outfile.write(newchrom+"\t"+newstart+"\t"+newend+"\t"+midline+"\t"+newstart+"\t"+newend+"\t"+endline)
        else:
            print("Error! splitter only works with bed3,bed4 or bed12 formats")
            break
    except KeyError:
        outfile.write(line)

outfile.close()
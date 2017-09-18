"""

Converts 23AndMe to a quantstudio-ish format and handles some genotyping edge cases I need to make it work for my own
purposes

"""

import sys
import re

raw = sys.argv[1]  # users raw 23AndMe download
snps = sys.argv[2]  # the snps in question
gene_table = sys.argv[3] # an index of gene-rsid relations
name = sys.argv[4] # the sample ID to assign 


#load genes into a dict
def load_genes(g):
    outDict = {}
    with open(g) as  file:
        for line in file:
            line = line.strip()
            temp = line.split('\t')
            outDict[temp[1]] = temp[0]
    return outDict


# assets for the quantstudio file
header = """
# Exported By : GUEST																										
# Export Date : 07/06/2017 14:03:33 CDT																										
# Study Name : 2017-7-6 CH PGx 14 SNP																										
# Experiment Type : Real-time																										
# Instrument Type : QuantStudio 12K Flex Real-Time PCR System																										
# Software Version Number : 1.3.1																										
# Creation Date : 07/06/2017 13:16:50 CDT																										
# Created By : GUEST																										
# Last Modified Date : 07/06/2017 14:00:26 CDT																										
# Last Modified By : GUEST																										
# Template File Name : N/A																										
# Template Originating Study Name : N/A																										
# Template Creation Date : N/A																										
# Template Created By User ID : N/A																										
# Template Software Version Number :  N/A																										
																										
																										
Assay Name	Assay ID	Gene Symbol	NCBI SNP Reference	Sample ID	Call	Manual	Quality	VIC(Rn)	FAM(Rn)	ROX	Task	Gender	Population	Well	Experiment Name	Plate Barcode	Annotation	Low ROX Intensity	NTC FAM Intensity High	NTC VIC Intensity High	Genotype Quality Low	Failed Control	Reference Sample Discordance	Replicate Sample Discordance	Chromosome #	Position
"""
snp_line = [
    "C_25625804D_20",
    "C_25625804D_20",
    "CYP2C9",
    "rs7900194",
    "012345",
    "G/G",
    "FALSE",
    "0.99",
    "138.995",
    "\"3,623.04\"",
    "865.382",
    "Unknown",
    "Unknown",
    "",
    "B1c4",
    "WTI53_2017_07_06_092106.eds",
    "WTI53",
    "",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "FALSE",
    "10",
    ""
]
tail = """










																										
Assay ID	Assay Name	Population	Allele 1 Freq	Allele 2 Freq	1/1 Freq	1/2 Freq	2/2 Freq	Chi-Squared	P-Value																	
C_25625804D_20	C_25625804D_20	All	0%	100%	0%	0%	100%	0	1																	
C_30634117C_K0	C_30634117C_K0	All	0%	100%	0%	0%	100%	0	1																	
C_30634117D_M0	C_30634117D_M0	All	0%	100%	0%	0%	100%	0	1																	
C__11484460_40	C__11484460_40	All	18.90%	81.10%	3.30%	31.10%	65.60%	0.02	0.888																	
C__11975250_10	hCV11975250	All	98.40%	1.60%	96.70%	3.30%	0%	0.017	0.896																	
C__15859191_30	C__15859191_30	All	87.10%	12.90%	77.40%	19.40%	3.20%	1.196	0.274																	
C__16017734_10	C__16017734_10	All	1.60%	98.40%	0%	3.20%	96.80%	0.017	0.896																	
C__25625804_10	C__25625804_10	All	0%	100%	0%	0%	100%	0	1																	
C__25625805_10	C__25625805_10	All	93.50%	6.50%	87.10%	12.90%	0%	0.295	0.587																	
C__25745302_30	C__25745302_30	All	0%	100%	0%	0%	100%	0	1																	
C__25746809_50	C__25746809_50	All	47.50%	52.50%	18%	59%	23%	2.047	0.153																	
C__25986767_70	C__25986767_70	All	8.10%	91.90%	0%	16.10%	83.90%	0.477	0.49																	
C__26201809_30	C__26201809_30	All	6.50%	93.50%	3.20%	6.50%	90.30%	13.436	0																	
C__27028164_10	C__27028164_10	All	48.40%	51.60%	26.20%	44.30%	29.50%	0.79	0.374																	
C__27102425_10	C__27102425_10	All	40.20%	59.80%	23%	34.40%	42.60%	4.912	0.027																	
C__27102431_D0	C__27102431_D0	All	79%	21%	64.50%	29%	6.50%	0.954	0.329																	
C__27104892_10	C__27104892_10	All	4.80%	95.20%	0%	9.70%	90.30%	0.16	0.689																	
C__27531552_A0	C__27531552_A0	All	100%	0%	100%	0%	0%	0	1																	
C__27531918_10	C__27531918_10	All	0%	100%	0%	0%	100%	0	1																	
C__27535825_20	C__27535825_20	All	98.40%	1.60%	96.80%	3.20%	0%	0.017	0.896																	
C__27830964_10	C__27830964_10	All	0%	100%	0%	0%	100%	0	1																	
C__27859817_40	C__27859817_40	All	99.20%	0.80%	98.40%	1.60%	0%	0.004	0.95																	
C__27861809_10	C__27861809_10	All	0%	100%	0%	0%	100%	0	1																	
C__27861810_10	C__27861810_10	All	100%	0%	100%	0%	0%	0	1																	
C__30203950_10	C__30203950_10	All	100%	0%	100%	0%	0%	0	1																	
C__30403261_20	C__30403261_20	All	56.60%	43.40%	36.10%	41%	23%	1.681	0.195																	
C__30633862_10	C__30633862_10	All	100%	0%	100%	0%	0%	0	1																	
C__30633871_50	C__30633871_50	All	100%	0%	100%	0%	0%	0	1																	
C__30633906_10	C__30633906_10	All	22.10%	77.90%	3.30%	37.70%	59%	0.539	0.463																	
C__30634127_10	C__30634127_10	All	0%	100%	0%	0%	100%	0	1																	
C__30634128_10	C__30634128_10	All	100%	0%	100%	0%	0%	0	1																	
C__30634130_30	C__30634130_30	All	0%	100%	0%	0%	100%	0	1																	
C__30634131_20	C__30634131_20	All	0%	100%	0%	0%	100%	0	1																	
C__30634132_70	C__30634132_70	All	0%	100%	0%	0%	100%	0	1																	
C__30634136_10	C__30634136_10	All	100%	0%	100%	0%	0%	0	1																	
C__30634146_10	C__30634146_10	All	98.40%	1.60%	96.80%	3.20%	0%	0.017	0.896																	
C__32287188_10	C__32287188_10	All	3.20%	96.80%	3.20%	0%	96.80%	62	0																	
C__32287221_20	C__32287221_20	All	0%	100%	0%	0%	100%	0	1																	
C__32388575_A0	C__32388575_A0	All	0%	100%	0%	0%	100%	0	1																	
C__32407229_60	C__32407229_60	All	100%	0%	100%	0%	0%	0	1																	
C__32407232_50	C__32407232_50	All	98.40%	1.60%	96.80%	3.20%	0%	0.017	0.896																	
C__32407243_20	C__32407243_20	All	99.20%	0.80%	98.30%	1.70%	0%	0.004	0.95																	
C__33845811_20	C__33845811_20	All	100%	0%	100%	0%	0%	0	1																	
C__34816113_20	C__34816113_20	All	98.40%	1.60%	96.70%	3.30%	0%	0.017	0.896																	
C__34816116_20	C__34816116_20	All	89.80%	10.20%	83.10%	13.60%	3.40%	3.923	0.048																	
C__59013445_10	hCV59013445	All	94.30%	5.70%	88.50%	11.50%	0%	0.226	0.635																	
C__60142977_10	C__60142977_10	All	84.40%	15.60%	75.40%	18%	6.60%	6.024	0.014																	
C___1202883_20	hCV1202883	All	61.30%	38.70%	41.90%	38.70%	19.40%	2.104	0.147																	
C___2138789_10	hCV2138789	All	56.50%	43.50%	32.30%	48.40%	19.40%	0.016	0.899																	
C___2222771_A0	C___2222771_A0	All	98.40%	1.60%	96.80%	3.20%	0%	0.017	0.896																	
C___3084793_20	hCV3084793	All	20.50%	79.50%	0%	41%	59%	4.052	0.044																	
C___7486676_10	hCV7486676	All	17.70%	82.30%	3.20%	29%	67.70%	0.002	0.964																	
C___7817765_60	C___7817765_60	All	75.90%	24.10%	53.40%	44.80%	1.70%	2.911	0.088																	
C___8303531_40	C___8303531_40	All	93.50%	6.50%	90.30%	6.50%	3.20%	13.436	0																	
C___8726802_20	hCV8726802	All	0%	100%	0%	0%	100%	0	1																	
C___8881221_40	C___8881221_40	All	25.80%	74.20%	6.50%	38.70%	54.80%	0.007	0.933																	
C___8950074_1_	hCV8950074	All	87.10%	12.90%	77.40%	19.40%	3.20%	1.196	0.274																	
C____469857_10	C____469857_10	All	81.10%	18.90%	65.60%	31.10%	3.30%	0.02	0.888																	
C____850486_20	hCV850486	All	28.70%	71.30%	6.60%	44.30%	49.20%	0.408	0.523																	
C____904973_10	hCV904973	All	89.30%	10.70%	78.70%	21.30%	0%	0.868	0.352																	
"""


# build a structure for 23AndMe
def build(file):
    outDict = {}
    with open(file) as f:
        for line in f:
            line = line.strip()
            if re.match("^#", line):
                continue
            else:
                temp = line.split('\t')
                vars = temp[3]
                try:
                    varsOut = vars[0] + "/" + vars[1]
                except:
                    varsOut = vars  # catches haploid calls (MT/Y) 
                outDict[temp[0]] = varsOut
    return outDict


# search the 23AndMe data for the relevant markers 
def search(search, seek):
    outDict = {}
    with open(seek) as f:
        for line in f:
            rsid = line.strip()
            try:
                outDict[rsid] = search[rsid]
            except:
                outDict[rsid] = "UND"
            if outDict[rsid] == "-/-":
                outDict[rsid] = "UND"  # catches null genotypes in 23AndMe data
    return outDict


# resolve edge cases
def edges(raws):
    outDict = {}
    outDict["rs2069514"] = "UND" # an rsid not seen in this version of 23AndMe
    outDict["rs35694136"] = "UND" # an rsid not seen in this version of 23AndMe
    outDict["rs35599367"] = "UND" # an rsid not seen in this version of 23AndMe
    outDict["rs1799963"] = "UND" # an rsid not seen in this version of 23AndMe, this is actually almost fixed so maybe we can assign it anyways
    outDict["rs5030862"] = "UND" # an rsid not seen in this version of 23AndMe, probably could be fixed as T though due to freqs. UND for now because can't get CNV for CYP2D6 anyways.
    
    for key, value in raws.items():
        if key == "rs41303343":  # an insertion event
            if value == "I/I":
                value = "T/T"
            elif value == "D/I":
                value = "-/T"
            elif value == "D/D":
                value = "-/-"
        if key == "rs9332131":  # a deletion event
            if value == "I/I":
                value = "A/A"
            elif value == "D/I":
                value = "-/A"
            elif value == "D/D":
                value = "-/-"
        if key == "rs5030656":  # this needs to be updated to rs869035800, but functional consequence similar so swap it in for now
            if value == "I/I":
                value = "TCT/TCT"
            elif value == "D/I":
                value = "-/TCT"
            elif value == "D/D":
                value = "-/-"
        if key == "rs35742686":  # this is a multiallelic loci that has relevant point or indel changes, 23AndMe seems to only consider indel
            if value == "I/I":
                value = "T/T"
            elif value == "D/I":
                value = "-/T"
            elif value == "D/D":
                value = "-/-"
        if key == "rs5030655":  # this is a multiallelic loci that has relevant point or indel changes, 23AndMe seems to only consider indel
            if value == "I/I":
                value = "A/A"
            elif value == "D/I":
                value = "-/A"
            elif value == "D/D":
                value = "-/-"  
        outDict[key] = value
    return outDict


# convert the results to the quant studio format    
def convert(edged,genes):
    outsnps = header
    for key,value in edged.items():
        line = snp_line
        line[2] = genes[key]
        line[3] = key
        line[4] = name
        line[5] = value
        outline = '\t'.join(line)
        outsnps += outline+"\n"
    outsnps += tail
    return outsnps
        


if __name__ == "__main__":
    genes = load_genes(gene_table)
    searcher = build(raw)
    raw_results = search(searcher, snps)
    print(raw_results)
    edged = edges(raw_results)
    print(edged)
    snp_write = convert(edged,genes)
    outfile = open(name+".snps.txt",'w')
    outfile.write(snp_write)
    
    

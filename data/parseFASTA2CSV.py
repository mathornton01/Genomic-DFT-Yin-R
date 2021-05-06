filename = "gisaid_hcov-19_2020_08_20_02.fasta"
ofilename = "gisaid_hcov-19_2020_08_20_02_unaligned.csv"
headers = []
sequences = [] 
newseq = False
seq = ""
with open(filename,'r') as fastain: 
    for line in fastain: 
        if line[0] == ">": 
            headers.append(line.strip())
            if (seq!=""):
                sequences.append(seq)
            newseq = True
        else:
            if newseq:
                newseq = False
                seq = ""
            else: 
                seq = seq+line.strip()

sequences.append(seq)

with open(ofilename,'w') as csvout: 
    csvout.write("Header,Sequence\n");
    for i in range(len(headers)):
        csvline = headers[i] + "," + sequences[i] + "\n";
        csvout.write(csvline)

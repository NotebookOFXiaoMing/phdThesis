import sys

fr = open(sys.argv[1],'r')
fw = open(sys.argv[2],'w')

for line in fr:
    if line.split()[2] == "transcript":
        temp_list = line.strip().split()
        temp_list[2] = "gene"
        fw.write("%s\n"%('\t'.join(temp_list)))
        
    else:
        fw.write(line)
        
fr.close()
fw.close()
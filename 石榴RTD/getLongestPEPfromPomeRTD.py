import sys

gene_dict = {}

#a = 0
with open(sys.argv[1],'r') as fr:
    for line in fr:
        #a += 1
        if line.startswith("pome"):
            temp_list = line.strip().split(",")
            #print(a)
            if temp_list[2] == "Coding":
                if temp_list[0] not in gene_dict.keys():
                    gene_dict[temp_list[0]] = temp_list[-1]
                if temp_list[0] in gene_dict.keys():
                    if len(temp_list[-1]) > len(gene_dict[temp_list[0]]):
                        gene_dict[temp_list[0]] = temp_list[-1]
                        

fw = open(sys.argv[2],'w')
for i,j in gene_dict.items():
    fw.write(">%s\n%s\n"%(i,j))
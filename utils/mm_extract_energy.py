def generate_list(file_path: str) -> str:
    handle = open(file_path,'r')
    line_no = 0
    primer_2mer =""
    tmp_2mer =""
    list_str = "mm_dh = deque([\n"
    gt_pair = 0
    for line in handle:
        
        line_no += 1
        # 16行一个矩阵 遇到GT pair 时跳过
        no = line_no % 14
        
        if no == 2:
            primer_2mer = line.rstrip().replace("\t","").replace(" ","")
            # print(primer_2mer,end="\n")
        if no == 3:
            tmp_2mer = line.rstrip().replace("\t","").replace(" ","")
            # print(tmp_2mer,end="\n")
            # skip GT pair 71 + 84
            for base1, base2 in zip(primer_2mer, tmp_2mer):
                if ord(base1) + ord(base2) == 155:
                    gt_pair = 1
            if not gt_pair:
                list_str += "\t# "
                list_str += " ".join(list(primer_2mer))
                list_str += "\n"
                list_str += "\t# "
                list_str += " ".join(list(tmp_2mer))
                list_str += "\n"
                list_str += "\t[\n"
        if not gt_pair and 9 <= no <= 12:
            list_str += "\t\t"
            list_str += line[2:].replace("\t",",").rstrip()
            list_str += ",\n"
            # print(line[2:],end="")
        if not gt_pair and no == 13:
            list_str += "\t],\n"
        if no == 13:
            gt_pair = 0
    list_str += "])"
    print(list_str)
                    
                    
        


def main() -> None:
    generate_list("/mnt/ntc_data/wayne/Repositories/primer3_opt/utils/dna.tstackm.dh")
    

if __name__ == "__main__":
    main()
import pandas as pd


def base2int(base:str) -> int:
        trantab = str.maketrans('ACGTN', '01234')
        return int(base.upper().translate(trantab), base=5)


def fill_N(df:pd.DataFrame) -> dict:
    dh_dict = {}
    dh_li_m = [79, 84, 78, 72, 72, 85, 80, 106, 78, 78, 82, 98, 80, 84, 80, 72, 82, 85, 79, 72, 72, 80, 78, 72, 72]
    # 每四个碱基根据平均值补一个N
    base_li = ['A','C','G','T']
    basen_li = base_li + ['N']
    for i in base_li:
        for j in base_li:
            sum_2 = 0
            for k in base_li:
                # 补N
                sum_1 = 0
                for l in base_li:
                    duplex = f"{i}{j}/{k}{l}"
                    dh = df.loc[df[0] == duplex][1]
                    ds = df.loc[df[0] == duplex][2]
                    dh_data = int(dh.values[0]) if len(dh) else dh_li_m[base2int(duplex[:2])] * -100
                    sum_1 += dh_data
                    dh_dict[duplex] = dh_data
                    # print(f"{duplex}\t{dh_data}")
                
                sum_2 += int(sum_1/4)
                dh_dict[f"{duplex[:4]}N"] = int(sum_1/4)
                
            # /NA /NC /NG /NT /NN
            for l2 in base_li:
                # NA = (AA + CA + GA + TA )/4
                sum = 0
                for l3 in base_li:
                    sum += dh_dict[f"{i}{j}/{l3}{l2}"]
                dh_dict[f"{i}{j}/N{l2}"] = int(sum/4)
            dh_dict[f"{i}{j}/NN"] = int(sum_2/4)
    
        # XN/XX
        for x in basen_li:
            for y in basen_li:
                dh_dict[f"{i}N/{x}{y}"] = int((dh_dict[f"{i}A/{x}{y}"] + dh_dict[f"{i}C/{x}{y}"] + dh_dict[f"{i}G/{x}{y}"] + dh_dict[f"{i}T/{x}{y}"])/4)
        
    #NX/XX
    for x1 in basen_li:
        for y1 in basen_li:
            for z1 in basen_li:
                dh_dict[f"N{x1}/{y1}{z1}"] = int((dh_dict[f"A{x1}/{y1}{z1}"] + dh_dict[f"C{x1}/{y1}{z1}"] + dh_dict[f"T{x1}/{y1}{z1}"] + dh_dict[f"G{x1}/{y1}{z1}"])/4)
    
    return dh_dict

def print_format_li(hy_dict:dict) -> None:
    ct = 1
    comment = "# "
    data_li = []
    print("[",end="")
    for k,v in hy_dict.items():
        comment += f"{k} "
        data_li.append(str(v))
        if ct % 5 == 0:
            print(f",".join(data_li),end="")
            print(f", {comment}")
            comment = "# "
            data_li = []
        ct += 1
    print("]")
        
    
                    

def main():
    df = pd.read_csv("mm.txt",header=None,sep=r'\s+')
    dict1 = fill_N(df)
    print_format_li(dict1)
    # print(dict1)


if __name__ == "__main__":
    main()
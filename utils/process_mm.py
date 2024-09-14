import pandas as pd
from decimal import Decimal


def base2int(base:str) -> int:
        trantab = str.maketrans('ACGTN', '01234')
        return int(base.upper().translate(trantab), base=5)


def fill_N(df:pd.DataFrame) -> dict:
    dh_dict = {}
    dh_li_m = [79, 84, 78, 72, 72, 85, 80, 106, 78, 78, 82, 98, 80, 84, 80, 72, 82, 85, 79, 72, 72, 80, 78, 72, 72]
    ds_li_m = [222, 224, 210, 204, 224, 227, 199, 272, 210, 272, 222, 244, 199, 224, 244, 213, 222, 227, 222, 227, 168, 210, 220, 215, 220]
    # 每四个碱基根据平均值补一个N
    base_li = ['A','C','G','T']
    basen_li = base_li + ['N','-']
    for i in base_li:
        for j in base_li:
            t_max = 0
            for k in base_li:
                # 补N
                max_1 = 0
                for l in base_li:
                    duplex = f"{i}{j}/{k}{l}"
                    dh = df.loc[df[0] == duplex][1]
                    dh_data = dh.values[0] if len(dh) else dh_li_m[base2int(duplex[:2])] * -100
                    # dh = df.loc[df[0] == duplex][2]
                    # dh_data = dh.values[0] if len(dh) else ds_li_m[base2int(duplex[:2])] * -0.1
                    if dh_data > max_1: max_1 = round(dh_data, 1)
                    dh_dict[duplex] = round(dh_data, 1)
                    # print(f"{duplex}\t{dh_data}")
                if max_1 > t_max:t_max = round(max_1, 1)              
                dh_dict[f"{duplex[:4]}N"] = round(max_1, 1)
                # 补-
                dh_dict[f"{duplex[:4]}-"] = dh_dict[f"{duplex[:4]}N"]
                
            # /NA /NC /NG /NT /NN
            for l2 in base_li:
                # NA = (AA + CA + GA + TA )/4
                sum = 0
                max_0 = 0
                for l3 in base_li:
                    sum += dh_dict[f"{i}{j}/{l3}{l2}"]
                    if dh_dict[f"{i}{j}/{l3}{l2}"] > max_0: max_0 = round(dh_dict[f"{i}{j}/{l3}{l2}"],1)
                if max_0 > t_max:t_max = round(max_0, 1)
                dh_dict[f"{i}{j}/N{l2}"] = max_0
            dh_dict[f"{i}{j}/NN"] = t_max
            dh_dict[f"{i}{j}/N-"] = t_max
            for l2 in base_li:
                dh_dict[f"{i}{j}/-{l2}"] = dh_dict[f"{i}{j}/N{l2}"]
            dh_dict[f"{i}{j}/-N"] = t_max           
            dh_dict[f"{i}{j}/--"] = t_max
    
        # XN/XX
        for x in basen_li:
            for y in basen_li:
                dh_dict[f"{i}N/{x}{y}"] = round(min((dh_dict[f"{i}A/{x}{y}"], dh_dict[f"{i}C/{x}{y}"], dh_dict[f"{i}G/{x}{y}"], dh_dict[f"{i}T/{x}{y}"])), 1)
        for x in basen_li:
            for y in basen_li:
                dh_dict[f"{i}-/{x}{y}"] = round(max((dh_dict[f"{i}A/{x}{y}"], dh_dict[f"{i}C/{x}{y}"], dh_dict[f"{i}G/{x}{y}"], dh_dict[f"{i}T/{x}{y}"])), 1)
        
    #NX/XX
    for x1 in basen_li:
        for y1 in basen_li:
            for z1 in basen_li:
                dh_dict[f"N{x1}/{y1}{z1}"] = round(min(dh_dict[f"A{x1}/{y1}{z1}"],dh_dict[f"C{x1}/{y1}{z1}"], dh_dict[f"T{x1}/{y1}{z1}"],dh_dict[f"G{x1}/{y1}{z1}"]), 1)
    for x1 in basen_li:
        for y1 in basen_li:
            for z1 in basen_li:
                dh_dict[f"-{x1}/{y1}{z1}"] = round(max((dh_dict[f"{i}A/{x}{y}"], dh_dict[f"{i}C/{x}{y}"], dh_dict[f"{i}G/{x}{y}"], dh_dict[f"{i}T/{x}{y}"])), 1)
    
    return dh_dict

def print_format_li(hy_dict:dict) -> None:
    ct = 1
    comment = "# "
    data_li = []
    print("\t\tself.dX_mm = [")
    for k,v in hy_dict.items():
        if ct % 6 == 1:print("\t\t",end="")
        comment += f"{k} "
        data_li.append(str(v))
        if ct % 6 == 0:
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
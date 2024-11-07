from collections import defaultdict, deque


int_loop_energy_dict = defaultdict(deque)


int_loop_energy_dict["11dh"] = deque([
        # A A
        # T T
        [
            1.3, 2.2, 0.9, 1.0, # AA AC AG AT
            1.7, 2.4, 1.0, 1.4, # CA CC CG CT
            1.7, 2.4, 1.0, 1.4, # GA GC GG GT
            1.7, 2.4, 1.0, 1.4, # TA TC TG TT           
        ],
        # A C
        # T G
        [
            0.8,  1.4, -0.4, 1.0, # AA AC AG AT
            1.6,  2.1,  1.0, 1.6, # CA CC CG CT
            -0.2, 1.0, -1.2, 2.0, # GA GC GG GT
            1.0,  1.4,  2.0, 1.1, # TA TC TG TT    
        ],
        # A G
        # T C
        [
            1.0, 1.7,  0.3, 1.0, # AA AC AG AT
            1.5, 2.0,  1.0, 1.0, # CA CC CG CT
            0.1, 1.0, -0.3, 2.0, # GA GC GG GT
            1.0, 1.4,  2.0, 0.6, # TA TC TG TT 
        ],
        # A T
        # T A
        [
            1.2, 1.7,  0.2, 1.0, # AA AC AG AT
            1.7, 2.7,  1.0, 1.4, # CA CC CG CT
            0.2, 1.0, -0.3, 2.0, # GA GC GG GT
            1.0, 1.4,  2.0, 1.4, # TA TC TG TT 
        ],
        # C A
        # G T
        [
            1.1, 2.1, 0.8, 1.0, # AA AC AG AT
            1.7, 1.8, 1.0, 1.4, # CA CC CG CT
            0.5, 1.0, 0.3, 2.0, # GA GC GG GT
            1.0, 1.4, 2.0, 0.6, # TA TC TG TT 
        ],
        # C C
        # G G
        [
            0.6,  1.2, -0.5, 1.0, # AA AC AG AT
            1.6,  1.5,  1.0, 1.6, # CA CC CG CT
            -0.1, 1.0, -1.3, 2.0, # GA GC GG GT
            1.0,  1.0,  2.0, 0.3, # TA TC TG TT 
        ],
        # C G
        # G C
        [
            0.9, 1.5,  0.1,  1.0, # AA AC AG AT
            1.5, 1.4,  1.0,  1.0, # CA CC CG CT
            0.1, 1.0, -0.3,  2.0, # GA GC GG GT
            1.0, 1.0,  2.0, -0.2, # TA TC TG TT 
        ],
        # C T
        # G A
        [
            1.0, 1.5,  0.1, 1.0, # AA AC AG AT
            1.7, 2.0,  1.0, 1.4, # CA CC CG CT
            0.3, 1.0, -0.3, 2.0, # GA GC GG GT
            1.0, 1.0,  2.0, 0.6, # TA TC TG TT 
        ],
        # G A
        # C T
        [
            0.9,  2.1,  0.5, 1.0, # AA AC AG AT
            1.4,  1.8,  1.0, 1.4, # CA CC CG CT
            -0.1, 1.0, -0.7, 2.0, # GA GC GG GT
            1.0,  2.0,  2.0, 1.1, # TA TC TG TT 
        ],
        # G C
        # C G
        [
            0.3,  1.3, -0.8, 1.0, # AA AC AG AT
            1.3,  1.6,  1.0, 1.6, # CA CC CG CT
            -0.8, 1.0, -2.2, 2.0, # GA GC GG GT
            1.0,  1.6,  2.0, 0.9, # TA TC TG TT 
        ],
        # G G
        # C C
        [
            0.6,  1.6, -0.1, 1.0, # AA AC AG AT
            1.2,  1.5,  1.0, 1.0, # CA CC CG CT
            -0.5, 1.0, -1.3, 2.0, # GA GC GG GT
            1.0,  1.6,  2.0, 0.3, # TA TC TG TT 
        ],
        # G T
        # C A
        [
            0.8,  1.6, -0.2, 1.0, # AA AC AG AT
            1.4,  2.1,  1.0, 1.4, # CA CC CG CT
            -0.4, 1.0, -1.2, 2.0, # GA GC GG GT
            1.0,  1.6,  2.0, 1.1, # TA TC TG TT 
        ],
        # T A
        # A T
        [
            1.4, 2.3, 1.2, 1.0, # AA AC AG AT
            2.3, 2.1, 1.0, 1.7, # CA CC CG CT
            1.2, 1.0, 0.9, 2.0, # GA GC GG GT
            1.0, 1.7, 2.0, 1.4, # TA TC TG TT 
        ],
        # T C
        # A G
        [
            0.9, 1.4, -0.1, 1.0, # AA AC AG AT
            2.1, 1.8,  1.0, 2.0, # CA CC CG CT
            0.5, 1.0, -0.7, 2.0, # GA GC GG GT
            1.0, 1.4,  2.0, 1.1, # TA TC TG TT 
        ],
        # T G
        # A C
        [
            1.1, 1.7, 0.5, 1.0, # AA AC AG AT
            2.1, 1.8, 1.0, 1.4, # CA CC CG CT
            0.8, 1.0, 0.3, 2.0, # GA GC GG GT
            1.0, 1.4, 2.0, 0.6, # TA TC TG TT 
        ],
        # T T
        # A A
        [
            1.3, 1.7, 0.4, 1.0, # AA AC AG AT
            2.2, 2.4, 1.0, 1.7, # CA CC CG CT
            0.9, 1.0, 0.3, 2.0, # GA GC GG GT
            1.0, 1.4, 2.0, 1.4, # TA TC TG TT 
        ]      
    ])
    
matrix = int_loop_energy_dict["11dh"]


def base2int(base: str) -> int:
    """将base(without N)转化为索引号"""
    trantab = str.maketrans("ACGTN", "01234")
    return int(base.upper().translate(trantab), base=5)
def index_intl11(bases: str) -> int:
    """            
            X
           A  A  ---- AAA
           T  T
            YA
    """
    upstream_base, downstream_base = bases[:]
    external_index = base2int(upstream_base + downstream_base)
    return external_index

def obtain_N_sub_idx(bases: str) -> list:
    if index_intl11(bases) < 20 :
        sub_li = [bases[0] + x  for x in ["A", "C", "G", "T"]]
        # print([index_intl21(x) for x in sub_li])  
        return  [index_intl11(x) for x in sub_li]
    else:
        sub_li = [ x + bases[1]  for x in ["A", "C", "G", "T"]]
        return  [index_intl11(x) for x in sub_li]
    
def obtain_mean(bases: str) -> list:
    list2append = []
    idx_li = obtain_N_sub_idx(bases)
    
    for i in range(16):
        ele = round(sum([matrix[x][i] for x in idx_li]) / 4 , 1)
        list2append.append(ele)
    return list2append

# insert matrix
insert_idx = 4
# AN
int_loop_energy_dict["11dh"].insert(insert_idx, obtain_mean("AN"))
insert_idx += 5
int_loop_energy_dict["11dh"].insert(insert_idx, obtain_mean("CN"))
insert_idx += 5
int_loop_energy_dict["11dh"].insert(insert_idx, obtain_mean("GN"))
insert_idx += 5
int_loop_energy_dict["11dh"].insert(insert_idx, obtain_mean("TN"))
insert_idx += 1
# NA
int_loop_energy_dict["11dh"].append(obtain_mean("NA"))
int_loop_energy_dict["11dh"].append(obtain_mean("NC"))
int_loop_energy_dict["11dh"].append(obtain_mean("NG"))
int_loop_energy_dict["11dh"].append(obtain_mean("NT"))
int_loop_energy_dict["11dh"].append(obtain_mean("NN"))

# format print
basen_li = ["A", "C", "G", "T", "N"]
for i, x in enumerate(int_loop_energy_dict["11dh"]):
    first_base = basen_li[i // 5]
    second_base = basen_li[i % 5]
    print(f"\t\t# {first_base} {second_base}")
    
    print("\t\t[")
    print("\t\t\t",end="")
    print(", ".join([str(y) for y in x]))
    print("\t\t],")
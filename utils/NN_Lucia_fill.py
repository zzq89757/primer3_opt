from collections import defaultdict, deque
from math import log, sqrt


def complement_base(base: str) -> str:
    trantab = str.maketrans("ACGTNRYMKSWVBHDZ", "TGCANYRKMSWBVDHZ")
    return base.upper().translate(trantab)



def is_complement(base1: str, base2: str) -> bool:
    return complement_base(base1) == base2.upper()


def modify_pair_base(primer_degenerate_base: str, template_degenerate_base: str) -> list:
    base_li = ["A", "T", "C", "G"]
    if (primer_degenerate_base in base_li or primer_degenerate_base == "N") and template_degenerate_base in base_li:
        return [primer_degenerate_base, template_degenerate_base]
    

    # degenerate bases table
    degenerate_base_dict = defaultdict(deque)
    def fill_dict(i: str) -> None:
        degenerate_base_dict[i] = deque([i])
    [fill_dict(i) for i in base_li]
        
    
    degenerate_base_dict['R'] = deque(["A", "G"])
    degenerate_base_dict['Y'] = deque(["C", "T"])
    degenerate_base_dict['M'] = deque(["C", "A"])
    degenerate_base_dict['K'] = deque(["T", "G"])
    degenerate_base_dict['S'] = deque(["C", "G"])
    degenerate_base_dict['W'] = deque(["T", "A"])
    
    degenerate_base_dict['B'] = deque(["G", "T", "C"])
    degenerate_base_dict['V'] = deque(["G", "A", "C"])
    degenerate_base_dict['D'] = deque(["G", "A", "T"])
    degenerate_base_dict['H'] = deque(["C", "A", "T"])
    
    degenerate_base_dict['N'] = deque(["C", "A", "T", "G"])
    
    degenerate_base_dict['Z'] = deque([])
    # print(degenerate_base_dict)
    
    # is primer_degenerate_base contain complement to template_degenerate_base
    primer_base_li = degenerate_base_dict[primer_degenerate_base]
    template_base_li = degenerate_base_dict[template_degenerate_base]
    
    if len(primer_base_li) < len(template_base_li):
        if len(primer_base_li) == 1:
            return [primer_degenerate_base, "Z"]
        else:
            return ["Z", "Z"]
    
    primer_base_set = set(primer_base_li)
    tmp_base_set = set([complement_base(x) for x in template_base_li])
    if len(template_base_li) == 1:
        if tmp_base_set.issubset(primer_base_set):
            return ["N", template_degenerate_base]
        else:
            return ["Z", template_degenerate_base]
    else:
        if tmp_base_set.issubset(primer_base_set):
            if tmp_base_set:
                return ["N", "N"]
            else:
                if len(primer_degenerate_base) == 1:
                    # return ['X']
                    return [primer_degenerate_base, template_degenerate_base]
                else:
                    return ['X']
        else:
            return ["Z", "Z"]
    
    
# for i in "ACGTNRYMKSWVBHD":
#     for j in "ACGTNRYMKSWVBHD":
#         print([i,j],end=" ===>> ")
#         print(modify_pair_base(i, j))



def modified_duplex_str(duplex_str: str) -> str:
    seq1, seq2 = duplex_str.split("\n")
    modified_seq1 = modified_seq2 = ""
    for primer_base, tmp_base in zip(seq1.upper(),seq2.upper()):
        # primer seq,
        modified_primer_base, modified_tmp_base = modify_pair_base(primer_base, tmp_base)
        modified_seq1 += modified_primer_base
        modified_seq2 += modified_tmp_base           
        
    return f"{modified_seq1}\n{modified_seq2}"

ds = "AYNTCGATNBYGT\nTTNAGCYAYDCGT"
print(modified_duplex_str(ds))



def pair_type(ord_sum: int) -> int:
    # gap for -1
    if ord_sum % 45 == 0:
        return -1
    # match for 1
    if ord_sum in [5460, 4757] or ord_sum % 73 == 0:
        return 1
    # 0 for mismatch
    return 0


def loop_detective(duplex_str: str) -> defaultdict:
    """
    检测分子杂交(duplex)中的loop结构(internal loop 和 bulge loop)
    """
    seq1, seq2 = duplex_str.split("\n")
    # ord sum is  AT: 65 + 84 = 149   CG: 67 + 71 = 138
    match_asc_li = [149, 138]
    gap_asc_li = [90, 129, 116, 112, 110] # -- -T -G -C -A
    ord_sum_li = [ord(x) + ord(y) for x, y in zip(seq1, seq2)]
    # record mismatch and gaps,then recogonize bulge loop and internal loop
    loop_region_dict = defaultdict(
        deque
    )  # region_pos: start and end index, region_type: 0 for bulge,1 for int loop
    region_type_flag = 0  # 0 for bulge,1 for int loop
    flag = 0  # switch to record region start and end index
    # print(ord_sum_li)
    for i, v in enumerate(ord_sum_li):
        # mismatch or gap,record region start
        if v not in match_asc_li and not flag:
            loop_region_dict["region_pos"].append(i)
            flag = 1
        # mismatch found,record int loop
        if v not in gap_asc_li and v not in match_asc_li:
            region_type_flag = 1
        # record region end
        if v in match_asc_li and flag:
            loop_region_dict["region_pos"].append(i - 1)
            loop_region_dict["region_type"].append(region_type_flag)
            flag = 0
            region_type_flag = 0
    # print(loop_region_dict)
    return loop_region_dict
duplex_str = "CAGACG\nGTAGGC"
duplex_str = "CA--G--CG\nGTGAAAGGC"
duplex_str = "CA-GACG\nGTGAGGC"
duplex_str = "GCCCG\nCGG-C"
duplex_str = "GAACAG\nCT---C"
duplex_str = "TAAACTGCC-----GGCAGTACATC\nATTTGACGGGTGAACCGTCATGTAG"
# duplex_str = "CCCTATT---------------------------ATTGACGTCAATA\nGGGATAACCGCAATGATACCCTTGTATGCAGTAATAACTGCAGTTAT"
# duplex_str = "GGCCGGAGTAAGCTGACAT\nCCGGCCTCATTCGACTGTA"
# duplex_str = "AAAAAAAAAA\nTTTTTTTTTT"
region_dict = loop_detective(duplex_str)
print(region_dict)
exit()

def base2int(base: str) -> int:
    """将base(without N)转化为索引号"""
    trantab = str.maketrans("ACGT", "0123")
    return int(base.upper().translate(trantab), base=4)


def stack_energy(segment: str) -> list:
    """近邻法计算stack的总能量值"""
    delta_h = [
        -7.9, -8.4, -7.8,  -7.2, # AA AC AG AT
        -8.5, -8.0, -10.6, -7.8, # CA CC CG CT
        -8.2, -9.8, -8.0,  -8.4, # GA GC GG GT
        -7.2, -8.2, -8.5,  -7.9, # TA TC TG TT
    ]
    delta_s = [
        -22.2, -22.4, -21.0, -20.4, # AA AC AG AT
        -22.7, -19.9, -27.2, -21.0, # CA CC CG CT
        -22.2, -24.4, -19.9, -22.4, # GA GC GG GT
        -21.3, -22.2, -22.7, -22.2, # TA TC TG TT
    ]
    stack_dh = 0
    stack_ds = 0
    stack_dg = 0
    # print(delta_g)
    for i in range(len(segment) - 1):
        two_mer = segment[i] + segment[i + 1]
        d_index = base2int(two_mer)
        stack_dh += delta_h[d_index]
        stack_ds += delta_s[d_index]
    stack_dg = stack_dh - 310.15 * stack_ds / 1000
    return [stack_dh, stack_dg]


def intermolecular_initiation_energy() -> list:
    ii_dh = -7.2  # intermolecular initiation dh
    ii_dg = 1.0  # intermolecular initiation dg
    return [ii_dh, ii_dg]


def states_correction_dg(seq: str, bulge_pos: int) -> float:
    state_num = 1
    # get number of states
    bulge_base = seq[bulge_pos]
    forward_idx = bulge_pos
    afterward_idx = bulge_pos
    while forward_idx != 0 and afterward_idx != len(seq) - 1:
        if forward_idx != 0 and seq[forward_idx] == bulge_base:
            state_num += 1
            forward_idx += 1
        else:
            forward_idx = 0
        
        if afterward_idx != len(seq) - 1 and seq[afterward_idx] == bulge_base :
            state_num += 1
            afterward_idx += 1
        else:
            afterward_idx = len(seq) - 1  
    return -0.616 * log(state_num)


def bulge_energy(segment: str, bulge_length: int, seq: str, bulge_pos: int) -> list:
    """计算不同长度bulge的总能量值"""
    
    # bulge loop inition
    bulge_loop_dict = defaultdict(deque)
    bulge_loop_dict["dh"] = deque([18.9, -0.6, -2.3] + [-14.1] * 27)
    bulge_loop_dict["dg"] = deque(
        [
            2.9,
            2.3,
            2.5,
            2.7,
            3.0,
            3.2,
            3.4,
            3.5,
            3.6,
            3.7,
            3.9,
            3.9,
            4.0,
            4.1,
            4.2,
            4.3,
            4.3,
            4.4,
            4.4,
            4.5,
            4.5,
            4.6,
            4.6,
            4.7,
            4.7,
            4.7,
            4.8,
            4.9,
            4.9,
            4.9,
        ]
    )
    
    bulge_loop_dh = bulge_loop_dict["dh"][bulge_length - 1] if bulge_length <= 30 else bulge_loop_dict["dh"][-1]
    bulge_loop_dg = bulge_loop_dict["dg"][bulge_length - 1] if bulge_length <= 30 else bulge_loop_dict["dg"][-1]
    
    # bulge length == 1,extra stack energy(bulge in 2 or -2) and - RT ln(number of states) and no AT panalty
    if bulge_length == 1:
        # extra stack energy 
        ex_stack_dh, ex_stack_dg = stack_energy(segment[0] + segment[-1])
        # states correction
        states_dg = states_correction_dg(seq, bulge_pos)
        # print(f"{ex_stack_dg} is ex_stack_dg")
        bulge_loop_dh += ex_stack_dh
        bulge_loop_dg += ex_stack_dg + states_dg
    else:
        # AT closure
        at_dh, at_dg = ATclosure_energy(segment)
        bulge_loop_dh += at_dh
        bulge_loop_dg += at_dg
        
        # n > 6,bulge_loop_dg += 1.75 RT ln(n/6)
        if bulge_length > 6:
            bulge_loop_dg += 1.75 * 0.616 * log(bulge_length / 6)
    
    return [bulge_loop_dh, bulge_loop_dg]


def index_intl11(upstream_base: str, downstream_base: str, x_base: str, y_base: str) -> list:
    """
            X
           A T
           T A
            Y
    """ 
    external_index = base2int(upstream_base + downstream_base)
    internal_index = base2int(x_base + y_base)
    return [external_index, internal_index]


def index_intl21(upstream_base: str, downstream_base: str, x_base: str, y_base: str, y_neibor_base: str) -> list:
    """            
            X
           A  A
           T  T
            YA
    """
    external_index = base2int(y_neibor_base + upstream_base + downstream_base)
    internal_index = base2int(x_base + y_base)
    return [external_index, internal_index]


def index_intl22(upstream_base: str, downstream_base: str, x_base1: str, y_base1: str, x_base2: str, y_base2: str) -> list:
    """
    	A X1 Y1 A
	    T X2 Y2 T
    """
    external_index = base2int(upstream_base + downstream_base)
    internal_index = base2int(x_base1 + x_base2 + y_base1 + y_base2)
    return [external_index, internal_index]
    
    
def symmetric_int_loop_energy(segment1: str, segment2: str, loop_sum: int, int_loop_energy_dict: defaultdict) -> list:
    """对称的int loop结构直接读取矩阵数据累加"""
    symmetric_int_loop_dh = symmetric_int_loop_dg = 0
    # base init
    upstream_base = segment1[0] 
    downstream_base = segment1[-1]
    x_base = segment1[1]
    y_base = segment2[1]
    
    # 1 * 1
    if loop_sum == 2:
       external_index, internal_index = index_intl11(upstream_base, downstream_base, x_base, y_base)
       symmetric_int_loop_dh += int_loop_energy_dict["11dh"][external_index][internal_index]
       symmetric_int_loop_dg += int_loop_energy_dict["11dg"][external_index][internal_index]
    # 1 * 2
    elif loop_sum == 3:
        y_neibor_base = segment2[2]
        external_index, internal_index = index_intl21(upstream_base, downstream_base, x_base, y_base, y_neibor_base)
        symmetric_int_loop_dh += int_loop_energy_dict["21dh"][external_index][internal_index]
        symmetric_int_loop_dg += int_loop_energy_dict["21dg"][external_index][internal_index]
    # 2 * 2
    else:
        x_base1 = x_base
        x_base2 = y_base
        y_base1 = segment1[2]
        y_base2 = segment2[2]
        external_index, internal_index = index_intl22(upstream_base, downstream_base, x_base1, y_base1, x_base2, y_base2)
        symmetric_int_loop_dh += int_loop_energy_dict["22dh"][external_index][internal_index]
        symmetric_int_loop_dg += int_loop_energy_dict["22dg"][external_index][internal_index]
    print(f"symmetric_int_loop_dg is {symmetric_int_loop_dg}")
    
    return [symmetric_int_loop_dh, symmetric_int_loop_dg]


def asymmetry_correct_energy(loop_diff_abs: int) -> list:
    asymmetry_dg = 0.4
    asymmetry_dh = 0
    asymmetry_correct_dh = asymmetry_dh * loop_diff_abs
    asymmetry_correct_dg = asymmetry_dg * loop_diff_abs
    return [asymmetry_correct_dh, asymmetry_correct_dg]


def asymmetric_int_loop_initiation_energy(loop_sum: int) -> list:
    initiation_dg = [
        3.1,
        3.5,
        3.9,
        4.1,
        4.2,
        4.3,
        4.5,
        4.6,
        4.6,
        4.7,
        4.8,
        4.9,
        5.0,
        5.0,
        5.1,
        5.1,
        5.2,
        5.3,
        5.3,
        5.3,
        5.4,
        5.4,
        5.5,
        5.5,
        5.6,
        5.6,
        5.6,
    ]
    initiation_dh = [0.0] * 27

    return [initiation_dh[loop_sum - 4], initiation_dg[loop_sum - 4]]


def asymmetric_int_loop_mismatch_energy(
    segment1: str, segment2: str, int_loop_type: list
) -> list:
    """只有非对称且min loop > 1 才考虑mismatch"""
    mm_dg = deque([
        # A X
        # T Y
        [
            -0.7,-0.3,-0.5,-0.0,
            -0.6,-0.2,0.0,-0.3,
            -0.6,0.0,-0.4,-0.5,
            0.0,-0.3,-0.5,-0.4,
        ],
        # C X
        # G Y
        [
            -1.0,-0.8,-0.9,0.0,
            -0.8,-0.5,0.0,-0.7,
            -1.0,0.0,-0.9,-1.0,
            0.0,-0.6,-0.9,-0.9,
        ],
        # G X
        # C Y
        [
            -1.0,-0.7,-0.8,0.0,
            -1.0,-0.6,0.0,-0.7,
            -1.0,0.0,-1.0,-0.8,
            0.0,-0.6,-0.9,-0.9,
        ],
        # T X
        # A Y
        [
            -0.6,-0.4,-0.5,0.0,
            -0.5,-0.2,0.0,-0.5,
            -0.6,-0.0,-0.4,-0.5,
            0.0,-0.3,-0.6,-0.3,
        ],
])
    mm_dh = deque([
        # A X
        # T Y
        [
            4.0,5.2,4.7,4.3,
            4.0,5.2,4.7,4.3,
            4.0,5.2,4.7,4.3,
            4.0,5.2,4.7,4.3,
        ],
        # C X
        # G Y
        [
            -4.6,-4.2,-4.7,-5.0,
            -4.6,-4.2,-4.7,-5.0,
            -4.6,-4.2,-4.7,-5.0,
            -4.6,-4.2,-4.7,-5.0,
        ],
        # G X
        # C Y
        [
            -6.2,-3.3,-5.1,-0.9,
            -6.2,-3.3,-5.1,-0.9,
            -6.2,-3.3,-5.1,-0.9,
            -6.2,-3.3,-5.1,-0.9,
        ],
        # T X
        # A Y
        [
            2.3,1.7,0.7,-5.8,
            2.3,1.7,0.7,-5.8,
            2.3,1.7,0.7,-5.8,
            2.3,1.7,0.7,-5.8,
        ],
])

    mismatch_dh = mismatch_dg = 0
    mismatch1_dh = mismatch1_dg = 0
    mismatch2_dh = mismatch2_dg = 0
    # min loop length == 1,do NOT consider mismatch energy
    if int_loop_type[0] == 1:return [0.0, 0.0]
    # mismatch energy calc by segment
    external_index = base2int(segment1[0])
    internal_index = base2int(segment1[1] + segment2[1])
    mismatch1_dh = mm_dh[external_index][internal_index]
    mismatch1_dg = mm_dg[external_index][internal_index]
    
    external_index = base2int(segment2[-1])
    internal_index = base2int(segment2[-2] + segment1[-2])
    mismatch2_dh = mm_dh[external_index][internal_index]
    mismatch2_dg = mm_dg[external_index][internal_index]
    
    mismatch_dh = mismatch1_dh + mismatch2_dh
    mismatch_dg = mismatch1_dg + mismatch2_dg
    
    return [mismatch_dh, mismatch_dg]


def ATclosure_energy(segment1: str) -> list:
    closure_AT_dh = closure_AT_dg = 0
    if segment1[0] == "A" or segment1[0] == "T":
        closure_AT_dh += 3.2
    if segment1[-1] == "A" or segment1[-1] == "T":
        closure_AT_dh += 3.2
    return [closure_AT_dh, closure_AT_dg]
        
        

def asymmetric_int_loop_energy(
    segment1: str, segment2: str, loop_sum: int, loop_diff_abs: int, loop_type: list
) -> list:
    asymmetric_int_loop_dh = asymmetric_int_loop_dg = 0
    # loop sum initiation energy
    li_dh, li_dg = asymmetric_int_loop_initiation_energy(loop_sum)
    # asymmetry correct energy
    asymmetry_correct_dh, asymmetry_correct_dg = asymmetry_correct_energy(loop_diff_abs)
    # mismatch energy
    mm_dh, mm_dg = asymmetric_int_loop_mismatch_energy(segment1, segment2, loop_type)
    # AT closure 
    at_dh, at_dg = ATclosure_energy(segment1)
    asymmetric_int_loop_dh = li_dh + asymmetry_correct_dh + mm_dh + at_dh
    asymmetric_int_loop_dg = li_dg + asymmetry_correct_dg + mm_dg + at_dg
    return [asymmetric_int_loop_dh, asymmetric_int_loop_dg]

def int_loop_energy(segment1: str, segment2: str, int_loop_energy_dict: defaultdict) -> list:
    """对intloop进行分型并计算能量"""
    int_loop_dh = 0
    int_loop_dg = 0


    
    first_loop_length = len(segment1) - segment1.count("-") - 2
    second_loop_length = len(segment2) - segment2.count("-") - 2 
    loop_sum = first_loop_length + second_loop_length
    max_loop_length = max(first_loop_length, second_loop_length)
    min_loop_length = loop_sum - max_loop_length
    loop_type = [min_loop_length, max_loop_length]
    is_symmetric = loop_sum <= 4 and max_loop_length < 3

    # 2 x 1, reverse
    if loop_sum == 3:
        if first_loop_length == 2:
            tmp = segment1[::-1]
            segment1 = segment2[::-1]
            segment2 = tmp
        # ensure - in -2
        if segment1[-2] != "-":
            segment1 = segment1[0] + segment1[-2] + "-" + segment1[-1]
    
    # 1×1, 1×2, 2×2 Internal Loops
    if is_symmetric:
        symmetric_int_loop_dh, symmetric_int_loop_dg = symmetric_int_loop_energy(segment1, segment2, loop_sum, int_loop_energy_dict)
        int_loop_dh += symmetric_int_loop_dh
        int_loop_dg += symmetric_int_loop_dg
    # Other Internal Loops
    else:
        loop_diff_abs = abs(first_loop_length - second_loop_length)
        asymmetric_int_loop_dh, asymmetric_int_loop_dg = asymmetric_int_loop_energy(
            segment1, segment2, loop_sum, loop_diff_abs, loop_type
        )
        int_loop_dh += asymmetric_int_loop_dh
        int_loop_dg += asymmetric_int_loop_dg

    return [int_loop_dh, int_loop_dg]


def symmetry(seq):
    '''
    Return 1 if string is symmetrical, 0 otherwise.
    '''
    seq_len = len(seq)
    mp = seq_len // 2
    if seq_len % 2 == 1:return 0
    for i in range(mp):
        s = seq[i]
        e = seq[seq_len - i - 1]
        terminal_base_set={s,e}
        if terminal_base_set != {"A","T"} and terminal_base_set != {"C","G"}:return 0
    return 1


def divalent_to_monovalent(divalent,dntp):
    '''
    Convert divalent salt concentration to monovalent
    '''
    if divalent == 0:dntp = 0
    if divalent < dntp:divalent = dntp
    return 120 * sqrt((divalent - dntp))


def calc_Tm_by_NN(duplex_str: str, loop_region_dict: defaultdict) -> float:
    # int loop energy matrix
    int_loop_energy_dict = defaultdict(deque)
    int_loop_energy_dict["11dg"] = deque([
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
    
    int_loop_energy_dict["11dh"] = deque([
        # A A
        # T T
        [
            12.9,1.5,10.7,1.5,
            10.4,14.7,10.4,-1.6,
            1.8,1.5,1.6,1.5,
            10.4,3.5,10.4,-4.8,
        ],
        # A C
        # T G
        [
            10.0,10.9,8.9,10.9,
            4.1,12.8,4.1,1.3,
            4.5,10.9,-5.0,10.9,
            4.1,11.6,4.1,0.5,
        ],
        # A G
        # T C
        [
            5.0,3.9,-6.3,3.9,
            5.3,2.6,5.3,-4.9,
            -2.2,3.9,-9.9,3.9,
            5.3,-12.3,5.3,-7.5,
        ],
        # A T
        # T A
        [
            4.9,9.8,7.1,9.8,
            9.8,7.5,9.8,4.3,
            7.1,9.8,4.1,9.8,
            9.8,4.3,9.8,5.7,
        ],
        # C A
        # G T
        [
            10.0,13.2,11.3,13.2,
            9.6,5.9,9.6,0.1,
            6.9,13.2,-1.1,13.2,
            9.6,3.0,9.6,0.7,
        ],
        # C C
        # G G
        [
            -0.6,4.6,1.3,4.6,
            9.3,7.5,9.3,1.9,
            3.9,4.6,-0.9,4.6,
            9.3,0.5,9.3,2.1,
        ],
        # C G
        # G C
        [
            4.0,5.1,6.1,5.1,
            5.1,4.1,5.1,1.6,
            6.1,5.1,0.1,5.1,
            5.1,1.6,5.1,-2.7,
        ],
        # C T
        # G A
        [
            5.0,5.3,-2.2,5.3,
            3.9,2.6,3.9,-12.3,
            -6.3,5.3,-9.9,5.3,
            3.9,-4.9,3.9,-7.5,
        ],
        # G A
        # C T
        [
            -1.4,13.2,11.3,13.2,
            4.2,12.9,4.2,-0.4,
            6.9,13.2,3.4,13.2,
            4.2,3.7,4.2,-10.0,
        ],
        # G C
        # C G
        [
            2.6,3.9,2.4,3.9,
            3.9,6.0,3.9,-0.6,
            2.4,3.9,-2.8,3.9,
            3.9,-0.6,3.9,0.3,
        ],
        # G G
        # C C
        [
            -0.6,9.3,3.9,9.3,
            4.6,7.5,4.6,0.5,
            1.3,9.3,-0.9,9.3,
            4.6,1.9,4.6,2.1,
        ],
        # G T
        # C A
        [
            10.0,4.1,4.5,4.1,
            10.9,12.8,10.9,11.6,
            8.9,4.1,-5.0,4.1,
            10.9,1.3,10.9,0.5,
        ],
        # T A
        # A T
        [
            12.1,14.3,11.9,14.3,
            14.3,17.4,14.3,2.5,
            11.9,14.3,8.4,14.3,
            14.3,2.5,14.3,6.2,
        ],
        # T C
        # A G
        [
            -1.4,4.2,6.9,4.2,
            13.2,12.9,13.2,3.7,
            11.3,4.2,3.4,4.2,
            13.2,-0.4,13.2,-10.0,
        ],
        # T G
        # A C
        [
            10.0,9.6,6.9,9.6,
            13.2,5.9,13.2,3.0,
            11.3,9.6,-1.1,9.6,
            13.2,0.1,13.2,0.7,
        ],
        # T T
        # A A
        [
            12.9,10.4,1.8,10.4,
            1.5,14.7,1.5,3.5,
            10.7,10.4,1.6,10.4,
            1.5,-1.6,1.5,-4.8,
        ],
])
    
    int_loop_energy_dict["21dg"] = deque([
        # A  A
        # T AT
        [
        	1.8,2.0,1.9,2.5,
        	2.1,2.3,2.2,2.2,
        	1.9,2.3,1.9,2.5,
        	3.1,2.7,3.1,2.5,
        ],
        # A  C
        # T AG
        [
        	1.3,1.5,1.4,2.0,
        	1.8,2.0,1.9,1.9,
        	1.5,1.9,1.5,2.1,
        	2.6,2.2,2.6,2.0,
        ],
        # A  G
        # T AC
        [
        	1.4,1.6,1.5,2.1,
        	1.9,2.1,2.0,2.0,
        	1.6,2.0,1.6,2.2,
        	2.7,2.3,2.7,2.1,
        ],
        # A  T
        # T AA
        [
        	1.7,1.9,1.8,2.4,
        	2.1,2.3,2.2,2.2,
        	1.9,2.3,1.9,2.5,
        	3.1,2.7,3.1,2.5,
        ],
        # A  A
        # T CT
        [
        	1.9,2.1,2.0,2.6,
        	2.3,2.5,2.4,2.4,
        	2.0,2.4,2.0,2.6,
        	2.7,2.3,2.7,2.1,
        ],
        # A  C
        # T CG
        [
        	1.6,1.8,1.7,2.3,
        	2.0,2.2,2.1,2.1,
        	1.2,1.6,1.2,1.8,
        	2.4,2.0,2.4,1.8,
        ],
        # A  G
        # T CC
        [
        	1.7,1.9,1.8,2.4,
        	2.1,2.3,2.2,2.2,
        	1.7,2.1,1.7,2.3,
        	2.4,2.0,2.4,1.8,
        ],
        # A  T
        # T CA
        [
        	1.9,2.1,2.0,2.6,
        	2.3,2.5,2.4,2.4,
        	2.1,2.5,2.1,2.7,
        	2.7,2.3,2.7,2.1,
        ],
        # A  A
        # T GT
        [
        	1.8,2.0,1.9,2.5,
        	2.6,2.8,2.7,2.7,
        	1.9,2.3,1.9,2.5,
        	3.1,2.7,3.1,2.5,
        ],
        # A  C
        # T GG
        [
        	1.4,1.6,1.5,2.1,
        	1.8,2.0,1.9,1.9,
        	1.5,1.9,1.5,2.1,
        	2.8,2.4,2.8,2.2,
        ],
        # A  G
        # T GC
        [
        	1.5,1.7,1.6,2.2,
        	1.5,1.7,1.6,1.6,
        	1.6,2.0,1.6,2.2,
        	2.6,2.2,2.6,2.0,
        ],
        # A  T
        # T GA
        [
        	1.8,2.0,1.9,2.5,
        	2.4,2.6,2.5,2.5,
        	1.9,2.3,1.9,2.5,
        	3.1,2.7,3.1,2.5,
        ],
        # A  A
        # T TT
        [
        	2.4,2.6,2.5,3.1,
        	2.2,2.4,2.3,2.3,
        	2.4,2.8,2.4,3.0,
        	2.5,2.1,2.5,1.9,
        ],
        # A  C
        # T TG
        [
        	2.0,2.2,2.1,2.7,
        	1.9,2.1,2.0,2.0,
        	2.1,2.5,2.1,2.7,
        	2.2,1.8,2.2,1.6,
        ],
        # A  G
        # T TC
        [
        	2.2,2.4,2.3,2.9,
        	1.9,2.1,2.0,2.0,
        	2.1,2.5,2.1,2.7,
        	2.2,1.8,2.2,1.6,
        ],
        # A  T
        # T TA
        [
        	2.4,2.6,2.5,3.1,
        	2.2,2.4,2.3,2.3,
        	2.5,2.9,2.5,3.1,
        	2.5,2.1,2.5,1.9,
        ],
        # C  A
        # G AT
        [
        	1.5,1.8,1.6,2.1,
        	1.9,2.1,1.8,1.9,
        	1.6,1.4,1.6,2.0,
        	2.9,2.4,2.7,2.2,
        ],
        # C  C
        # G AG
        [
        	1.0,1.3,1.1,1.6,
        	1.6,1.8,1.5,1.6,
        	1.2,1.0,1.2,1.6,
        	2.4,1.9,2.2,1.7,
        ],
        # C  G
        # G AC
        [
        	1.1,1.4,1.2,1.7,
        	1.7,1.9,1.6,1.7,
        	1.3,1.1,1.3,1.7,
        	2.5,2.0,2.3,1.8,
        ],
        # C  T
        # G AA
        [
        	1.4,1.7,1.5,2.0,
        	1.9,2.1,1.8,1.9,
        	1.6,1.4,1.6,2.0,
        	2.9,2.4,2.7,2.2,
        ],
        # C  A
        # G CT
        [
        	1.6,1.9,1.7,2.2,
        	2.1,2.3,2.0,2.1,
        	1.7,1.5,1.7,2.1,
        	2.5,2.0,2.3,1.8,
        ],
        # C  C
        # G CG
        [
        	1.3,1.6,1.4,1.9,
        	1.8,2.0,1.7,1.8,
        	0.9,0.7,0.9,1.3,
        	2.2,1.7,2.0,1.5,
        ],
        # C  G
        # G CC
        [
        	1.4,1.7,1.5,2.0,
        	1.9,2.1,1.8,1.9,
        	1.4,1.2,1.4,1.8,
        	2.2,1.7,2.0,1.5,
        ],
        # C  T
        # G CA
        [
        	1.6,1.9,1.7,2.2,
        	2.1,2.3,2.0,2.1,
        	1.8,1.6,1.8,2.2,
        	2.5,2.0,2.3,1.8,
        ],
        # C  A
        # G GT
        [
        	1.5,1.8,1.6,2.1,
        	2.4,2.6,2.3,2.4,
        	1.6,1.4,1.6,2.0,
        	2.9,2.4,2.7,2.2,
        ],
        # C  C
        # G GG
        [
        	1.1,1.4,1.2,1.7,
        	1.6,1.8,1.5,1.6,
        	1.2,1.0,1.2,1.6,
        	2.6,2.1,2.4,1.9,
        ],
        # C  G
        # G GC
        [
        	1.2,1.5,1.3,1.8,
        	1.3,1.5,1.2,1.3,
        	1.3,1.1,1.3,1.7,
        	2.4,1.9,2.2,1.7,
        ],
        # C  T
        # G GA
        [
        	1.5,1.8,1.6,2.1,
        	2.2,2.4,2.1,2.2,
        	1.6,1.4,1.6,2.0,
        	2.9,2.4,2.7,2.2,
        ],
        # C  A
        # G TT
        [
        	2.1,2.4,2.2,2.7,
        	2.0,2.2,1.9,2.0,
        	2.1,1.9,2.1,2.5,
        	2.3,1.8,2.1,1.6,
        ],
        # C  C
        # G TG
        [
        	1.7,2.0,1.8,2.3,
        	1.7,1.9,1.6,1.7,
        	1.8,1.6,1.8,2.2,
        	2.0,1.5,1.8,1.3,
        ],
        # C  G
        # G TC
        [
        	1.9,2.2,2.0,2.5,
        	1.7,1.9,1.6,1.7,
        	1.8,1.6,1.8,2.2,
        	2.0,1.5,1.8,1.3,
        ],
        # C  T
        # G TA
        [
        	2.1,2.4,2.2,2.7,
        	2.0,2.2,1.9,2.0,
        	2.2,2.0,2.2,2.6,
        	2.3,1.8,2.1,1.6,
        ],
        # G  A
        # C AT
        [
        	1.4,1.7,1.5,2.0,
        	1.8,2.0,1.3,1.9,
        	1.5,1.7,1.5,2.2,
        	2.7,2.4,2.7,2.2,
        ],
        # G  C
        # C AG
        [
        	0.9,1.2,1.0,1.5,
        	1.5,1.7,1.0,1.6,
        	1.1,1.3,1.1,1.8,
        	2.2,1.9,2.2,1.7,
        ],
        # G  G
        # C AC
        [
        	1.0,1.3,1.1,1.6,
        	1.6,1.8,1.1,1.7,
        	1.2,1.4,1.2,1.9,
        	2.3,2.0,2.3,1.8,
        ],
        # G  T
        # C AA
        [
        	1.3,1.6,1.4,1.9,
        	1.8,2.0,1.3,1.9,
        	1.5,1.7,1.5,2.2,
        	2.7,2.4,2.7,2.2,
        ],
        # G  A
        # C CT
        [
        	1.5,1.8,1.6,2.1,
        	2.0,2.2,1.5,2.1,
        	1.6,1.8,1.6,2.3,
        	2.3,2.0,2.3,1.8,
        ],
        # G  C
        # C CG
        [
        	1.2,1.5,1.3,1.8,
        	1.7,1.9,1.2,1.8,
        	0.8,1.0,0.8,1.5,
        	2.0,1.7,2.0,1.5,
        ],
        # G  G
        # C CC
        [
        	1.3,1.6,1.4,1.9,
        	1.8,2.0,1.3,1.9,
        	1.3,1.5,1.3,2.0,
        	2.0,1.7,2.0,1.5,
        ],
        # G  T
        # C CA
        [
        	1.5,1.8,1.6,2.1,
        	2.0,2.2,1.5,2.1,
        	1.7,1.9,1.7,2.4,
        	2.3,2.0,2.3,1.8,
        ],
        # G  A
        # C GT
        [
        	1.4,1.7,1.5,2.0,
        	2.3,2.5,1.8,2.4,
        	1.5,1.7,1.5,2.2,
        	2.7,2.4,2.7,2.2,
        ],
        # G  C
        # C GG
        [
        	1.0,1.3,1.1,1.6,
        	1.5,1.7,1.0,1.6,
        	1.1,1.3,1.1,1.8,
        	2.4,2.1,2.4,1.9,
        ],
        # G  G
        # C GC
        [
        	1.1,1.4,1.2,1.7,
        	1.2,1.4,0.7,1.3,
        	1.2,1.4,1.2,1.9,
        	2.2,1.9,2.2,1.7,
        ],
        # G  T
        # C GA
        [
        	1.4,1.7,1.5,2.0,
        	2.1,2.3,1.6,2.2,
        	1.5,1.7,1.5,2.2,
        	2.7,2.4,2.7,2.2,
        ],
        # G  A
        # C TT
        [
        	2.0,2.3,2.1,2.6,
        	1.9,2.1,1.4,2.0,
        	2.0,2.2,2.0,2.7,
        	2.1,1.8,2.1,1.6,
        ],
        # G  C
        # C TG
        [
        	1.6,1.9,1.7,2.2,
        	1.6,1.8,1.1,1.7,
        	1.7,1.9,1.7,2.4,
        	1.8,1.5,1.8,1.3,
        ],
        # G  G
        # C TC
        [
        	1.8,2.1,1.9,2.4,
        	1.6,1.8,1.1,1.7,
        	1.7,1.9,1.7,2.4,
        	1.8,1.5,1.8,1.3,
        ],
        # G  T
        # C TA
        [
        	2.0,2.3,2.1,2.6,
        	1.9,2.1,1.4,2.0,
        	2.1,2.3,2.1,2.8,
        	2.1,1.8,2.1,1.6,
        ],
        # T  A
        # A AT
        [
        	1.9,2.0,1.9,2.5,
        	2.1,2.3,2.1,2.2,
        	1.9,2.5,1.9,2.5,
        	3.1,2.7,3.0,2.5,
        ],
        # T  C
        # A AG
        [
        	1.4,1.5,1.4,2.0,
        	1.8,2.0,1.8,1.9,
        	1.5,2.1,1.5,2.1,
        	2.6,2.2,2.5,2.0,
        ],
        # T  G
        # A AC
        [
        	1.5,1.6,1.5,2.1,
        	1.9,2.1,1.9,2.0,
        	1.6,2.2,1.6,2.2,
        	2.7,2.3,2.6,2.1,
        ],
        # T  T
        # A AA
        [
        	1.8,1.9,1.8,2.4,
        	2.1,2.3,2.1,2.2,
        	1.9,2.5,1.9,2.5,
        	3.1,2.7,3.0,2.5,
        ],
        # T  A
        # A CT
        [
        	2.0,2.1,2.0,2.6,
        	2.3,2.5,2.3,2.4,
        	2.0,2.6,2.0,2.6,
        	2.7,2.3,2.6,2.1,
        ],
        # T  C
        # A CG
        [
        	1.7,1.8,1.7,2.3,
        	2.0,2.2,2.0,2.1,
        	1.2,1.8,1.2,1.8,
        	2.4,2.0,2.3,1.8,
        ],
        # T  G
        # A CC
        [
        	1.8,1.9,1.8,2.4,
        	2.1,2.3,2.1,2.2,
        	1.7,2.3,1.7,2.3,
        	2.4,2.0,2.3,1.8,
        ],
        # T  T
        # A CA
        [
        	2.0,2.1,2.0,2.6,
        	2.3,2.5,2.3,2.4,
        	2.1,2.7,2.1,2.7,
        	2.7,2.3,2.6,2.1,
        ],
        # T  A
        # A GT
        [
        	1.9,2.0,1.9,2.5,
        	2.6,2.8,2.6,2.7,
        	1.9,2.5,1.9,2.5,
        	3.1,2.7,3.0,2.5,
        ],
        # T  C
        # A GG
        [
        	1.5,1.6,1.5,2.1,
        	1.8,2.0,1.8,1.9,
        	1.5,2.1,1.5,2.1,
        	2.8,2.4,2.7,2.2,
        ],
        # T  G
        # A GC
        [
        	1.6,1.7,1.6,2.2,
        	1.5,1.7,1.5,1.6,
        	1.6,2.2,1.6,2.2,
        	2.6,2.2,2.5,2.0,
        ],
        # T  T
        # A GA
        [
        	1.9,2.0,1.9,2.5,
        	2.4,2.6,2.4,2.5,
        	1.9,2.5,1.9,2.5,
        	3.1,2.7,3.0,2.5,
        ],
        # T  A
        # A TT
        [
        	2.5,2.6,2.5,3.1,
        	2.2,2.4,2.2,2.3,
        	2.4,3.0,2.4,3.0,
        	2.5,2.1,2.4,1.9,
        ],
        # T  C
        # A TG
        [
        	2.1,2.2,2.1,2.7,
        	1.9,2.1,1.9,2.0,
        	2.1,2.7,2.1,2.7,
        	2.2,1.8,2.1,1.6,
        ],
        # T  G
        # A TC
        [
        	2.3,2.4,2.3,2.9,
        	1.9,2.1,1.9,2.0,
        	2.1,2.7,2.1,2.7,
        	2.2,1.8,2.1,1.6,
        ],
        # T  T
        # A TA
        [
        	2.5,2.6,2.5,3.1,
        	2.2,2.4,2.2,2.3,
        	2.5,3.1,2.5,3.1,
        	2.5,2.1,2.4,1.9,
        ],
])

    int_loop_energy_dict["21dh"] = deque([
        # A  A
        # T AT
        [
        	31.8,31.8,29.6,31.8,
        	29.3,29.3,29.3,29.3,
        	20.7,20.7,20.7,20.7,
        	29.3,29.3,29.3,14.1,
        ],
        # A  C
        # T AG
        [
        	28.9,29.8,27.8,29.8,
        	23.0,23.0,23.0,20.2,
        	23.4,23.4,13.9,23.4,
        	23.0,23.0,23.0,19.4,
        ],
        # A  G
        # T AC
        [
        	23.9,23.9,12.6,23.9,
        	24.2,24.2,24.2,14.0,
        	16.7,16.7,16.7,16.7,
        	24.2,6.6,24.2,11.4,
        ],
        # A  T
        # T AA
        [
        	23.8,28.7,26.0,23.8,
        	28.7,28.7,28.7,23.2,
        	26.0,26.0,26.0,26.0,
        	28.7,23.2,28.7,28.7,
        ],
        # A  A
        # T CT
        [
        	31.8,20.4,29.6,20.4,
        	29.3,33.6,29.3,17.3,
        	20.7,20.4,20.5,20.4,
        	29.3,22.4,29.3,14.1,
        ],
        # A  C
        # T CG
        [
        	29.8,29.8,27.8,29.8,
        	23.0,31.7,23.0,20.2,
        	23.4,29.8,13.9,29.8,
        	23.0,30.5,23.0,19.4,
        ],
        # A  G
        # T CC
        [
        	22.8,22.8,12.6,22.8,
        	24.2,21.5,24.2,14.0,
        	16.7,22.8,9.0,22.8,
        	6.6,6.6,6.6,6.6,
        ],
        # A  T
        # T CA
        [
        	28.7,28.7,26.0,28.7,
        	28.7,26.4,28.7,23.2,
        	26.0,28.7,23.0,28.7,
        	23.2,23.2,23.2,23.2,
        ],
        # A  A
        # T GT
        [
        	29.6,29.6,29.6,29.6,
        	29.3,29.3,29.3,29.3,
        	20.7,20.5,20.5,20.5,
        	29.3,29.3,29.3,14.1,
        ],
        # A  C
        # T GG
        [
        	27.8,27.8,27.8,27.8,
        	23.0,23.0,23.0,20.2,
        	13.9,13.9,13.9,13.9,
        	23.0,23.0,23.0,19.4,
        ],
        # A  G
        # T GC
        [
        	12.6,12.6,12.6,12.6,
        	24.2,24.2,24.2,14.0,
        	16.7,9.0,9.0,9.0,
        	24.2,6.6,24.2,11.4,
        ],
        # A  T
        # T GA
        [
        	26.0,26.0,26.0,26.0,
        	28.7,28.7,28.7,23.2,
        	26.0,23.0,23.0,23.0,
        	28.7,23.2,28.7,28.7,
        ],
        # A  A
        # T TT
        [
        	31.8,20.4,29.6,20.4,
        	17.3,17.3,17.3,17.3,
        	20.7,20.4,20.5,20.4,
        	14.1,14.1,14.1,14.1,
        ],
        # A  C
        # T TG
        [
        	29.8,29.8,27.8,29.8,
        	20.2,20.2,20.2,20.2,
        	23.4,29.8,13.9,29.8,
        	19.4,19.4,19.4,19.4,
        ],
        # A  G
        # T TC
        [
        	22.8,22.8,12.6,22.8,
        	14.0,14.0,14.0,14.0,
        	16.7,22.8,9.0,22.8,
        	11.4,6.6,11.4,11.4,
        ],
        # A  T
        # T TA
        [
        	28.7,28.7,26.0,28.7,
        	23.2,23.2,23.2,23.2,
        	26.0,28.7,23.0,28.7,
        	28.7,23.2,28.7,24.6,
        ],
        # C  A
        # G AT
        [
        	28.9,28.9,30.2,28.9,
        	28.5,28.5,28.5,19.0,
        	25.8,25.8,17.8,25.8,
        	28.5,21.9,28.5,19.6,
        ],
        # C  C
        # G AG
        [
        	18.3,23.5,20.2,23.5,
        	28.2,28.2,28.2,20.8,
        	22.8,22.8,18.0,22.8,
        	28.2,19.4,28.2,21.0,
        ],
        # C  G
        # G AC
        [
        	22.9,22.9,25.0,22.9,
        	24.0,24.0,24.0,20.5,
        	25.0,25.0,25.0,25.0,
        	24.0,20.5,24.0,16.2,
        ],
        # C  T
        # G AA
        [
        	23.9,24.2,16.7,24.2,
        	22.8,22.8,22.8,6.6,
        	12.6,12.6,9.0,12.6,
        	22.8,14.0,22.8,11.4,
        ],
        # C  A
        # G CT
        [
        	28.9,32.1,30.2,32.1,
        	28.5,24.8,28.5,19.0,
        	25.8,32.1,17.8,32.1,
        	21.9,21.9,21.9,19.6,
        ],
        # C  C
        # G CG
        [
        	23.5,23.5,20.2,23.5,
        	28.2,26.4,28.2,20.8,
        	22.8,23.5,18.0,23.5,
        	19.4,19.4,19.4,19.4,
        ],
        # C  G
        # G CC
        [
        	22.9,24.0,25.0,24.0,
        	24.0,23.0,24.0,20.5,
        	25.0,24.0,19.0,24.0,
        	20.5,20.5,20.5,16.2,
        ],
        # C  T
        # G CA
        [
        	24.2,24.2,16.7,24.2,
        	22.8,21.5,22.8,6.6,
        	12.6,24.2,9.0,24.2,
        	14.0,14.0,14.0,11.4,
        ],
        # C  A
        # G GT
        [
        	30.2,30.2,30.2,30.2,
        	28.5,28.5,28.5,19.0,
        	17.8,17.8,17.8,17.8,
        	28.5,21.9,28.5,19.6,
        ],
        # C  C
        # G GG
        [
        	20.2,20.2,20.2,20.2,
        	28.2,28.2,28.2,20.8,
        	18.0,18.0,18.0,18.0,
        	28.2,19.4,28.2,21.0,
        ],
        # C  G
        # G GC
        [
        	25.0,25.0,25.0,25.0,
        	24.0,24.0,24.0,20.5,
        	25.0,19.0,19.0,19.0,
        	24.0,20.5,24.0,16.2,
        ],
        # C  T
        # G GA
        [
        	16.7,16.7,16.7,16.7,
        	22.8,22.8,22.8,6.6,
        	9.0,9.0,9.0,9.0,
        	22.8,14.0,22.8,11.4,
        ],
        # C  A
        # G TT
        [
        	28.9,32.1,30.2,32.1,
        	19.0,19.0,19.0,19.0,
        	25.8,32.1,17.8,32.1,
        	19.6,19.6,19.6,19.6,
        ],
        # C  C
        # G TG
        [
        	23.5,23.5,20.2,23.5,
        	20.8,20.8,20.8,20.8,
        	22.8,23.5,18.0,23.5,
        	21.0,21.0,21.0,21.0,
        ],
        # C  G
        # G TC
        [
        	22.9,24.0,25.0,24.0,
        	20.5,20.5,20.5,20.5,
        	25.0,24.0,19.0,24.0,
        	16.2,16.2,16.2,16.2,
        ],
        # C  T
        # G TA
        [
        	24.2,24.2,16.7,24.2,
        	6.6,6.6,6.6,6.6,
        	12.6,24.2,9.0,24.2,
        	11.4,11.4,11.4,11.4,
        ],
        # G  A
        # C AT
        [
        	17.5,17.5,17.5,17.5,
        	23.1,23.1,23.1,23.1,
        	25.8,25.8,25.8,25.8,
        	23.1,23.1,23.1,8.9,
        ],
        # G  C
        # C AG
        [
        	21.5,21.5,21.3,21.5,
        	22.8,22.8,22.8,22.8,
        	21.3,21.3,16.1,21.3,
        	22.8,22.8,22.8,22.8,
        ],
        # G  G
        # C AC
        [
        	18.3,18.3,22.8,18.3,
        	23.5,23.5,23.5,23.5,
        	20.2,20.2,20.2,20.2,
        	23.5,23.5,23.5,23.5,
        ],
        # G  T
        # C AA
        [
        	28.9,23.0,23.4,23.0,
        	29.8,29.8,29.8,29.8,
        	27.8,27.8,13.9,27.8,
        	29.8,20.2,29.8,19.4,
        ],
        # G  A
        # C CT
        [
        	17.5,32.1,30.2,32.1,
        	23.1,31.8,23.1,18.5,
        	25.8,32.1,22.3,32.1,
        	23.1,22.6,23.1,8.9,
        ],
        # G  C
        # C CG
        [
        	21.5,22.8,21.3,22.8,
        	22.8,24.9,22.8,18.3,
        	21.3,22.8,16.1,22.8,
        	22.8,18.3,22.8,18.3,
        ],
        # G  G
        # C CC
        [
        	18.3,28.2,22.8,28.2,
        	23.5,26.4,23.5,19.4,
        	20.2,28.2,18.0,28.2,
        	23.5,20.8,23.5,21.0,
        ],
        # G  T
        # C CA
        [
        	23.0,23.0,23.4,23.0,
        	29.8,31.7,29.8,30.5,
        	27.8,23.0,13.9,23.0,
        	20.2,20.2,20.2,20.2,
        ],
        # G  A
        # C GT
        [
        	30.2,30.2,30.2,30.2,
        	23.1,23.1,23.1,23.1,
        	25.8,22.3,22.3,22.3,
        	23.1,23.1,23.1,8.9,
        ],
        # G  C
        # C GG
        [
        	21.3,21.3,21.3,21.3,
        	22.8,22.8,22.8,22.8,
        	16.1,16.1,16.1,16.1,
        	22.8,22.8,22.8,22.8,
        ],
        # G  G
        # C GC
        [
        	22.8,22.8,22.8,22.8,
        	23.5,23.5,23.5,23.5,
        	20.2,18.0,18.0,18.0,
        	23.5,23.5,23.5,23.5,
        ],
        # G  T
        # C GA
        [
        	23.4,23.4,23.4,23.4,
        	29.8,29.8,29.8,29.8,
        	13.9,13.9,13.9,13.9,
        	29.8,20.2,29.8,19.4,
        ],
        # G  A
        # C TT
        [
        	17.5,32.1,30.2,32.1,
        	23.1,18.5,23.1,18.5,
        	25.8,32.1,22.3,32.1,
        	8.9,8.9,8.9,8.9,
        ],
        # G  C
        # C TG
        [
        	21.5,22.8,21.3,22.8,
        	22.8,18.3,22.8,18.3,
        	21.3,22.8,16.1,22.8,
        	22.8,18.3,22.8,19.2,
        ],
        # G  G
        # C TC
        [
        	18.3,28.2,22.8,28.2,
        	23.5,19.4,23.5,19.4,
        	20.2,28.2,18.0,28.2,
        	23.5,21.0,23.5,21.0,
        ],
        # G  T
        # C TA
        [
        	23.0,23.0,23.4,23.0,
        	29.8,30.5,29.8,30.5,
        	27.8,23.0,13.9,23.0,
        	19.4,19.4,19.4,19.4,
        ],
        # T  A
        # A AT
        [
        	31.0,31.0,30.8,31.0,
        	33.2,33.2,33.2,21.4,
        	30.8,30.8,30.8,30.8,
        	33.2,21.4,33.2,25.1,
        ],
        # T  C
        # A AG
        [
        	17.5,17.5,25.8,17.5,
        	32.1,32.1,32.1,32.1,
        	30.2,30.2,22.3,30.2,
        	32.1,18.5,32.1,8.9,
        ],
        # T  G
        # A AC
        [
        	28.9,28.9,25.8,28.9,
        	32.1,32.1,32.1,21.9,
        	30.2,30.2,17.8,30.2,
        	32.1,19.0,32.1,19.6,
        ],
        # T  T
        # A AA
        [
        	31.8,29.3,20.7,29.3,
        	20.4,20.4,20.4,22.4,
        	29.6,29.6,20.5,29.6,
        	20.4,17.3,20.4,14.1,
        ],
        # T  A
        # A CT
        [
        	31.0,33.2,30.8,33.2,
        	33.2,36.3,33.2,21.4,
        	30.8,33.2,27.3,33.2,
        	21.4,21.4,21.4,21.4,
        ],
        # T  C
        # A CG
        [
        	17.5,23.1,25.8,23.1,
        	32.1,31.8,32.1,31.8,
        	30.2,23.1,22.3,23.1,
        	18.5,18.5,18.5,8.9,
        ],
        # T  G
        # A CC
        [
        	28.5,28.5,25.8,28.5,
        	32.1,24.8,32.1,21.9,
        	30.2,28.5,17.8,28.5,
        	19.0,19.0,19.0,19.6,
        ],
        # T  T
        # A CA
        [
        	29.3,29.3,20.7,29.3,
        	20.4,33.6,20.4,22.4,
        	29.6,29.3,20.5,29.3,
        	17.3,17.3,17.3,14.1,
        ],
        # T  A
        # A GT
        [
        	30.8,30.8,30.8,30.8,
        	33.2,33.2,33.2,21.4,
        	30.8,27.3,27.3,27.3,
        	33.2,21.4,33.2,25.1,
        ],
        # T  C
        # A GG
        [
        	25.8,25.8,25.8,25.8,
        	32.1,32.1,32.1,32.1,
        	22.3,22.3,22.3,22.3,
        	32.1,18.5,32.1,8.9,
        ],
        # T  G
        # A GC
        [
        	25.8,25.8,25.8,25.8,
        	32.1,32.1,32.1,21.9,
        	17.8,17.8,17.8,17.8,
        	32.1,19.0,32.1,19.6,
        ],
        # T  T
        # A GA
        [
        	20.7,20.7,20.7,20.7,
        	20.4,20.4,20.4,22.4,
        	20.5,20.5,20.5,20.5,
        	20.4,17.3,20.4,14.1,
        ],
        # T  A
        # A TT
        [
        	31.0,33.2,30.8,33.2,
        	21.4,21.4,21.4,21.4,
        	30.8,33.2,27.3,33.2,
        	25.1,21.4,25.1,25.1,
        ],
        # T  C
        # A TG
        [
        	17.5,23.1,25.8,23.1,
        	32.1,22.6,32.1,22.6,
        	30.2,23.1,22.3,23.1,
        	8.9,8.9,8.9,8.9,
        ],
        # T  G
        # A TC
        [
        	28.5,28.5,25.8,28.5,
        	21.9,21.9,21.9,21.9,
        	30.2,28.5,17.8,28.5,
        	19.6,19.6,19.6,19.6,
        ],
        # T  T
        # A TA
        [
        	29.3,29.3,20.7,29.3,
        	22.4,22.4,22.4,22.4,
        	29.6,29.3,20.5,29.3,
        	14.1,14.1,14.1,14.1,
        ],
])

    int_loop_energy_dict["22dg"] = deque([
        # A A
        # T T
        [
        	1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	2.2,2.3,2.2,2.8,2.3,2.5,2.8,2.4,2.2,2.3,2.2,2.7,2.8,2.4,2.8,2.2,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	2.3,2.4,2.3,2.9,2.4,2.6,2.9,2.5,2.3,2.4,2.3,2.8,2.9,2.5,2.9,2.3,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        ],
        # A C
        # T G
        [
        	1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.7,2.0,1.8,2.4,2.0,2.2,2.0,2.1,1.8,1.5,1.8,2.4,2.3,2.1,2.5,1.9,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	1.8,2.1,1.9,2.5,2.1,2.3,2.1,2.2,1.9,1.6,1.9,2.5,2.4,2.2,2.6,2.0,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        ],
        # A G
        # T C
        [
        	1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.8,2.1,1.9,2.6,2.1,2.3,1.7,2.1,1.9,2.0,1.9,2.4,2.4,2.1,2.3,1.9,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	1.9,2.2,2.0,2.7,2.2,2.4,1.8,2.2,2.0,2.1,2.0,2.5,2.5,2.2,2.4,2.0,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        ],
        # A T
        # T A
        [
        	1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	2.1,2.3,2.2,2.8,2.3,2.5,2.6,2.4,2.2,2.4,2.2,2.8,2.8,2.4,2.8,2.2,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	2.2,2.4,2.3,2.9,2.4,2.6,2.7,2.5,2.3,2.5,2.3,2.9,2.9,2.5,2.9,2.3,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        ],
        # C A
        # G T
        [
        	1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
        	1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
        	1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
        	1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
        	1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
        	1.4,1.5,1.4,2.0,1.5,1.7,2.0,1.6,1.4,1.5,1.4,1.9,2.0,1.6,2.0,1.4,
        	1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	2.3,2.4,2.3,2.9,2.4,2.6,2.9,2.5,2.3,2.4,2.3,2.8,2.9,2.5,2.9,2.3,
        	1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
        ],
        # C C
        # G G
        [
        	1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
        	1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
        	1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
        	1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
        	1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
        	0.9,1.2,1.0,1.6,1.2,1.4,1.2,1.3,1.0,0.7,1.0,1.6,1.5,1.3,1.7,1.1,
        	1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.8,2.1,1.9,2.5,2.1,2.3,2.1,2.2,1.9,1.6,1.9,2.5,2.4,2.2,2.6,2.0,
        	1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
        ],
        # C G
        # G C
        [
        	1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
        	1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
        	1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
        	1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
        	1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
        	1.0,1.3,1.1,1.8,1.3,1.5,0.9,1.3,1.1,1.2,1.1,1.6,1.6,1.3,1.5,1.1,
        	1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.9,2.2,2.0,2.7,2.2,2.4,1.8,2.2,2.0,2.1,2.0,2.5,2.5,2.2,2.4,2.0,
        	1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
        ],
        # C T
        # G A
        [
        	1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
        	1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
        	1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
        	1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
        	1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
        	1.3,1.5,1.4,2.0,1.5,1.7,1.8,1.6,1.4,1.6,1.4,2.0,2.0,1.6,2.0,1.4,
        	1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	2.2,2.4,2.3,2.9,2.4,2.6,2.7,2.5,2.3,2.5,2.3,2.9,2.9,2.5,2.9,2.3,
        	1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
        ],
        # G A
        # C T
        [
        	1.4,1.5,1.4,2.0,1.5,1.7,2.0,1.6,1.4,1.5,1.4,1.9,2.0,1.6,2.0,1.4,
        	1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
        	1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	1.2,1.3,1.2,1.8,1.3,1.5,1.8,1.4,1.2,1.3,1.2,1.7,1.8,1.4,1.8,1.2,
        	1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
        	1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
        	1.7,1.8,1.7,2.3,1.8,2.0,2.3,1.9,1.7,1.8,1.7,2.2,2.3,1.9,2.3,1.7,
        	1.5,1.6,1.5,2.1,1.6,1.8,2.1,1.7,1.5,1.6,1.5,2.0,2.1,1.7,2.1,1.5,
        	2.2,2.3,2.2,2.8,2.3,2.5,2.8,2.4,2.2,2.3,2.2,2.7,2.8,2.4,2.8,2.2,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	1.8,1.9,1.8,2.4,1.9,2.1,2.4,2.0,1.8,1.9,1.8,2.3,2.4,2.0,2.4,1.8,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	1.6,1.7,1.6,2.2,1.7,1.9,2.2,1.8,1.6,1.7,1.6,2.1,2.2,1.8,2.2,1.6,
        ],
        # G C
        # C G
        [
        	0.9,1.2,1.0,1.6,1.2,1.4,1.2,1.3,1.0,0.7,1.0,1.6,1.5,1.3,1.7,1.1,
        	1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
        	1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	0.7,1.0,0.8,1.4,1.0,1.2,1.0,1.1,0.8,0.5,0.8,1.4,1.3,1.1,1.5,0.9,
        	1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
        	1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
        	1.2,1.5,1.3,1.9,1.5,1.7,1.5,1.6,1.3,1.0,1.3,1.9,1.8,1.6,2.0,1.4,
        	1.0,1.3,1.1,1.7,1.3,1.5,1.3,1.4,1.1,0.8,1.1,1.7,1.6,1.4,1.8,1.2,
        	1.7,2.0,1.8,2.4,2.0,2.2,2.0,2.1,1.8,1.5,1.8,2.4,2.3,2.1,2.5,1.9,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.3,1.6,1.4,2.0,1.6,1.8,1.6,1.7,1.4,1.1,1.4,2.0,1.9,1.7,2.1,1.5,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.1,1.4,1.2,1.8,1.4,1.6,1.4,1.5,1.2,0.9,1.2,1.8,1.7,1.5,1.9,1.3,
        ],
        # G G
        # C C
        [
        	1.0,1.3,1.1,1.8,1.3,1.5,0.9,1.3,1.1,1.2,1.1,1.6,1.6,1.3,1.5,1.1,
        	1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
        	1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	0.8,1.1,0.9,1.6,1.1,1.3,0.7,1.1,0.9,1.0,0.9,1.4,1.4,1.1,1.3,0.9,
        	1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
        	1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
        	1.3,1.6,1.4,2.1,1.6,1.8,1.2,1.6,1.4,1.5,1.4,1.9,1.9,1.6,1.8,1.4,
        	1.1,1.4,1.2,1.9,1.4,1.6,1.0,1.4,1.2,1.3,1.2,1.7,1.7,1.4,1.6,1.2,
        	1.8,2.1,1.9,2.6,2.1,2.3,1.7,2.1,1.9,2.0,1.9,2.4,2.4,2.1,2.3,1.9,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	1.4,1.7,1.5,2.2,1.7,1.9,1.3,1.7,1.5,1.6,1.5,2.0,2.0,1.7,1.9,1.5,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	1.2,1.5,1.3,2.0,1.5,1.7,1.1,1.5,1.3,1.4,1.3,1.8,1.8,1.5,1.7,1.3,
        ],
        # G T
        # C A
        [
        	1.3,1.5,1.4,2.0,1.5,1.7,1.8,1.6,1.4,1.6,1.4,2.0,2.0,1.6,2.0,1.4,
        	1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
        	1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	1.1,1.3,1.2,1.8,1.3,1.5,1.6,1.4,1.2,1.4,1.2,1.8,1.8,1.4,1.8,1.2,
        	1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
        	1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
        	1.6,1.8,1.7,2.3,1.8,2.0,2.1,1.9,1.7,1.9,1.7,2.3,2.3,1.9,2.3,1.7,
        	1.4,1.6,1.5,2.1,1.6,1.8,1.9,1.7,1.5,1.7,1.5,2.1,2.1,1.7,2.1,1.5,
        	2.1,2.3,2.2,2.8,2.3,2.5,2.6,2.4,2.2,2.4,2.2,2.8,2.8,2.4,2.8,2.2,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	1.7,1.9,1.8,2.4,1.9,2.1,2.2,2.0,1.8,2.0,1.8,2.4,2.4,2.0,2.4,1.8,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	1.5,1.7,1.6,2.2,1.7,1.9,2.0,1.8,1.6,1.8,1.6,2.2,2.2,1.8,2.2,1.6,
        ],
        # T A
        # A T
        [
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	2.2,2.3,2.2,2.8,2.3,2.5,2.8,2.4,2.2,2.3,2.2,2.7,2.8,2.4,2.8,2.2,
        	2.0,2.1,2.0,2.6,2.1,2.3,2.6,2.2,2.0,2.1,2.0,2.5,2.6,2.2,2.6,2.0,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	2.5,2.6,2.5,3.1,2.6,2.8,3.1,2.7,2.5,2.6,2.5,3.0,3.1,2.7,3.1,2.5,
        	2.1,2.2,2.1,2.7,2.2,2.4,2.7,2.3,2.1,2.2,2.1,2.6,2.7,2.3,2.7,2.1,
        	2.4,2.5,2.4,3.0,2.5,2.7,3.0,2.6,2.4,2.5,2.4,2.9,3.0,2.6,3.0,2.4,
        	1.9,2.0,1.9,2.5,2.0,2.2,2.5,2.1,1.9,2.0,1.9,2.4,2.5,2.1,2.5,1.9,
        ],
        # T C
        # A G
        [
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.7,2.0,1.8,2.4,2.0,2.2,2.0,2.1,1.8,1.5,1.8,2.4,2.3,2.1,2.5,1.9,
        	1.5,1.8,1.6,2.2,1.8,2.0,1.8,1.9,1.6,1.3,1.6,2.2,2.1,1.9,2.3,1.7,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	2.0,2.3,2.1,2.7,2.3,2.5,2.3,2.4,2.1,1.8,2.1,2.7,2.6,2.4,2.8,2.2,
        	1.6,1.9,1.7,2.3,1.9,2.1,1.9,2.0,1.7,1.4,1.7,2.3,2.2,2.0,2.4,1.8,
        	1.9,2.2,2.0,2.6,2.2,2.4,2.2,2.3,2.0,1.7,2.0,2.6,2.5,2.3,2.7,2.1,
        	1.4,1.7,1.5,2.1,1.7,1.9,1.7,1.8,1.5,1.2,1.5,2.1,2.0,1.8,2.2,1.6,
        ],
        # T G
        # A C
        [
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.8,2.1,1.9,2.6,2.1,2.3,1.7,2.1,1.9,2.0,1.9,2.4,2.4,2.1,2.3,1.9,
        	1.6,1.9,1.7,2.4,1.9,2.1,1.5,1.9,1.7,1.8,1.7,2.2,2.2,1.9,2.1,1.7,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	2.1,2.4,2.2,2.9,2.4,2.6,2.0,2.4,2.2,2.3,2.2,2.7,2.7,2.4,2.6,2.2,
        	1.7,2.0,1.8,2.5,2.0,2.2,1.6,2.0,1.8,1.9,1.8,2.3,2.3,2.0,2.2,1.8,
        	2.0,2.3,2.1,2.8,2.3,2.5,1.9,2.3,2.1,2.2,2.1,2.6,2.6,2.3,2.5,2.1,
        	1.5,1.8,1.6,2.3,1.8,2.0,1.4,1.8,1.6,1.7,1.6,2.1,2.1,1.8,2.0,1.6,
        ],
        # T T
        # A A
        [
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	2.1,2.3,2.2,2.8,2.3,2.5,2.6,2.4,2.2,2.4,2.2,2.8,2.8,2.4,2.8,2.2,
        	1.9,2.1,2.0,2.6,2.1,2.3,2.4,2.2,2.0,2.2,2.0,2.6,2.6,2.2,2.6,2.0,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	2.4,2.6,2.5,3.1,2.6,2.8,2.9,2.7,2.5,2.7,2.5,3.1,3.1,2.7,3.1,2.5,
        	2.0,2.2,2.1,2.7,2.2,2.4,2.5,2.3,2.1,2.3,2.1,2.7,2.7,2.3,2.7,2.1,
        	2.3,2.5,2.4,3.0,2.5,2.7,2.8,2.6,2.4,2.6,2.4,3.0,3.0,2.6,3.0,2.4,
        	1.8,2.0,1.9,2.5,2.0,2.2,2.3,2.1,1.9,2.1,1.9,2.5,2.5,2.1,2.5,1.9,
        ],
])

    int_loop_energy_dict["22dh"] = deque([
        # A A
        # T T
        [
        	6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
        	7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
        	7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
        	6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
        	6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
        	7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
        	7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
        	6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
        	6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
        	7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
        	7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
        	6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
        	6.3,6.3,6.3,6.3,6.6,6.6,6.6,6.6,4.7,4.7,4.7,4.7,-1.8,-1.8,-1.8,-1.8,
        	7.5,7.5,7.5,7.5,7.8,7.8,7.8,7.8,5.9,5.9,5.9,5.9,-0.6,-0.6,-0.6,-0.6,
        	7.0,7.0,7.0,7.0,7.3,7.3,7.3,7.3,5.4,5.4,5.4,5.4,-1.1,-1.1,-1.1,-1.1,
        	6.6,6.6,6.6,6.6,6.9,6.9,6.9,6.9,5.0,5.0,5.0,5.0,-1.5,-1.5,-1.5,-1.5,
        ],
        # A C
        # T G
        [
        	-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
        	-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
        	-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
        	-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
        	-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
        	-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
        	-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
        	-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
        	-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
        	-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
        	-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
        	-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
        	-2.2,-2.2,-2.2,-2.2,0.7,0.7,0.7,0.7,-1.1,-1.1,-1.1,-1.1,3.1,3.1,3.1,3.1,
        	-1.0,-1.0,-1.0,-1.0,1.9,1.9,1.9,1.9,0.1,0.1,0.1,0.1,4.3,4.3,4.3,4.3,
        	-1.5,-1.5,-1.5,-1.5,1.4,1.4,1.4,1.4,-0.4,-0.4,-0.4,-0.4,3.8,3.8,3.8,3.8,
        	-1.9,-1.9,-1.9,-1.9,1.0,1.0,1.0,1.0,-0.8,-0.8,-0.8,-0.8,3.4,3.4,3.4,3.4,
        ],
        # A G
        # T C
        [
        	-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
        	0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
        	0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
        	-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
        	-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
        	0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
        	0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
        	-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
        	-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
        	0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
        	0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
        	-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
        	-0.6,-0.6,-0.6,-0.6,-0.2,-0.2,-0.2,-0.2,-0.7,-0.7,-0.7,-0.7,-1.0,-1.0,-1.0,-1.0,
        	0.6,0.6,0.6,0.6,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.2,0.2,0.2,0.2,
        	0.1,0.1,0.1,0.1,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.3,-0.3,-0.3,-0.3,
        	-0.3,-0.3,-0.3,-0.3,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.7,-0.7,-0.7,-0.7,
        ],
        # A T
        # T A
        [
        	8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
        	9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
        	8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
        	8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
        	8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
        	9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
        	8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
        	8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
        	8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
        	9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
        	8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
        	8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
        	8.0,8.0,8.0,8.0,9.2,9.2,9.2,9.2,8.7,8.7,8.7,8.7,8.3,8.3,8.3,8.3,
        	9.2,9.2,9.2,9.2,10.4,10.4,10.4,10.4,9.9,9.9,9.9,9.9,9.5,9.5,9.5,9.5,
        	8.7,8.7,8.7,8.7,9.9,9.9,9.9,9.9,9.4,9.4,9.4,9.4,9.0,9.0,9.0,9.0,
        	8.3,8.3,8.3,8.3,9.5,9.5,9.5,9.5,9.0,9.0,9.0,9.0,8.6,8.6,8.6,8.6,
        ],
        # C A
        # G T
        [
        	-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
        	-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
        	-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
        	-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
        	-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
        	-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
        	-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
        	-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
        	-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
        	-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
        	-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
        	-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
        	-2.3,-2.3,-2.3,-2.3,-2.0,-2.0,-2.0,-2.0,-3.9,-3.9,-3.9,-3.9,-10.4,-10.4,-10.4,-10.4,
        	-1.9,-1.9,-1.9,-1.9,-1.6,-1.6,-1.6,-1.6,-3.5,-3.5,-3.5,-3.5,-10.0,-10.0,-10.0,-10.0,
        	-2.4,-2.4,-2.4,-2.4,-2.1,-2.1,-2.1,-2.1,-4.0,-4.0,-4.0,-4.0,-10.5,-10.5,-10.5,-10.5,
        	-2.7,-2.7,-2.7,-2.7,-2.4,-2.4,-2.4,-2.4,-4.3,-4.3,-4.3,-4.3,-10.8,-10.8,-10.8,-10.8,
        ],
        # C C
        # G G
        [
        	-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
        	-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
        	-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
        	-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
        	-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
        	-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
        	-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
        	-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
        	-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
        	-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
        	-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
        	-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
        	-10.8,-10.8,-10.8,-10.8,-7.9,-7.9,-7.9,-7.9,-9.7,-9.7,-9.7,-9.7,-5.5,-5.5,-5.5,-5.5,
        	-10.4,-10.4,-10.4,-10.4,-7.5,-7.5,-7.5,-7.5,-9.3,-9.3,-9.3,-9.3,-5.1,-5.1,-5.1,-5.1,
        	-10.9,-10.9,-10.9,-10.9,-8.0,-8.0,-8.0,-8.0,-9.8,-9.8,-9.8,-9.8,-5.6,-5.6,-5.6,-5.6,
        	-11.2,-11.2,-11.2,-11.2,-8.3,-8.3,-8.3,-8.3,-10.1,-10.1,-10.1,-10.1,-5.9,-5.9,-5.9,-5.9,
        ],
        # C G
        # G C
        [
        	-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
        	-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
        	-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
        	-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
        	-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
        	-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
        	-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
        	-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
        	-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
        	-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
        	-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
        	-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
        	-9.2,-9.2,-9.2,-9.2,-8.8,-8.8,-8.8,-8.8,-9.3,-9.3,-9.3,-9.3,-9.6,-9.6,-9.6,-9.6,
        	-8.8,-8.8,-8.8,-8.8,-8.4,-8.4,-8.4,-8.4,-8.9,-8.9,-8.9,-8.9,-9.2,-9.2,-9.2,-9.2,
        	-9.3,-9.3,-9.3,-9.3,-8.9,-8.9,-8.9,-8.9,-9.4,-9.4,-9.4,-9.4,-9.7,-9.7,-9.7,-9.7,
        	-9.6,-9.6,-9.6,-9.6,-9.2,-9.2,-9.2,-9.2,-9.7,-9.7,-9.7,-9.7,-10.0,-10.0,-10.0,-10.0,
        ],
        # C T
        # G A
        [
        	-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
        	-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
        	-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
        	-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
        	-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
        	-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
        	-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
        	-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
        	-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
        	-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
        	-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
        	-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
        	-0.6,-0.6,-0.6,-0.6,0.6,0.6,0.6,0.6,0.1,0.1,0.1,0.1,-0.3,-0.3,-0.3,-0.3,
        	-0.2,-0.2,-0.2,-0.2,1.0,1.0,1.0,1.0,0.5,0.5,0.5,0.5,0.1,0.1,0.1,0.1,
        	-0.7,-0.7,-0.7,-0.7,0.5,0.5,0.5,0.5,0.0,0.0,0.0,0.0,-0.4,-0.4,-0.4,-0.4,
        	-1.0,-1.0,-1.0,-1.0,0.2,0.2,0.2,0.2,-0.3,-0.3,-0.3,-0.3,-0.7,-0.7,-0.7,-0.7,
        ],
        # G A
        # C T
        [
        	-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
        	-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
        	-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
        	1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
        	-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
        	-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
        	-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
        	1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
        	-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
        	-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
        	-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
        	1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
        	-3.9,-3.9,-3.9,-3.9,-3.6,-3.6,-3.6,-3.6,-5.5,-5.5,-5.5,-5.5,-12.0,-12.0,-12.0,-12.0,
        	-1.0,-1.0,-1.0,-1.0,-0.7,-0.7,-0.7,-0.7,-2.6,-2.6,-2.6,-2.6,-9.1,-9.1,-9.1,-9.1,
        	-2.8,-2.8,-2.8,-2.8,-2.5,-2.5,-2.5,-2.5,-4.4,-4.4,-4.4,-4.4,-10.9,-10.9,-10.9,-10.9,
        	1.4,1.4,1.4,1.4,1.7,1.7,1.7,1.7,-0.2,-0.2,-0.2,-0.2,-6.7,-6.7,-6.7,-6.7,
        ],
        # G C
        # C G
        [
        	-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
        	-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
        	-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
        	-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
        	-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
        	-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
        	-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
        	-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
        	-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
        	-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
        	-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
        	-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
        	-12.4,-12.4,-12.4,-12.4,-9.5,-9.5,-9.5,-9.5,-11.3,-11.3,-11.3,-11.3,-7.1,-7.1,-7.1,-7.1,
        	-9.5,-9.5,-9.5,-9.5,-6.6,-6.6,-6.6,-6.6,-8.4,-8.4,-8.4,-8.4,-4.2,-4.2,-4.2,-4.2,
        	-11.3,-11.3,-11.3,-11.3,-8.4,-8.4,-8.4,-8.4,-10.2,-10.2,-10.2,-10.2,-6.0,-6.0,-6.0,-6.0,
        	-7.1,-7.1,-7.1,-7.1,-4.2,-4.2,-4.2,-4.2,-6.0,-6.0,-6.0,-6.0,-1.8,-1.8,-1.8,-1.8,
        ],
        # G G
        # C C
        [
        	-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
        	-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
        	-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
        	-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
        	-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
        	-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
        	-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
        	-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
        	-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
        	-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
        	-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
        	-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
        	-10.8,-10.8,-10.8,-10.8,-10.4,-10.4,-10.4,-10.4,-10.9,-10.9,-10.9,-10.9,-11.2,-11.2,-11.2,-11.2,
        	-7.9,-7.9,-7.9,-7.9,-7.5,-7.5,-7.5,-7.5,-8.0,-8.0,-8.0,-8.0,-8.3,-8.3,-8.3,-8.3,
        	-9.7,-9.7,-9.7,-9.7,-9.3,-9.3,-9.3,-9.3,-9.8,-9.8,-9.8,-9.8,-10.1,-10.1,-10.1,-10.1,
        	-5.5,-5.5,-5.5,-5.5,-5.1,-5.1,-5.1,-5.1,-5.6,-5.6,-5.6,-5.6,-5.9,-5.9,-5.9,-5.9,
        ],
        # G T
        # C A
        [
        	-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
        	0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
        	-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
        	3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
        	-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
        	0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
        	-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
        	3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
        	-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
        	0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
        	-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
        	3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
        	-2.2,-2.2,-2.2,-2.2,-1.0,-1.0,-1.0,-1.0,-1.5,-1.5,-1.5,-1.5,-1.9,-1.9,-1.9,-1.9,
        	0.7,0.7,0.7,0.7,1.9,1.9,1.9,1.9,1.4,1.4,1.4,1.4,1.0,1.0,1.0,1.0,
        	-1.1,-1.1,-1.1,-1.1,0.1,0.1,0.1,0.1,-0.4,-0.4,-0.4,-0.4,-0.8,-0.8,-0.8,-0.8,
        	3.1,3.1,3.1,3.1,4.3,4.3,4.3,4.3,3.8,3.8,3.8,3.8,3.4,3.4,3.4,3.4,
        ],
        # T A
        # A T
        [
        	4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
        	4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
        	3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
        	-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
        	4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
        	4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
        	3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
        	-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
        	4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
        	4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
        	3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
        	-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
        	4.6,4.6,4.6,4.6,4.9,4.9,4.9,4.9,3.0,3.0,3.0,3.0,-3.5,-3.5,-3.5,-3.5,
        	4.9,4.9,4.9,4.9,5.2,5.2,5.2,5.2,3.3,3.3,3.3,3.3,-3.2,-3.2,-3.2,-3.2,
        	3.0,3.0,3.0,3.0,3.3,3.3,3.3,3.3,1.4,1.4,1.4,1.4,-5.1,-5.1,-5.1,-5.1,
        	-3.5,-3.5,-3.5,-3.5,-3.2,-3.2,-3.2,-3.2,-5.1,-5.1,-5.1,-5.1,-11.6,-11.6,-11.6,-11.6,
        ],
        # T C
        # A G
        [
        	-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
        	-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
        	-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
        	-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
        	-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
        	-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
        	-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
        	-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
        	-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
        	-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
        	-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
        	-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
        	-3.9,-3.9,-3.9,-3.9,-1.0,-1.0,-1.0,-1.0,-2.8,-2.8,-2.8,-2.8,1.4,1.4,1.4,1.4,
        	-3.6,-3.6,-3.6,-3.6,-0.7,-0.7,-0.7,-0.7,-2.5,-2.5,-2.5,-2.5,1.7,1.7,1.7,1.7,
        	-5.5,-5.5,-5.5,-5.5,-2.6,-2.6,-2.6,-2.6,-4.4,-4.4,-4.4,-4.4,-0.2,-0.2,-0.2,-0.2,
        	-12.0,-12.0,-12.0,-12.0,-9.1,-9.1,-9.1,-9.1,-10.9,-10.9,-10.9,-10.9,-6.7,-6.7,-6.7,-6.7,
        ],
        # T G
        # A C
        [
        	-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
        	-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
        	-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
        	-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
        	-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
        	-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
        	-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
        	-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
        	-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
        	-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
        	-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
        	-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
        	-2.3,-2.3,-2.3,-2.3,-1.9,-1.9,-1.9,-1.9,-2.4,-2.4,-2.4,-2.4,-2.7,-2.7,-2.7,-2.7,
        	-2.0,-2.0,-2.0,-2.0,-1.6,-1.6,-1.6,-1.6,-2.1,-2.1,-2.1,-2.1,-2.4,-2.4,-2.4,-2.4,
        	-3.9,-3.9,-3.9,-3.9,-3.5,-3.5,-3.5,-3.5,-4.0,-4.0,-4.0,-4.0,-4.3,-4.3,-4.3,-4.3,
        	-10.4,-10.4,-10.4,-10.4,-10.0,-10.0,-10.0,-10.0,-10.5,-10.5,-10.5,-10.5,-10.8,-10.8,-10.8,-10.8,
        ],
        # T T
        # A A
        [
        	6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
        	6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
        	4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
        	-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
        	6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
        	6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
        	4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
        	-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
        	6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
        	6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
        	4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
        	-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
        	6.3,6.3,6.3,6.3,7.5,7.5,7.5,7.5,7.0,7.0,7.0,7.0,6.6,6.6,6.6,6.6,
        	6.6,6.6,6.6,6.6,7.8,7.8,7.8,7.8,7.3,7.3,7.3,7.3,6.9,6.9,6.9,6.9,
        	4.7,4.7,4.7,4.7,5.9,5.9,5.9,5.9,5.4,5.4,5.4,5.4,5.0,5.0,5.0,5.0,
        	-1.8,-1.8,-1.8,-1.8,-0.6,-0.6,-0.6,-0.6,-1.1,-1.1,-1.1,-1.1,-1.5,-1.5,-1.5,-1.5,
        ],
])

    # init variable
    dH = dG = 0
    # base = 4000000000
    # T_KELVIN = 273.15
    # DNA_nM = 50

    seq1, seq2 = duplex_str.split("\n")
    region_pos_li = loop_region_dict["region_pos"]
    region_type_li = loop_region_dict["region_type"]
    region_idx = 0
    while region_idx <= len(region_pos_li):
        start = region_pos_li[region_idx - 1] if region_idx != 0 else -1
        end = (
            region_pos_li[region_idx] if region_idx != len(region_pos_li) else len(seq1)
        )
        # stack
        if region_idx % 2 == 0 or region_idx == len(region_pos_li):
            start += 1
            end -= 1
            stack_length = end - start + 1
            # if stack length > 1, participate energy calc
            if stack_length > 1:
                segment = seq1[start : end + 1]
                stack_dh, stack_dg = stack_energy(segment)
                dH += stack_dh
                dG += stack_dg
        # loop
        else:
            # intermolecular initiation energy
            ii_dh, ii_dg = intermolecular_initiation_energy()
            dH += ii_dh
            dG += ii_dg
            loop_idx = region_idx // 2
            loop_type = region_type_li[loop_idx]
            loop_length = end - start + 1
            segment1 = seq1[start - 1:end + 2]
            segment2 = seq2[start - 1:end + 2]
            segment = segment2 if segment1.count("-") else segment1
            seq = seq2[::-1] if segment1.count("-") else seq1
            if not loop_type:  # bulge
                bulge_dh, bulge_dg = bulge_energy(segment, loop_length, seq, start)
                dH += bulge_dh
                dG += bulge_dg
            else:  # int loop
                int_loop_dh, int_loop_dg = int_loop_energy(segment1, segment2, int_loop_energy_dict)
                dH += int_loop_dh
                dG += int_loop_dg 
        region_idx += 1       
        
    # calc Tm
    print(f"\ndG is {dG}")
    dS = (dH - dG) / 310.15 * 1000
    dH *= 1000
    Tm = 0
    print(f"dH is {dH}")
    print(f"dS is {dS}")
    
    # init values
    T_KELVIN = 273.15
    K_mM = 50
    base = 4000000000
    
    # salt params  
    DNA_nM = 50
    dmso_conc = 0
    dmso_fact = 0.6
    formamide_conc = 0.0
    divalent = 1.5
    dntp = 0.6
    
    # symmetry correction if seq is symmetrical
    sym = symmetry(seq1)
    if sym:
        dS += -1.4
        base /= 4
    
    # Terminal AT penalty 
    for i in [seq1[0],seq1[-1]]:
        if i in ["A","T"]:
            dS += 4.1
            dH += 2300
        else:
            dS += -2.8
            dH += 100
    
    
    GC_count = 0 if formamide_conc == 0.0 else (str.count(seq1,"C") + str.count(seq1,"G") + str.count(seq2,"C") + str.count(seq2,"G")) / 2
    K_mM += divalent_to_monovalent(divalent,dntp)
    dS += 0.368 * (len(seq1) - 1) * log(K_mM / 1000.0 )

    Tm = dH / (dS + 1.987 * log(DNA_nM / base)) - T_KELVIN
    Tm -= dmso_conc * dmso_fact
    Tm += (0.453 * GC_count / len(seq1) - 2.88) * formamide_conc
    print(f"Tm is {Tm}")


def main():
    duplex_str = "CAGACG\nGTAGGC"
    duplex_str = "CA--G--CG\nGTGAAAGGC"
    duplex_str = "CA-GACG\nGTGAGGC"
    duplex_str = "GCCCG\nCGG-C"
    # duplex_str = "GAACAG\nCT---C"
    # duplex_str = "TAAACTGCC-----GGCAGTACATC\nATTTGACGGGTGAACCGTCATGTAG"
    # duplex_str = "CCCTATT---------------------------ATTGACGTCAATA\nGGGATAACCGCAATGATACCCTTGTATGCAGTAATAACTGCAGTTAT"
    # duplex_str = "GGCCGGAGTAAGCTGACAT\nCCGGCCTCATTCGACTGTA"
    # duplex_str = "AAAAAAAAAA\nTTTTTTTTTT"
    duplex_str = "CCTGATTCTGTGGATAACC--ATTACCGCCTTTGAGTGAGCTN\nGGACTAAGACACCTATTGGCATAATGGCGGAAACTCACTCGAN"
    loop_region_dict = loop_detective(duplex_str)
    calc_Tm_by_NN(duplex_str, loop_region_dict)
    
    
if __name__ == "__main__":
    main()

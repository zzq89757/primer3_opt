from collections import defaultdict, deque


def loop_detective(duplex_str: str) -> defaultdict:
    """
    检测分子杂交(duplex)中的loop结构(internal loop 和 bulge loop)
    """
    seq1, seq2 = duplex_str.split("\n")
    # ord sum is  AT: 65 + 84 = 149   CG: 67 + 71 = 138
    match_asc_li = [149, 138]
    gap_asc_li = [90, 129, 116, 112, 110]
    ord_sum_li = [ord(x) + ord(y) for x, y in zip(seq1, seq2)]
    # record mismatch and gaps,then recogonize bulge loop and internal loop
    loop_region_dict = defaultdict(
        deque
    )  # region_pos: start and end index, region_type: 0 for bulge,1 for int loop
    region_type_flag = 0  # 0 for bulge,1 for int loop
    flag = 0  # switch to record region start and end index
    print(ord_sum_li)
    for i, v in enumerate(ord_sum_li):
        if v not in match_asc_li and not flag:
            loop_region_dict["region_pos"].append(i)
            flag = 1
        if v not in gap_asc_li and v not in match_asc_li:
            region_type_flag = 1
        if v in match_asc_li and flag:
            loop_region_dict["region_pos"].append(i - 1)
            loop_region_dict["region_type"].append(region_type_flag)
            flag = 0
            region_type_flag = 0
    print(loop_region_dict)
    return loop_region_dict


def basen2int(base: str) -> int:
    """将base(with N)转化为索引号"""
    trantab = str.maketrans("ACGTN", "01234")
    return int(base.upper().translate(trantab), base=5)


def base2int(base: str) -> int:
    """将base(without N)转化为索引号"""
    trantab = str.maketrans("ACGT", "0123")
    return int(base.upper().translate(trantab), base=4)


def stack_energy(segment: str) -> list:
    """近邻法计算stack的总能量值"""
    delta_h = [
        79,
        84,
        78,
        72,
        72,
        85,
        80,
        106,
        78,
        78,
        82,
        98,
        80,
        84,
        80,
        72,
        82,
        85,
        79,
        72,
        72,
        80,
        78,
        72,
        72,
    ]
    delta_s = [
        222,
        224,
        210,
        204,
        224,
        227,
        199,
        272,
        210,
        272,
        222,
        244,
        199,
        224,
        244,
        213,
        222,
        227,
        222,
        227,
        168,
        210,
        220,
        215,
        220,
    ]
    delta_g = [(x - 310.15 * y) / 1000 for x, y in zip(delta_h, delta_s)]
    stack_dh = 0
    stack_ds = 0
    stack_dg = 0
    # print(delta_g)
    for i in range(len(segment) - 1):
        two_mer = segment[i] + segment[i + 1]
        d_index = basen2int(two_mer)
        stack_dh += delta_h[d_index]
        stack_ds += delta_s[d_index]
    return [stack_dh, stack_ds]


def intermolecular_initiation_energy() -> list:
    ii_dh = -7.2  # intermolecular initiation dh
    ii_dg = -1.0  # intermolecular initiation dg
    ii_ds = (ii_dh - ii_dg) / 310.15 * 1000  # intermolecular initiation ds
    return [ii_dh, ii_ds]


def bulge_energy(bulge_length: int) -> list:
    """计算不同长度bulge的总能量值"""
    # bulge loop dh dg
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
    bulge_loop_dict["ds"] = deque(
        [
            ((_dh - _dg) / 310.15) * 1000
            for _dh, _dg in zip(bulge_loop_dict["dh"], bulge_loop_dict["dg"])
        ]
    )
    return [
        bulge_loop_dict["dh"][bulge_length - 1],
        bulge_loop_dict["ds"][bulge_length - 1],
    ]


def symmetric_int_loop_energy(loop_type: list) -> list:
    """对称的int loop结构直接读取矩阵数据累加即可"""
    symmetric_int_loop_dh = symmetric_int_loop_ds = 0
    
    return [symmetric_int_loop_dh, symmetric_int_loop_ds]


def asymmetry_correct_energy(loop_diff_abs: int) -> list:
    asymmetry_dg = 0.4
    asymmetry_dh = 0
    asymmetry_ds = ((asymmetry_dh - asymmetry_dg) / 310.15) * 1000
    asymmetry_correct_dh = asymmetry_dh * loop_diff_abs
    asymmetry_correct_ds = asymmetry_ds * asymmetry_ds
    return [asymmetry_correct_dh, asymmetry_correct_ds]


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
    initiation_ds = [
        ((x - y) / 310.15) * 1000 for x, y in zip(initiation_dh, initiation_dg)
    ]
    return [initiation_dh[loop_sum - 4], initiation_ds[loop_sum - 4]]


def asymmetric_int_loop_mismatch_energy(
    segment1: str, segmnet2: str, int_loop_type: list
) -> list:
    """只有非对称且min loop > 1 才考虑mismatch"""
    # min loop length == 1,do NOT consider mismatch
    if int_loop_type[0] == 1:
        return [0.0, 0.0]
    # mismatch energy calc by segment
    


def asymmetric_int_loop_energy(
    segment1: str, segment2: str, loop_sum: int, loop_diff_abs: int, loop_type: list
) -> list:
    asymmetric_int_loop_dh = asymmetric_int_loop_ds = 0
    # intermolecular initiation energy
    ii_dh, ii_ds = intermolecular_initiation_energy()
    # loop sum initiation energy
    li_dh, li_ds = asymmetric_int_loop_initiation_energy(loop_sum)
    # asymmetry correct energy
    asymmetry_correct_dh, asymmetry_correct_ds = asymmetry_correct_energy(loop_diff_abs)
    # mismatch energy
    mm_dh, mm_ds = asymmetric_int_loop_mismatch_energy(segment1, segment2, loop_type)


def int_loop_energy(segment1: str, segment2: str) -> list:
    """对intloop进行分型并计算能量"""
    int_loop_dh = 0
    int_loop_ds = 0
    int_loop_dg = 0

    first_loop_length = len(segment1) - segment1.count("-")
    second_loop_length = len(segment2) - segment2.count("-")
    loop_sum = first_loop_length + second_loop_length
    max_loop_length = max(first_loop_length, second_loop_length)
    min_loop_length = loop_sum - max_loop_length
    loop_type = [min_loop_length, max_loop_length]
    is_symmetric = loop_sum <= 4 and max_loop_length < 3

    # 1×1, 1×2, 2×2 Internal Loops
    if is_symmetric:
        symmetric_int_loop_dh, symmetric_int_loop_ds = symmetric_int_loop_energy(loop_type)
        int_loop_dh += symmetric_int_loop_dh
        int_loop_ds += symmetric_int_loop_ds
    # Other Internal Loops
    else:
        loop_diff_abs = abs(first_loop_length - second_loop_length)
        asymmetric_int_loop_dh, asymmetric_int_loop_ds = asymmetric_int_loop_energy(
            segment1, segment2, loop_sum, loop_diff_abs, loop_type
        )
        int_loop_dh += asymmetric_int_loop_dh
        int_loop_ds += asymmetric_int_loop_ds

    return [int_loop_dh, int_loop_ds]


def calc_Tm_by_NN(duplex_str: str, loop_region_dict: defaultdict) -> float:
    # int loop energy matrix
    int_loop_energy_dict = defaultdict(deque)
    int_loop_energy_dict['11dg'] = deque([
        # A
        # T 
        [
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
                0.8, 1.4, -0.4, 1.0, # AA AC AG AT
                1.6, 2.1, 1.0, 1.6, # CA CC CG CT
                -0.2, 1.0, -1.2, 2.0, # GA GC GG GT
                1.0, 1.4, 2.0, 1.1, # TA TC TG TT    
            ],
            # A G
            # T C
            [
                1.0, 1.7, 0.3, 1.0, # AA AC AG AT
                1.5, 2.0, 1.0, 1.0, # CA CC CG CT
                0.1, 1.0, -0.3, 2.0, # GA GC GG GT
                1.0, 1.4, 2.0, 0.6, # TA TC TG TT 
            ],
            # A T
            # T A
            [
                1.2, 1.7, 0.2, 1.0, # AA AC AG AT
                1.7, 2.7, 1.0, 1.4, # CA CC CG CT
                0.2, 1.0, -0.3, 2.0, # GA GC GG GT
                1.0, 1.4, 2.0, 1.4, # TA TC TG TT 
            ]
            
        ],
        # C
        # G
        [
            # C A
            # G T
            [
                
            ]
        ]
        
    ])
    
    # init variable
    dH = dS = 0
    
    

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
                stack_dh, stack_ds = stack_energy(segment)
                dH += stack_dh
                dS += stack_ds
            print(f"stack {start}->{end}")
        # loop
        else:
            loop_idx = region_idx // 2
            loop_type = region_type_li[loop_idx]
            loop_length = end - start + 1
            if not loop_type:  # bulge
                bulge_dh, bulge_ds = bulge_energy(loop_length)
                dH += bulge_dh
                dS += bulge_ds
                print(f"bulge {start} -> {end}")
            else:  # int loop
                int_loop_dh, int_loop_ds = int_loop_energy()
                dH += int_loop_dh
                dS += int_loop_ds
                print(f"intloop {start} -> {end}")

        region_idx += 1


if __name__ == "__main__":
    loop_region_dict = loop_detective("GCTAGCATCGTA--GCTCGA\nCGTAGCTGATGCTTGTAGCT")
    calc_Tm_by_NN("GCTAGCATCGTA--GCTCGA\nCGTAGCTGATGCTTGTAGCT", loop_region_dict)

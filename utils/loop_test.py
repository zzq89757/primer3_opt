from collections import defaultdict, deque


def loop_detective(duplex_str: str) -> defaultdict:
    """
    检测分子杂交(duplex)中的loop结构(internal loop 和 bulge loop)

    """
    seq1, seq2 = duplex_str.split("\n")
    match_asc_li = [149, 138]  # ord sum is  AT: 65 + 84 = 149   CG: 67 + 71 = 138
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


def int_loop_type():
    """对intloop进行分型"""
    ...


def delta_bulge():
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


def delta_stack(segment: str, stack_dict: defaultdict, thermodynamics_params: str): ...


def calc_Tm_by_NN(duplex_str: str, loop_region_dict: defaultdict) -> float:
    dS = 0
    dG = 0
    dH = 0
    ii_dh = -7.2  # intermolecular initiation dh
    ii_dg = -1.0  # intermolecular initiation dg

    seq1, seq2 = duplex_str.split("\n")
    region_pos_li = loop_region_dict["region_pos"]
    region_type_li = loop_region_dict["region_type"]
    region_idx = 0
    while region_idx <= len(region_pos_li):
        start = region_pos_li[region_idx - 1] if region_idx != 0 else -1
        end = (
            region_pos_li[region_idx] if region_idx != len(region_pos_li) else len(seq1)
        )
        # stack region
        if region_idx % 2 == 0 or region_idx == len(region_pos_li):
            start += 1
            end -= 1
            # for in start to end, do something
            dH += delta_stack(segment, stack_dict, "dh")
            dG += delta_stack(segment, stack_dict, "dg")
            print(f"stack {start}->{end}")
        # loop region
        else:
            loop_idx = region_idx // 2
            loop_type = region_type_li[loop_idx]
            # bulge
            if not loop_type:
                print(f"bulge {start} -> {end}")
            else:  # int loop
                print(f"intloop {start} -> {end}")

        region_idx += 1


def sum_test() -> None:
    asc_li = [65, 67, 71, 84, 45]
    base_li = ["A", "C", "G", "T", "-"]
    for i in range(len(asc_li)):
        for j in range(i, len(asc_li)):
            print(f"{base_li[i]}>{base_li[j]}", end="\t")
            print(asc_li[i] + asc_li[j])


# sum_test()
# exit()


def md_tag_detective(md_tag: str) -> list[str]:
    """
    md标签反应参考上的碱基情况
    当出现^时 read相较于参考有缺失 直接增加^后的碱基数 保证位置不偏移
    若read相较于参考有插入 需看CigarI标签 不添加碱基数
    """
    import re

    # md_tag_li = re.split(r'([ACTG^])', md_tag)
    # md_tag_li = [x for x in  re.split(r'(\d+)', md_tag) if x]
    md_tag_li = re.split(r"(\d+)", md_tag)
    md_tag_li = list(filter(None, md_tag_li))
    print(md_tag_li)


def cigar_detective(reference_start: int, cigar_str: str) -> defaultdict:
    import re

    insertion_dict = defaultdict()
    # extract to a indival function:has_insertion
    if cigar_str.find("I") == -1:
        return insertion_dict
    # only preserve cigar_str until last I
    suffix = cigar_str.rpartition("I")[2]
    cigar_str = cigar_str.removesuffix(suffix)
    cigar_li = re.split(r"(\d+)", cigar_str)
    cigar_li = list(filter(None, cigar_li))
    current_pos = reference_start
    for idx in range(0, len(cigar_li), 2):
        region_length = cigar_li[idx]
        tag = cigar_li[idx + 1]
        # I not add pos, D add
        if tag != "I":
            current_pos += int(region_length)
        else:
            insertion_dict[current_pos] = region_length
    print(insertion_dict)


def snp_check(reference_start: int, cigar_str: str, md_tag: str | None) -> None: ...


if __name__ == "__main__":
    # re_test("506^CA9C17A0C0A18^G9^C11")
    # cigar_test(100, "89M1I11M38S22M5I30M")
    loop_region_dict = loop_detective("GCTAGCATCGTA--GCTCGA\nCGTAGCTGATGCTTGTAGCT")
    calc_Tm_by_NN("GCTAGCATCGTA--GCTCGA\nCGTAGCTGATGCTTGTAGCT", loop_region_dict)

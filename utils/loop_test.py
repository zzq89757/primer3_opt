from collections import defaultdict, deque


def loop_detective(duplex_str: str) -> defaultdict:
    """
    检测分子杂交(duplex)中的loop结构(internal loop 和 bulge loop)

    """
    match_asc_li = [149, 138]  # ord sum is  AT: 65 + 84 = 149   CG: 67 + 71 = 138
    gap_asc_li = [90, 129, 116, 112, 110]
    seq1, seq2 = duplex_str.split("\n")
    ord_sum_li = [ord(x) + ord(y) for x, y in zip(seq1, seq2)]
    # record mismatch and gaps,then recogonize bulge loop and internal loop
    loop_region_dict = defaultdict(
        list
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


def calc_Tm_by_NN(duplex_str: str, loop_region_dict: defaultdict) -> float:
    seq1, seq2 = duplex_str.split("\n")
    for i in range(len(seq1)):
        ...


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
    loop_detective("GCTAGCATCGTA--GCT\nCGTAGCTGATGCTTGTA")

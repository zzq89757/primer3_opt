from collections import defaultdict, deque


def loop_detective(duplex_str: str) -> None:
    """
    检测分子杂交(duplex)中的loop结构(internal loop和bulge loop)

    """
    gap_index_li = deque()
    seq1, seq2 = duplex_str.split("\n")
    # record pairs and gaps,then recogonize bulge loop and internal loop
    for idx in range(len(seq1)):
        # gap found, record index
        if seq1[idx] == '-' or seq2[idx] == '-':
            gap_index_li.append(idx)
        
    # cycle 2, calc 
            


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
    loop_detective("GCTAGCATCGTAGCT\nCGTAGCTGATGCGTA")

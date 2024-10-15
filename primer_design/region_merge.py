def get_final_included_region(includes_left, includes_right, excludes_left, excludes_right):
    # 将包含和排除区域转换为 (start, end) 形式，方便处理
    def to_intervals(regions):
        return [(start, start + length) for start, length in regions]

    # 合并区间
    def merge_intervals(intervals):
        if not intervals:
            return []
        intervals.sort(key=lambda x: x[0])
        merged = [intervals[0]]
        for current in intervals[1:]:
            prev_start, prev_end = merged[-1]
            curr_start, curr_end = current
            if curr_start <= prev_end:
                merged[-1] = (prev_start, max(prev_end, curr_end))  # 合并重叠区间
            else:
                merged.append(current)
        return merged

    # 执行区间的减法：从包含的区间中减去排除的区间
    def subtract_intervals(included, excluded):
        result = []
        for inc_start, inc_end in included:
            temp_start = inc_start
            for exc_start, exc_end in excluded:
                # 如果没有重叠，跳过
                if exc_end <= temp_start or exc_start >= inc_end:
                    continue
                # 有重叠，分割区间
                if exc_start > temp_start:
                    result.append((temp_start, exc_start))  # 添加排除前的部分
                temp_start = max(temp_start, exc_end)  # 更新当前起始
            if temp_start < inc_end:
                result.append((temp_start, inc_end))  # 添加剩余的部分
        return result

    # 转换区域
    include_intervals_left = to_intervals(includes_left)
    include_intervals_right = to_intervals(includes_right)
    exclude_intervals_left = to_intervals(excludes_left)
    exclude_intervals_right = to_intervals(excludes_right)

    # 合并重叠的包含和排除区域
    merged_includes_left = merge_intervals(include_intervals_left)
    merged_includes_right = merge_intervals(include_intervals_right)
    merged_excludes_left = merge_intervals(exclude_intervals_left)
    merged_excludes_right = merge_intervals(exclude_intervals_right)

    # 减去排除区域，得到最终的包含区间
    final_included_left = subtract_intervals(merged_includes_left, merged_excludes_left)
    final_included_right = subtract_intervals(merged_includes_right, merged_excludes_right)

    return final_included_left, final_included_right

# 示例使用
includes_left = [[30, 10]]  # 左侧包含区域，(起始, 长度)
includes_right = [[50, 5]]  # 右侧包含区域，(起始, 长度)
excludes_left = [[35, 2]]    # 左侧排除区域，(起始, 长度)
excludes_right = [[52, 3]]   # 右侧排除区域，(起始, 长度)

final_forward_region, final_reverse_region = get_final_included_region(includes_left, includes_right, excludes_left, excludes_right)

print("最终包含区域（左）：", final_forward_region)
print("最终包含区域（右）：", final_reverse_region)

# fill if len uneq
# if len(final_forward_region) > len(final_reverse_region):
#     final_reverse_region += abs(len(final_forward_region) - len(final_reverse_region)) * tuple([final_reverse_region[0]])

def fill_region(final_forward_region: list, final_reverse_region: list):
    if len(final_forward_region) == len(final_reverse_region):
        return final_forward_region, final_reverse_region
    region_to_fill = final_forward_region if len(final_forward_region) < len(final_reverse_region) else final_reverse_region
    # final_reverse_region += abs(len(final_forward_region) - len(final_reverse_region)) * tuple([final_reverse_region[0]])
    region_to_fill += abs(len(final_forward_region) - len(final_reverse_region)) * tuple([region_to_fill[0]])
    return final_forward_region, final_reverse_region

final_forward_region,final_reverse_region = fill_region(final_forward_region,final_reverse_region)
    

print("最终包含区域（左）：", final_forward_region)
print("最终包含区域（右）：", final_reverse_region)

def generate_pair_ok_regions(final_forward_region: list, final_reverse_region: list) -> list:
    pairs_ok_region_li = []
    for forward_region, reverse_region in zip(final_forward_region, final_reverse_region):
        pairs_ok_region_li.append(forward_region + reverse_region)
    return pairs_ok_region_li

generate_pair_ok_regions(final_forward_region,final_reverse_region)
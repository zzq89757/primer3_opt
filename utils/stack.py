import numpy as np
from collections import deque
stack_dh = deque([
        -7.9, -8.4, -7.8, -7.2, # AA AC AG AT AN
        -8.5, -8.0, -10.6, -7.8, # CA CC CG CT CN
        -8.2, -9.8, -8.0, -8.4, # GA GC GG GT GN
        -7.2, -8.2, -8.5, -7.9, # TA TC TG TT TN
    ])
stack_ds = deque([
        -22.2, -22.4, -21.0, -20.4, # AA AC AG AT AN
        -22.7, -19.9, -27.2, -21.0, # CA CC CG CT CN
        -22.2, -24.4, -19.9, -22.4, # GA GC GG GT GN
        -21.3, -22.2, -22.7, -22.2, # TA TC TG TT TN
    ])

def append_stack(stack_table: list):
    arr = np.array(stack_table)

    # 将一维数组转换为4x4的二维数组
    matrix = arr.reshape(4, 4)

    # 计算每行的平均值
    row_avg = np.mean(matrix, axis=1)
    
    # 插入原列表
    insert_dix = 4
    for i in range(4):
        stack_table.insert(insert_dix,round(row_avg[i],1))
        insert_dix += 5

    # 计算每列的平均值
    col_avg = np.mean(matrix, axis=0)
    # 插入原列表
    for i in col_avg:
        stack_table.append(round(i,1))

    # 计算整个矩阵的平均值
    total_avg = np.mean(matrix)
    stack_table.append(round(total_avg,1))
    print(stack_table)



def main() -> None:
    append_stack(stack_ds)
    
if __name__ == "__main__":
    main()
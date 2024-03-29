import primer3
import pandas as pd
import yaml

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)

yaml_dict = yaml.safe_load(open("/home/waynezheng/Code/Bundle/Lab/Primer3/config/config.yaml"))

# Primer3 序列和设计参数，必须要有
seq_args = {
        'SEQUENCE_ID': 'BRCA2_SEGMENT',
        'SEQUENCE_TEMPLATE': "GTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGGTCGCCGCCGGGAGAAGCGTGAGGGGACAGATTTGTGACCGGCGCGGTTTTTGTCAGCTTACTCCGGCCAAAAAAGAACTGCACCTCTGGAGCGGGTTAGTGGTGGTGGTAGTGGGTGTGGCGCGAGCTTCTGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCTCTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGGTCGCCGCCGGGAGAAGCGTGAGGGGACAGATTTGTGACCGGCGCGGTTTTTGTCAGCTTACTCCGGCCAAAAAAGAACTGCACCTCTGGAGCGGGTTAGTGGTGGTGGTAGTGGGT",
        # SEQ range(第二个参数是长度)
        # 'SEQUENCE_INCLUDED_REGION': [20,180],
        # EXCLUDE range(第二个参数是长度)
        'SEQUENCE_EXCLUDED_REGIONSEQUENCE_EXCLUDED_REGION':[],
        # 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':[[50,50,150,50], [200,60,300,50]],
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':[[50,50,150,50]]
        # 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':[200,60,-1,-1]
        # 'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST':100,50,300,50 ; 900,60,, ; ,,930,100
        # 表明有引物设计有三种选择：
        # 左引物在100～150bp区间进行设计，右引物在300～350bp的区间进行设计；
        # 左引物在900～960bp区间进行设计，右引物随意；
        # 右引物在930～1030bp区间进行设计，左引物随意。
}
seq_args=yaml_dict['SEQ_ARG']
# primer3.designPrimers()
# Primer3 全局参数，这个可选
global_args = {
        # NUM return
        # 'PRIMER_NUM_RETURN':5,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        # SIZE range
        # 'PRIMER_OPT_SIZE': 20,
        # 'PRIMER_MIN_SIZE': 18,
        # 'PRIMER_MAX_SIZE': 25,
        # # TM range
        # 'PRIMER_OPT_TM': 60.0,
        # 'PRIMER_MIN_TM': 57.0,
        # 'PRIMER_MAX_TM': 63.0,
        # # GC range
        # 'PRIMER_MIN_GC': 20.0,
        # 'PRIMER_MAX_GC': 80.0,
        
        # 'PRIMER_MAX_POLY_X': 3,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        # FR primer TM diff
        'PRIMER_PAIR_MAX_DIFF_TM':0.5,
        # 数组前后顺序决定优先级 数组放相同数字即为最优产物长度设定
        # 'PRIMER_PRODUCT_SIZE_RANGE': [[150,150],[100,125],[125,150],[175,200]],
}
global_args.update(yaml_dict['GLOBAL_ARG'])
# 执行命令并返回结果
primer3_result = primer3.designPrimers(seq_args, global_args)
# print(primer3_result)

# parse result

# over view

num_return=primer3_result['PRIMER_PAIR_NUM_RETURNED']

# primer details
for i in range(num_return):
  # primer seq 
  l_seq=primer3_result[f"PRIMER_LEFT_{i}_SEQUENCE"]
  r_seq=primer3_result[f"PRIMER_RIGHT_{i}_SEQUENCE"]
  # primer tm
  l_tm=round(primer3_result[f"PRIMER_LEFT_{i}_TM"],2)
  r_tm=round(primer3_result[f"PRIMER_RIGHT_{i}_TM"],2)
  # gc percent
  l_gc=primer3_result[f"PRIMER_LEFT_{i}_GC_PERCENT"]
  r_gc=primer3_result[f"PRIMER_RIGHT_{i}_GC_PERCENT"]
  # prod
  prod_len=primer3_result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"]
  # filter
#   print(l_seq,end="\t")
#   print(l_tm,end="\t")
#   print(r_seq,end="\t")
#   print(prod_len,end="\t")
#   print(r_tm)
i=1
lsame_list=[f"PRIMER_LEFT_{i}_"]*3
rsame_list=[f"PRIMER_RIGHT_{i}_"]*3
para_li=["SEQUENCE","TM","GC_PERCENT"]
linfo_list=[lsame_list[i] + para_li[i] for i in range(len(lsame_list))]
rinfo_list=[rsame_list[i] + para_li[i] for i in range(len(rsame_list))]
info_list=linfo_list+rinfo_list
info_list.append(f"PRIMER_PAIR_{i}_PRODUCT_SIZE")
print(info_list)
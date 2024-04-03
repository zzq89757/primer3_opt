from dataclasses import dataclass
import primer3

@dataclass

class Primer3:
  seq_args:dict
  alternative_args:dict
  
  def execute(self):
    self.deal_paramas()
    self.deal_result()
  
  def deal_paramas(self):
    self.global_args={
        # NUM return
        # 'PRIMER_NUM_RETURN':5,
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
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,       
        'PRIMER_MAX_POLY_X': 4,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        # FR primer TM diff
        'PRIMER_PAIR_MAX_DIFF_TM':0.5
}
    self.global_args.update(self.alternative_args)
    
  def deal_result(self):
    self.primer3_result = primer3.design_primers(self.seq_args, self.global_args)
    # print(self.primer3_result)
    self.primer_pair_num=self.primer3_result['PRIMER_PAIR_NUM_RETURNED']
    self.result_dict = []
    # 设置需要的结果列表 并从primer3结果中提取
    for i in range(self.primer_pair_num):
      para_li=["","_SEQUENCE","_TM","_GC_PERCENT"]
      primer_res_len=len(para_li)
      lsame_list=[f"PRIMER_LEFT_{i}"] * primer_res_len
      rsame_list=[f"PRIMER_RIGHT_{i}"] * primer_res_len
      linfo_list=[lsame_list[i] + para_li[i] for i in range(len(lsame_list))]
      rinfo_list=[rsame_list[i] + para_li[i] for i in range(len(rsame_list))]
      info_list=linfo_list+rinfo_list
      info_list.append(f"PRIMER_PAIR_{i}_PRODUCT_SIZE")
      self.result_dict.append(dict([(key,self.primer3_result[key]) for key in info_list]))
    return self.result_dict
  
if __name__=="__main__":
  import yaml
  yaml_dict = yaml.safe_load(open("/home/wayne/Code/Bundle/Lab/Primer3/config/config.yaml"))
  input_seq_arg=yaml_dict['SEQ_ARG']
  user_arg=yaml_dict['GLOBAL_ARG']
  myres=Primer3(input_seq_arg,user_arg)
  myres.deal_paramas()
  # myres.deal_result()
  print(myres.deal_result())
  # print(myres.primer3_result)
from math import sqrt,log

# Tables of nearest-neighbor thermodynamics for DNA bases, from the
# paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
# and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
# Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]

delta_h={
"DH_A_A" : 79,
"DH_A_C" : 84,
"DH_A_G" : 78,
"DH_A_T" : 72,
"DH_A_N" : 72,

"DH_C_A" : 85,
"DH_C_C" : 80,
"DH_C_G" : 106,
"DH_C_T" : 78,
"DH_C_N" : 78,

"DH_G_A" : 82,
"DH_G_C" : 98,
"DH_G_G" : 80,
"DH_G_T" : 84,
"DH_G_N" : 80,

"DH_T_A" : 72,
"DH_T_C" : 82,
"DH_T_G" : 85,
"DH_T_T" : 79,
"DH_T_N" : 72,

"DH_N_A" : 72,
"DH_N_C" : 80,
"DH_N_G" : 78,
"DH_N_T" : 72,
"DH_N_N" : 72
}


delta_s={
"DS_A_A" : 222,
"DS_A_C" : 224,
"DS_A_G" : 210,
"DS_A_T" : 204,
"DS_A_N" : 224,

"DS_C_A" : 227,
"DS_C_C" : 199,
"DS_C_G" : 272,
"DS_C_T" : 210,
"DS_C_N" : 272,

"DS_G_A" : 222,
"DS_G_C" : 244,
"DS_G_G" : 199,
"DS_G_T" : 224,
"DS_G_N" : 244,

"DS_T_A" : 213,
"DS_T_C" : 222,
"DS_T_G" : 227,
"DS_T_T" : 222,
"DS_T_N" : 227,

"DS_N_A" : 168,
"DS_N_C" : 210,
"DS_N_G" : 220,
"DS_N_T" : 215,
"DS_N_N" : 220
}

# End of tables nearest-neighbor parameter ------------------------------

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



def calc_tm(seq:str):
  '''
  Calculate tm by SantaLucia method and correction
  '''
  # init values
  T_KELVIN = 273.15
  K_mM = 50
  ds = 0
  dh = 0
  Tm = 0
  base = 4000000000
  
  # primer3 default params  
  DNA_nM = 50
  dmso_conc = 0
  dmso_fact = 0.6
  formamide_conc = 0.8
  divalent = 1.5
  dntp = 0.6
  # monovalent = 50 

  # symmetry correction if seq is symmetrical
  sym = symmetry(seq)
  if sym:
    ds += 14
    base /= 4
  
  # Terminal AT penalty 
  for i in [seq[0],seq[-1]]:
    if i in ["A","T"]:
      ds += -41
      dh += -23
    else:
      ds += 28
      dh += -1
  

  # calculate delta by NN 
  for i,s in enumerate(seq):
    if  i == 0 :continue
    last_base = seq[i-1]
    dh_index_name = f"DH_{last_base}_{s}"
    dh += delta_h[dh_index_name]
    ds_index_name = f"DS_{last_base}_{s}"
    ds += delta_s[ds_index_name]
    
  
  # init value and salt corrections and calculate Tm finally
  dh *= -100
  ds *= -0.1

  GC_count = 0 if formamide_conc == 0.0 else str.count(seq,"C") + str.count(seq,"G")
  K_mM += divalent_to_monovalent(divalent,dntp)
  ds=ds + 0.368 * (len(seq) - 1) * log(K_mM / 1000.0 )

  Tm = dh / (ds + 1.987 * log(DNA_nM / base)) - T_KELVIN
  Tm -= dmso_conc * dmso_fact
  Tm += (0.453 * GC_count / len(seq) - 2.88) * formamide_conc
  return Tm


if __name__=="__main__":
  import primer3
  tm=calc_tm("AAAAAAAAAA")
  print(tm)
  print(primer3.calc_tm("AAAAAAAAAA"))
  tm=calc_tm("GGCCGGAGTAAGCTGACA")
  print(tm)
  print(primer3.calc_tm("GGCCGGAGTAAGCTGACA"))
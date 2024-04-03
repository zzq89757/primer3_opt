from math import log10
import primer3
print(primer3.bindings.calc_tm("GGCCGGAGTAAGCTGACA",))
print(primer3.bindings.calc_tm("CGGGCCGGAGTAAGCTGACAGGCCGGAGTAAGCTGACAGGCCGGAGTAAGCTGACAGGCCGGAGTAAGCTGACA"))

def cal_tm(seq,mv_conc):
  gc_percent=(str.count(seq,"C")+ str.count(seq,"G"))/len(seq) * 100
  # Tm = 81.5 + 16.6(log10([mv_conc])) + 0.41(%GC) - 600/length mv_conc default:50
  Tm=81.5+16.6*log10(50)+0.41*gc_percent-600/len(seq)
  # Tm=81.5+16.6*log10(50)+0.41*gc_percent-
  print(Tm)
  return 0

if __name__=="__main__":
  my_seq="GGCCGGAGTAAGCTGACA"
  my_seq="CGGGCCGGAGTAAGCTGACAGGCCGGAGTAAGCTGACAGGCCGGAGTAAGCTGACAGGCCGGAGTAAGCTGACA"
  cal_tm(my_seq,None)
from math import sqrt,log

# Tables of nearest-neighbor thermodynamics for DNA bases, from the
# paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
# and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
# Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]

delta_h = [79, 84, 78, 72, 72, 85, 80, 106, 78, 78, 82, 98, 80, 84, 80, 72, 82, 85, 79, 72, 72, 80, 78, 72, 72]

delta_s = [222, 224, 210, 204, 224, 227, 199, 272, 210, 272, 222, 244, 199, 224, 244, 213, 222, 227, 222, 227, 168, 210, 220, 215, 220]

# End of tables nearest-neighbor parameter ------------------------------

def base2int(base:str) -> int:
    trantab = str.maketrans('ACGTN', '01234')
    return int(base.upper().translate(trantab), base=5)

def base10to6(num):
    l = []
    while True:
        num, remainder = divmod(num, 5)
        l.append(str(remainder))
        if num == 0:return ''.join(l[::-1]) 
          
def int2base(int:int) -> str:
    trantab = str.maketrans('01234', 'ACGTN')
    int6 = base10to6(int).zfill(2)
    return int6.translate(trantab)


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
        two_mer = seq[i-1] + s
        d_index = base2int(two_mer)
        dh += delta_h[d_index]
        ds += delta_s[d_index]
        
    
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
    tm=calc_tm("GGCCGGAGTAAGCTGACAT")
    print(tm)
    print(primer3.calc_tm("GGCCGGAGTAAGCTGACAT"))
from math import sqrt,log


class Primer():
    
    def __init__(self, seq) -> None:
        # Tables of nearest-neighbor thermodynamics for DNA bases, from the
        # paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
        # and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
        # Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]

        self.delta_h = [79, 84, 78, 72, 72, 85, 80, 106, 78, 78, 82, 98, 80, 84, 80, 72, 82, 85, 79, 72, 72, 80, 78, 72, 72]

        self.delta_s = [222, 224, 210, 204, 224, 227, 199, 272, 210, 272, 222, 244, 199, 224, 244, 213, 222, 227, 222, 227, 168, 210, 220, 215, 220]

        # End of tables nearest-neighbor parameter ------------------------------
        
        # init values
        self.seq = seq
        self.T_KELVIN = 273.15
        self.K_mM = 50
        self.ds = 0
        self.dh = 0
        self.Tm = 0
        self.base = 4000000000
        self.GC_count = 0
        
        # primer3 default params  
        self.DNA_nM = 50
        self.dmso_conc = 0
        self.dmso_fact = 0.6
        self.formamide_conc = 0.8
        self.divalent = 1.5
        self.dntp = 0.6
        # self.monovalent = 50
        
        # cal Tm
        self.calc_tm()
        

    def base2int(self, base:str) -> int:
        trantab = str.maketrans('ACGTN', '01234')
        return int(base.upper().translate(trantab), base=5)

    def base10to6(self, num):
        l = []
        while True:
            num, remainder = divmod(num, 5)
            l.append(str(remainder))
            if num == 0:return ''.join(l[::-1]) 
            
    def int2base(self, int:int) -> str:
        trantab = str.maketrans('01234', 'ACGTN')
        int6 = self.base10to6(int).zfill(2)
        return int6.translate(trantab)


    def symmetry(self):
        '''
        Return 1 if string is symmetrical, 0 otherwise.
        '''
        seq_len = len(self.seq)
        mp = seq_len // 2
        if seq_len % 2 == 1:return 0
        for i in range(mp):
            s = self.seq[i]
            e = self.seq[seq_len - i - 1]
            terminal_base_set={s,e}
            if terminal_base_set != {"A","T"} and terminal_base_set != {"C","G"}:return 0
        return 1


    def divalent_to_monovalent(self):
        '''
        Convert divalent salt concentration to monovalent
        '''
        if self.divalent == 0:self.dntp = 0
        if self.divalent < self.dntp:self.divalent = self.dntp
        return 120 * sqrt((self.divalent - self.dntp))



    def calc_tm(self):
        '''
        Calculate tm by SantaLucia method and correction
        '''

        # symmetry correction if seq is symmetrical
        sym = self.symmetry()
        if sym:
            self.ds += 14
            self.base /= 4
        
        # Terminal AT penalty 
        for i in [self.seq[0],self.seq[-1]]:
            if i in ["A","T"]:
                self.ds += -41
                self.dh += -23
            else:
                self.ds += 28
                self.dh += -1
        

        # calculate delta by NN 
        for i,s in enumerate(self.seq):
            if  i == 0 :continue
            two_mer = self.seq[i-1] + s
            d_index = self.base2int(two_mer)
            self.dh += self.delta_h[d_index]
            self.ds += self.delta_s[d_index]
            
        
        # init value and salt corrections and calculate Tm finally
        self.dh *= -100
        self.ds *= -0.1

        self.GC_count = 0 if self.formamide_conc == 0.0 else str.count(self.seq,"C") + str.count(self.seq,"G")
        self.K_mM += self.divalent_to_monovalent()
        self.ds = self.ds + 0.368 * (len(self.seq) - 1) * log(self.K_mM / 1000.0 )

        self.Tm = self.dh / (self.ds + 1.987 * log(self.DNA_nM / self.base)) - self.T_KELVIN
        self.Tm -= self.dmso_conc * self.dmso_fact
        self.Tm += (0.453 * self.GC_count / len(self.seq) - 2.88) * self.formamide_conc
        return self.Tm


if __name__=="__main__":
    import primer3
    primer1=Primer("AAAAAAAAAA")
    print(primer1.Tm)
    print(primer3.calc_tm("AAAAAAAAAA"))
    primer2 = Primer("GGCCGGAGTAAGCTGACAT")
    print(primer2.Tm)
    print(primer3.calc_tm("GGCCGGAGTAAGCTGACAT"))
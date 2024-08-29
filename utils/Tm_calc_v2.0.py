from math import sqrt,log


class Primer():
    
    def __init__(self, seq, template) -> None:
        # Tables of nearest-neighbor thermodynamics for DNA bases, from the
        # paper [SantaLucia JR (1998) "A unified view of polymer, dumbbell
        # and oligonucleotide DNA nearest-neighbor thermodynamics", Proc Natl
        # Acad Sci 95:1460-65 http://dx.doi.org/10.1073/pnas.95.4.1460]

        # match table of dH dS
        self.dH_m = [79, 84, 78, 72, 72, 85, 80, 106, 78, 78, 82, 98, 80, 84, 80, 72, 82, 85, 79, 72, 72, 80, 78, 72, 72]

        self.dS_m = [222, 224, 210, 204, 224, 227, 199, 272, 210, 272, 222, 244, 199, 224, 244, 213, 222, 227, 222, 227, 168, 210, 220, 215, 220]

        # hybridition match table of dH dS
        self.dH_mm = [
        99999,99999,99999,4700,76174, # AA/AA AA/AC AA/AG AA/AT AA/AN 
        99999,99999,99999,7600,76899, # AA/CA AA/CC AA/CG AA/CT AA/CN 
        99999,99999,99999,3000,75749, # AA/GA AA/GC AA/GG AA/GT AA/GN 
        1200,2300,-600,-7900,-1250, # AA/TA AA/TC AA/TG AA/TT AA/TN 
        75299,75574,74849,1850,56893, # AA/NA AA/NC AA/NG AA/NT AA/NN 
        99999,99999,-2900,99999,74274, # AC/AA AC/AC AC/AG AC/AT AC/AN 
        99999,99999,-700,99999,74824, # AC/CA AC/CC AC/CG AC/CT AC/CN 
        99999,99999,500,99999,75124, # AC/GA AC/GC AC/GG AC/GT AC/GN 
        5300,0,-8400,700,-600, # AC/TA AC/TC AC/TG AC/TT AC/TN 
        76324,74999,-2875,75174,55905, # AC/NA AC/NC AC/NG AC/NT AC/NN 
        99999,-900,99999,99999,74774, # AG/AA AG/AC AG/AG AG/AT AG/AN 
        99999,600,99999,99999,75149, # AG/CA AG/CC AG/CG AG/CT AG/CN 
        99999,-4000,99999,99999,73999, # AG/GA AG/GC AG/GG AG/GT AG/GN 
        -700,-7800,-3100,1000,-2650, # AG/TA AG/TC AG/TG AG/TT AG/TN 
        74824,-3025,74224,75249,55318, # AG/NA AG/NC AG/NG AG/NT AG/NN 
        1200,99999,99999,99999,75299, # AT/AA AT/AC AT/AG AT/AT AT/AN 
        5300,99999,99999,99999,76324, # AT/CA AT/CC AT/CG AT/CT AT/CN 
        -700,99999,99999,99999,74824, # AT/GA AT/GC AT/GG AT/GT AT/GN 
        -7200,-1200,-2500,-2700,-3400, # AT/TA AT/TC AT/TG AT/TT AT/TN 
        -350,74699,74374,74324,55761, # AT/NA AT/NC AT/NG AT/NT AT/NN 
        75299,74774,74274,76174,75130, # AN/AA AN/AC AN/AG AN/AT AN/AN 
        76324,75149,74824,76899,75799, # AN/CA AN/CC AN/CG AN/CT AN/CN 
        74824,73999,75124,75749,74924, # AN/GA AN/GC AN/GG AN/GT AN/GN 
        -350,-1675,-3650,-2225,-1975, # AN/TA AN/TC AN/TG AN/TT AN/TN 
        56524,55561,55143,56649,55969, # AN/NA AN/NC AN/NG AN/NT AN/NN 
        99999,99999,99999,3400,75849, # CA/AA CA/AC CA/AG CA/AT CA/AN 
        99999,99999,99999,6100,76524, # CA/CA CA/CC CA/CG CA/CT CA/CN 
        -900,1900,-700,-8500,-2050, # CA/GA CA/GC CA/GG CA/GT CA/GN 
        99999,99999,99999,1000,75249, # CA/TA CA/TC CA/TG CA/TT CA/TN 
        74774,75474,74824,500,56393, # CA/NA CA/NC CA/NG CA/NT CA/NN 
        99999,99999,5200,99999,76299, # CC/AA CC/AC CC/AG CC/AT CC/AN 
        99999,99999,3600,99999,75899, # CC/CA CC/CC CC/CG CC/CT CC/CN 
        600,-1500,-8000,-800,-2425, # CC/GA CC/GC CC/GG CC/GT CC/GN 
        99999,99999,5200,99999,76299, # CC/TA CC/TC CC/TG CC/TT CC/TN 
        75149,74624,1500,74799,56518, # CC/NA CC/NC CC/NG CC/NT CC/NN 
        99999,1900,99999,99999,75474, # CG/AA CG/AC CG/AG CG/AT CG/AN 
        99999,-1500,99999,99999,74624, # CG/CA CG/CC CG/CG CG/CT CG/CN 
        -4000,-10600,-4900,-4100,-5900, # CG/GA CG/GC CG/GG CG/GT CG/GN 
        99999,-1500,99999,99999,74624, # CG/TA CG/TC CG/TG CG/TT CG/TN 
        73999,-2925,73774,73974,54705, # CG/NA CG/NC CG/NG CG/NT CG/NN 
        99999,99999,99999,99999,99999, # CT/AA CT/AC CT/AG CT/AT CT/AN 
        0,99999,99999,99999,74999, # CT/CA CT/CC CT/CG CT/CT CT/CN 
        -7800,-1500,-2800,-5000,-4275, # CT/GA CT/GC CT/GG CT/GT CT/GN 
        -1200,99999,99999,99999,74699, # CT/TA CT/TC CT/TG CT/TT CT/TN 
        22749,74624,74299,73749,61355, # CT/NA CT/NC CT/NG CT/NT CT/NN 
        99999,75474,76299,75849,81905, # CN/AA CN/AC CN/AG CN/AT CN/AN 
        74999,74624,75899,76524,75511, # CN/CA CN/CC CN/CG CN/CT CN/CN 
        -3025,-2925,-4100,-4600,-3662, # CN/GA CN/GC CN/GG CN/GT CN/GN 
        74699,74624,76299,75249,75217, # CN/TA CN/TC CN/TG CN/TT CN/TN 
        61667,55449,56099,55755,57242, # CN/NA CN/NC CN/NG CN/NT CN/NN 
        99999,99999,99999,700,75174, # GA/AA GA/AC GA/AG GA/AT GA/AN 
        -2900,5200,-600,-8200,-1625, # GA/CA GA/CC GA/CG GA/CT GA/CN 
        99999,99999,99999,1600,75399, # GA/GA GA/GC GA/GG GA/GT GA/GN 
        99999,99999,99999,-1300,74674, # GA/TA GA/TC GA/TG GA/TT GA/TN 
        74274,76299,74849,-1800,55905, # GA/NA GA/NC GA/NG GA/NT GA/NN 
        99999,99999,-600,99999,74849, # GC/AA GC/AC GC/AG GC/AT GC/AN 
        -700,3600,-9800,2300,-1150, # GC/CA GC/CC GC/CG GC/CT GC/CN 
        99999,99999,-6000,99999,73499, # GC/GA GC/GC GC/GG GC/GT GC/GN 
        99999,99999,-4400,99999,73899, # GC/TA GC/TC GC/TG GC/TT GC/TN 
        74824,75899,-5200,75574,55274, # GC/NA GC/NC GC/NG GC/NT GC/NN 
        99999,-700,99999,99999,74824, # GG/AA GG/AC GG/AG GG/AT GG/AN 
        500,-8000,-6000,3300,-2550, # GG/CA GG/CC GG/CG GG/CT GG/CN 
        99999,-4900,99999,99999,73774, # GG/GA GG/GC GG/GG GG/GT GG/GN 
        99999,-2800,99999,5800,50749, # GG/TA GG/TC GG/TG GG/TT GG/TN 
        75124,-4100,73499,52274,49199, # GG/NA GG/NC GG/NG GG/NT GG/NN 
        -600,99999,99999,99999,74849, # GT/AA GT/AC GT/AG GT/AT GT/AN 
        -8400,5200,-4400,-2200,-2450, # GT/CA GT/CC GT/CG GT/CT GT/CN 
        -3100,99999,99999,99999,74224, # GT/GA GT/GC GT/GG GT/GT GT/GN 
        -2500,99999,4100,99999,50399, # GT/TA GT/TC GT/TG GT/TT GT/TN 
        -3650,76299,49924,74449,49255, # GT/NA GT/NC GT/NG GT/NT GT/NN 
        74849,74824,74849,75174,74924, # GN/AA GN/AC GN/AG GN/AT GN/AN 
        -2875,1500,-5200,-1200,-1943, # GN/CA GN/CC GN/CG GN/CT GN/CN 
        74224,73774,73499,75399,74224, # GN/GA GN/GC GN/GG GN/GT GN/GN 
        74374,74299,49924,51124,62430, # GN/TA GN/TC GN/TG GN/TT GN/TN 
        55143,56099,48268,50124,52408, # GN/NA GN/NC GN/NG GN/NT GN/NN 
        4700,3400,700,-7200,400, # TA/AA TA/AC TA/AG TA/AT TA/AN 
        99999,99999,99999,1200,75299, # TA/CA TA/CC TA/CG TA/CT TA/CN 
        99999,99999,99999,-100,74974, # TA/GA TA/GC TA/GG TA/GT TA/GN 
        99999,99999,-100,200,50024, # TA/TA TA/TC TA/TG TA/TT TA/TN 
        76174,75849,50149,-1475,50174, # TA/NA TA/NC TA/NG TA/NT TA/NN 
        7600,6100,-8200,1200,1675, # TC/AA TC/AC TC/AG TC/AT TC/AN 
        99999,99999,2300,99999,75574, # TC/CA TC/CC TC/CG TC/CT TC/CN 
        99999,99999,3300,99999,75824, # TC/GA TC/GC TC/GG TC/GT TC/GN 
        99999,99999,-2200,99999,74449, # TC/TA TC/TC TC/TG TC/TT TC/TN 
        76899,76524,-1200,75299,56880, # TC/NA TC/NC TC/NG TC/NT TC/NN 
        3000,-8500,1600,-100,-1000, # TG/AA TG/AC TG/AG TG/AT TG/AN 
        99999,-800,99999,99999,74799, # TG/CA TG/CC TG/CG TG/CT TG/CN 
        99999,-4100,99999,-1400,48624, # TG/GA TG/GC TG/GG TG/GT TG/GN 
        99999,-5000,99999,99999,73749, # TG/TA TG/TC TG/TG TG/TT TG/TN 
        75749,-4600,75399,49624,49043, # TG/NA TG/NC TG/NG TG/NT TG/NN 
        -7900,1000,-1300,200,-2000, # TT/AA TT/AC TT/AG TT/AT TT/AN 
        700,99999,99999,99999,75174, # TT/CA TT/CC TT/CG TT/CT TT/CN 
        1000,99999,5800,99999,51699, # TT/GA TT/GC TT/GG TT/GT TT/GN 
        -2700,99999,99999,99999,74324, # TT/TA TT/TC TT/TG TT/TT TT/TN 
        -2225,75249,51124,75049,49799, # TT/NA TT/NC TT/NG TT/NT TT/NN 
        1850,500,-1800,-1475,-231, # TN/AA TN/AC TN/AG TN/AT TN/AN 
        75174,74799,75574,75299,75211, # TN/CA TN/CC TN/CG TN/CT TN/CN 
        75249,73974,52274,49624,62780, # TN/GA TN/GC TN/GG TN/GT TN/GN 
        74324,73749,49424,75049,68136, # TN/TA TN/TC TN/TG TN/TT TN/TN 
        56649,55755,43868,49624,51474, # TN/NA TN/NC TN/NG TN/NT TN/NN 
        76174,75849,75174,400,56899, # NA/AA NA/AC NA/AG NA/AT NA/AN 
        74274,76299,74849,1675,56774, # NA/CA NA/CC NA/CG NA/CT NA/CN 
        74774,75474,74824,-1000,56018, # NA/GA NA/GC NA/GG NA/GT NA/GN 
        75299,75574,49824,-2000,49674, # NA/TA NA/TC NA/TG NA/TT NA/TN 
        75130,75799,68667,-231,54841, # NA/NA NA/NC NA/NG NA/NT NA/NN 
        76899,76524,-1625,75299,56774, # NC/AA NC/AC NC/AG NC/AT NC/AN 
        74824,75899,-1150,75574,56286, # NC/CA NC/CC NC/CG NC/CT NC/CN 
        75149,74624,-2550,74799,55505, # NC/GA NC/GC NC/GG NC/GT NC/GN 
        76324,74999,-2450,75174,56011, # NC/TA NC/TC NC/TG NC/TT NC/TN 
        75799,75511,-1943,75211,56144, # NC/NA NC/NC NC/NG NC/NT NC/NN 
        75749,-2050,75399,74974,56018, # NG/AA NG/AC NG/AG NG/AT NG/AN 
        75124,-2425,73499,75824,55505, # NG/CA NG/CC NG/CG NG/CT NG/CN 
        73999,-5900,73774,48624,47624, # NG/GA NG/GC NG/GG NG/GT NG/GN 
        74824,-4275,74224,51699,49118, # NG/TA NG/TC NG/TG NG/TT NG/TN 
        74924,-3662,74224,62780,52066, # NG/NA NG/NC NG/NG NG/NT NG/NN 
        23174,75249,74674,75049,62036, # NT/AA NT/AC NT/AG NT/AT NT/AN 
        -600,76299,73899,74449,56011, # NT/CA NT/CC NT/CG NT/CT NT/CN 
        -2650,74624,50749,73749,49118, # NT/GA NT/GC NT/GG NT/GT NT/GN 
        -3400,74699,50399,74324,49005, # NT/TA NT/TC NT/TG NT/TT NT/TN 
        4131,75217,62430,74392,54042, # NT/NA NT/NC NT/NG NT/NT NT/NN 
        62999,56393,55905,56430,57932, # NN/AA NN/AC NN/AG NN/AT NN/AN 
        55905,56518,55274,56880,56144, # NN/CA NN/CC NN/CG NN/CT NN/CN 
        55318,54705,49199,49043,52066, # NN/GA NN/GC NN/GG NN/GT NN/GN 
        55761,55249,42999,49799,50952, # NN/TA NN/TC NN/TG NN/TT NN/TN 
        57495,55716,50844,53038,54273  # NN/NA NN/NC NN/NG NN/NT NN/NN 
        ]

        
        
        # End of tables nearest-neighbor parameter ------------------------------
        
        # init values
        self.seq = seq
        self.template = template
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
    
    @staticmethod
    def complement(seq:str) -> str:
        trantab = str.maketrans('ACGTN', 'TGCAN')
        return seq.upper().translate(trantab)
    
    @staticmethod
    def is_complement(seq1:str, seq2:str):
        return Primer.complement(seq1) == seq2.upper()

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

    def calc_match(self):
        ...
    
    
    def calc_mismatch(self):
        ...
    
    def calc_gap(self):
        ...

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
            two_mer_primer = self.seq[i-1] + s
            two_mer_tmp = self.template[i-1] + s
            # match 
            if Primer.is_complement(two_mer_primer, two_mer_tmp):
                d_index = self.base2int(two_mer_primer)
                self.dh += -100 * self.dH_m[d_index]
                self.ds += -0.1 * self.dS_m[d_index]
            else:  ## mismatch and gap
                # two gap
                if two_mer_primer == "--" or two_mer_tmp == "--":
                    self.ds -= 9.35
                elif two_mer_primer.find("-") > -1 or two_mer_tmp.find("-") > -1: # one gap
                    ...
                else: # mismatch
                    ...
                    
                    
            
        
        # init value and salt corrections and calculate Tm finally
        # self.dh *= -100
        # self.ds *= -0.1

        self.GC_count = 0 if self.formamide_conc == 0.0 else str.count(self.seq,"C") + str.count(self.seq,"G")
        self.K_mM += self.divalent_to_monovalent()
        self.ds = self.ds + 0.368 * (len(self.seq) - 1) * log(self.K_mM / 1000.0 )

        self.Tm = self.dh / (self.ds + 1.987 * log(self.DNA_nM / self.base)) - self.T_KELVIN
        self.Tm -= self.dmso_conc * self.dmso_fact
        self.Tm += (0.453 * self.GC_count / len(self.seq) - 2.88) * self.formamide_conc
        return self.Tm


if __name__=="__main__":
    import primer3
    primer1 = Primer("AAAAAAAAAA")
    print(primer1.Tm)
    print(primer3.calc_tm("AAAAAAAAAA"))
    primer2 = Primer("GGCCGGAGTAAGCTGACAT")
    print(primer2.Tm)
    print(primer3.calc_tm("GGCCGGAGTAAGCTGACAT"))
from collections import defaultdict, deque


int_loop_energy_dict = defaultdict(deque)

int_loop_energy_dict["21dg"] = deque([
		# A  A
		# T AT
		[
			31.8,31.8,29.6,31.8,
			29.3,29.3,29.3,29.3,
			20.7,20.7,20.7,20.7,
			29.3,29.3,29.3,14.1,
		],
		# A  C
		# T AG
		[
			28.9,29.8,27.8,29.8,
			23.0,23.0,23.0,20.2,
			23.4,23.4,13.9,23.4,
			23.0,23.0,23.0,19.4,
		],
		# A  G
		# T AC
		[
			23.9,23.9,12.6,23.9,
			24.2,24.2,24.2,14.0,
			16.7,16.7,16.7,16.7,
			24.2,6.6,24.2,11.4,
		],
		# A  T
		# T AA
		[
			23.8,28.7,26.0,23.8,
			28.7,28.7,28.7,23.2,
			26.0,26.0,26.0,26.0,
			28.7,23.2,28.7,28.7,
		],
		# A  A
		# T CT
		[
			31.8,20.4,29.6,20.4,
			29.3,33.6,29.3,17.3,
			20.7,20.4,20.5,20.4,
			29.3,22.4,29.3,14.1,
		],
		# A  C
		# T CG
		[
			29.8,29.8,27.8,29.8,
			23.0,31.7,23.0,20.2,
			23.4,29.8,13.9,29.8,
			23.0,30.5,23.0,19.4,
		],
		# A  G
		# T CC
		[
			22.8,22.8,12.6,22.8,
			24.2,21.5,24.2,14.0,
			16.7,22.8,9.0,22.8,
			6.6,6.6,6.6,6.6,
		],
		# A  T
		# T CA
		[
			28.7,28.7,26.0,28.7,
			28.7,26.4,28.7,23.2,
			26.0,28.7,23.0,28.7,
			23.2,23.2,23.2,23.2,
		],
		# A  A
		# T GT
		[
			29.6,29.6,29.6,29.6,
			29.3,29.3,29.3,29.3,
			20.7,20.5,20.5,20.5,
			29.3,29.3,29.3,14.1,
		],
		# A  C
		# T GG
		[
			27.8,27.8,27.8,27.8,
			23.0,23.0,23.0,20.2,
			13.9,13.9,13.9,13.9,
			23.0,23.0,23.0,19.4,
		],
		# A  G
		# T GC
		[
			12.6,12.6,12.6,12.6,
			24.2,24.2,24.2,14.0,
			16.7,9.0,9.0,9.0,
			24.2,6.6,24.2,11.4,
		],
		# A  T
		# T GA
		[
			26.0,26.0,26.0,26.0,
			28.7,28.7,28.7,23.2,
			26.0,23.0,23.0,23.0,
			28.7,23.2,28.7,28.7,
		],
		# A  A
		# T TT
		[
			31.8,20.4,29.6,20.4,
			17.3,17.3,17.3,17.3,
			20.7,20.4,20.5,20.4,
			14.1,14.1,14.1,14.1,
		],
		# A  C
		# T TG
		[
			29.8,29.8,27.8,29.8,
			20.2,20.2,20.2,20.2,
			23.4,29.8,13.9,29.8,
			19.4,19.4,19.4,19.4,
		],
		# A  G
		# T TC
		[
			22.8,22.8,12.6,22.8,
			14.0,14.0,14.0,14.0,
			16.7,22.8,9.0,22.8,
			11.4,6.6,11.4,11.4,
		],
		# A  T
		# T TA
		[
			28.7,28.7,26.0,28.7,
			23.2,23.2,23.2,23.2,
			26.0,28.7,23.0,28.7,
			28.7,23.2,28.7,24.6,
		],
		# C  A
		# G AT
		[
			28.9,28.9,30.2,28.9,
			28.5,28.5,28.5,19.0,
			25.8,25.8,17.8,25.8,
			28.5,21.9,28.5,19.6,
		],
		# C  C
		# G AG
		[
			18.3,23.5,20.2,23.5,
			28.2,28.2,28.2,20.8,
			22.8,22.8,18.0,22.8,
			28.2,19.4,28.2,21.0,
		],
		# C  G
		# G AC
		[
			22.9,22.9,25.0,22.9,
			24.0,24.0,24.0,20.5,
			25.0,25.0,25.0,25.0,
			24.0,20.5,24.0,16.2,
		],
		# C  T
		# G AA
		[
			23.9,24.2,16.7,24.2,
			22.8,22.8,22.8,6.6,
			12.6,12.6,9.0,12.6,
			22.8,14.0,22.8,11.4,
		],
		# C  A
		# G CT
		[
			28.9,32.1,30.2,32.1,
			28.5,24.8,28.5,19.0,
			25.8,32.1,17.8,32.1,
			21.9,21.9,21.9,19.6,
		],
		# C  C
		# G CG
		[
			23.5,23.5,20.2,23.5,
			28.2,26.4,28.2,20.8,
			22.8,23.5,18.0,23.5,
			19.4,19.4,19.4,19.4,
		],
		# C  G
		# G CC
		[
			22.9,24.0,25.0,24.0,
			24.0,23.0,24.0,20.5,
			25.0,24.0,19.0,24.0,
			20.5,20.5,20.5,16.2,
		],
		# C  T
		# G CA
		[
			24.2,24.2,16.7,24.2,
			22.8,21.5,22.8,6.6,
			12.6,24.2,9.0,24.2,
			14.0,14.0,14.0,11.4,
		],
		# C  A
		# G GT
		[
			30.2,30.2,30.2,30.2,
			28.5,28.5,28.5,19.0,
			17.8,17.8,17.8,17.8,
			28.5,21.9,28.5,19.6,
		],
		# C  C
		# G GG
		[
			20.2,20.2,20.2,20.2,
			28.2,28.2,28.2,20.8,
			18.0,18.0,18.0,18.0,
			28.2,19.4,28.2,21.0,
		],
		# C  G
		# G GC
		[
			25.0,25.0,25.0,25.0,
			24.0,24.0,24.0,20.5,
			25.0,19.0,19.0,19.0,
			24.0,20.5,24.0,16.2,
		],
		# C  T
		# G GA
		[
			16.7,16.7,16.7,16.7,
			22.8,22.8,22.8,6.6,
			9.0,9.0,9.0,9.0,
			22.8,14.0,22.8,11.4,
		],
		# C  A
		# G TT
		[
			28.9,32.1,30.2,32.1,
			19.0,19.0,19.0,19.0,
			25.8,32.1,17.8,32.1,
			19.6,19.6,19.6,19.6,
		],
		# C  C
		# G TG
		[
			23.5,23.5,20.2,23.5,
			20.8,20.8,20.8,20.8,
			22.8,23.5,18.0,23.5,
			21.0,21.0,21.0,21.0,
		],
		# C  G
		# G TC
		[
			22.9,24.0,25.0,24.0,
			20.5,20.5,20.5,20.5,
			25.0,24.0,19.0,24.0,
			16.2,16.2,16.2,16.2,
		],
		# C  T
		# G TA
		[
			24.2,24.2,16.7,24.2,
			6.6,6.6,6.6,6.6,
			12.6,24.2,9.0,24.2,
			11.4,11.4,11.4,11.4,
		],
		# G  A
		# C AT
		[
			17.5,17.5,17.5,17.5,
			23.1,23.1,23.1,23.1,
			25.8,25.8,25.8,25.8,
			23.1,23.1,23.1,8.9,
		],
		# G  C
		# C AG
		[
			21.5,21.5,21.3,21.5,
			22.8,22.8,22.8,22.8,
			21.3,21.3,16.1,21.3,
			22.8,22.8,22.8,22.8,
		],
		# G  G
		# C AC
		[
			18.3,18.3,22.8,18.3,
			23.5,23.5,23.5,23.5,
			20.2,20.2,20.2,20.2,
			23.5,23.5,23.5,23.5,
		],
		# G  T
		# C AA
		[
			28.9,23.0,23.4,23.0,
			29.8,29.8,29.8,29.8,
			27.8,27.8,13.9,27.8,
			29.8,20.2,29.8,19.4,
		],
		# G  A
		# C CT
		[
			17.5,32.1,30.2,32.1,
			23.1,31.8,23.1,18.5,
			25.8,32.1,22.3,32.1,
			23.1,22.6,23.1,8.9,
		],
		# G  C
		# C CG
		[
			21.5,22.8,21.3,22.8,
			22.8,24.9,22.8,18.3,
			21.3,22.8,16.1,22.8,
			22.8,18.3,22.8,18.3,
		],
		# G  G
		# C CC
		[
			18.3,28.2,22.8,28.2,
			23.5,26.4,23.5,19.4,
			20.2,28.2,18.0,28.2,
			23.5,20.8,23.5,21.0,
		],
		# G  T
		# C CA
		[
			23.0,23.0,23.4,23.0,
			29.8,31.7,29.8,30.5,
			27.8,23.0,13.9,23.0,
			20.2,20.2,20.2,20.2,
		],
		# G  A
		# C GT
		[
			30.2,30.2,30.2,30.2,
			23.1,23.1,23.1,23.1,
			25.8,22.3,22.3,22.3,
			23.1,23.1,23.1,8.9,
		],
		# G  C
		# C GG
		[
			21.3,21.3,21.3,21.3,
			22.8,22.8,22.8,22.8,
			16.1,16.1,16.1,16.1,
			22.8,22.8,22.8,22.8,
		],
		# G  G
		# C GC
		[
			22.8,22.8,22.8,22.8,
			23.5,23.5,23.5,23.5,
			20.2,18.0,18.0,18.0,
			23.5,23.5,23.5,23.5,
		],
		# G  T
		# C GA
		[
			23.4,23.4,23.4,23.4,
			29.8,29.8,29.8,29.8,
			13.9,13.9,13.9,13.9,
			29.8,20.2,29.8,19.4,
		],
		# G  A
		# C TT
		[
			17.5,32.1,30.2,32.1,
			23.1,18.5,23.1,18.5,
			25.8,32.1,22.3,32.1,
			8.9,8.9,8.9,8.9,
		],
		# G  C
		# C TG
		[
			21.5,22.8,21.3,22.8,
			22.8,18.3,22.8,18.3,
			21.3,22.8,16.1,22.8,
			22.8,18.3,22.8,19.2,
		],
		# G  G
		# C TC
		[
			18.3,28.2,22.8,28.2,
			23.5,19.4,23.5,19.4,
			20.2,28.2,18.0,28.2,
			23.5,21.0,23.5,21.0,
		],
		# G  T
		# C TA
		[
			23.0,23.0,23.4,23.0,
			29.8,30.5,29.8,30.5,
			27.8,23.0,13.9,23.0,
			19.4,19.4,19.4,19.4,
		],
		# T  A
		# A AT
		[
			31.0,31.0,30.8,31.0,
			33.2,33.2,33.2,21.4,
			30.8,30.8,30.8,30.8,
			33.2,21.4,33.2,25.1,
		],
		# T  C
		# A AG
		[
			17.5,17.5,25.8,17.5,
			32.1,32.1,32.1,32.1,
			30.2,30.2,22.3,30.2,
			32.1,18.5,32.1,8.9,
		],
		# T  G
		# A AC
		[
			28.9,28.9,25.8,28.9,
			32.1,32.1,32.1,21.9,
			30.2,30.2,17.8,30.2,
			32.1,19.0,32.1,19.6,
		],
		# T  T
		# A AA
		[
			31.8,29.3,20.7,29.3,
			20.4,20.4,20.4,22.4,
			29.6,29.6,20.5,29.6,
			20.4,17.3,20.4,14.1,
		],
		# T  A
		# A CT
		[
			31.0,33.2,30.8,33.2,
			33.2,36.3,33.2,21.4,
			30.8,33.2,27.3,33.2,
			21.4,21.4,21.4,21.4,
		],
		# T  C
		# A CG
		[
			17.5,23.1,25.8,23.1,
			32.1,31.8,32.1,31.8,
			30.2,23.1,22.3,23.1,
			18.5,18.5,18.5,8.9,
		],
		# T  G
		# A CC
		[
			28.5,28.5,25.8,28.5,
			32.1,24.8,32.1,21.9,
			30.2,28.5,17.8,28.5,
			19.0,19.0,19.0,19.6,
		],
		# T  T
		# A CA
		[
			29.3,29.3,20.7,29.3,
			20.4,33.6,20.4,22.4,
			29.6,29.3,20.5,29.3,
			17.3,17.3,17.3,14.1,
		],
		# T  A
		# A GT
		[
			30.8,30.8,30.8,30.8,
			33.2,33.2,33.2,21.4,
			30.8,27.3,27.3,27.3,
			33.2,21.4,33.2,25.1,
		],
		# T  C
		# A GG
		[
			25.8,25.8,25.8,25.8,
			32.1,32.1,32.1,32.1,
			22.3,22.3,22.3,22.3,
			32.1,18.5,32.1,8.9,
		],
		# T  G
		# A GC
		[
			25.8,25.8,25.8,25.8,
			32.1,32.1,32.1,21.9,
			17.8,17.8,17.8,17.8,
			32.1,19.0,32.1,19.6,
		],
		# T  T
		# A GA
		[
			20.7,20.7,20.7,20.7,
			20.4,20.4,20.4,22.4,
			20.5,20.5,20.5,20.5,
			20.4,17.3,20.4,14.1,
		],
		# T  A
		# A TT
		[
			31.0,33.2,30.8,33.2,
			21.4,21.4,21.4,21.4,
			30.8,33.2,27.3,33.2,
			25.1,21.4,25.1,25.1,
		],
		# T  C
		# A TG
		[
			17.5,23.1,25.8,23.1,
			32.1,22.6,32.1,22.6,
			30.2,23.1,22.3,23.1,
			8.9,8.9,8.9,8.9,
		],
		# T  G
		# A TC
		[
			28.5,28.5,25.8,28.5,
			21.9,21.9,21.9,21.9,
			30.2,28.5,17.8,28.5,
			19.6,19.6,19.6,19.6,
		],
		# T  T
		# A TA
		[
			29.3,29.3,20.7,29.3,
			22.4,22.4,22.4,22.4,
			29.6,29.3,20.5,29.3,
			14.1,14.1,14.1,14.1,
		],
])



# matrix construct

matrix = int_loop_energy_dict["21dg"]

aan_averages = []
acn_averages = []
agn_averages = []
atn_averages = []

ana_averages = []
anc_averages = []
ang_averages = []
ant_averages = []
ann_averages = []

can_averages = []
ccn_averages = []
cgn_averages = []
ctn_averages = []

cna_averages = []
cnc_averages = []
cng_averages = []
cnt_averages = []
cnn_averages = []

gan_averages = []
gcn_averages = []
ggn_averages = []
gtn_averages = []

gna_averages = []
gnc_averages = []
gng_averages = []
gnt_averages = []
gnn_averages = []

tan_averages = []
tcn_averages = []
tgn_averages = []
ttn_averages = []

tna_averages = []
tnc_averages = []
tng_averages = []
tnt_averages = []
tnn_averages = []

naa_averages = []
nac_averages = []
nag_averages = []
nat_averages = []
nan_averages = []

nca_averages = []
ncc_averages = []
ncg_averages = []
nct_averages = []
ncn_averages = []

nga_averages = []
ngc_averages = []
ngg_averages = []
ngt_averages = []
ngn_averages = []

nta_averages = []
ntc_averages = []
ntg_averages = []
ntt_averages = []
ntn_averages = []

nna_averages = []
nnc_averages = []
nng_averages = []
nnt_averages = []
nnn_averages = []

num_columns = len(matrix[0])


start = -1
end = 3

def append_avg(avg_li: list) -> None:
    global start, end
    avg_li.append(round(sum(row[col] for i, row in enumerate(matrix) if start < i <= end) / 4,1))  # 将平均值添加到结果列表
    start += 4
    end += 4
# print(num_columns)
def base2int(base: str) -> int:
    """将base(without N)转化为索引号"""
    trantab = str.maketrans("ACGTN", "01234")
    return int(base.upper().translate(trantab), base=5)
def index_intl21(bases: str) -> int:
    """            
            X
           A  A  ---- AAA
           T  T
            YA
    """
    upstream_base, y_neibor_base, downstream_base = bases[:]
    external_index = base2int(upstream_base + y_neibor_base + downstream_base)
    return external_index
# print(index_intl21("TNN"))

def obtain_N_sub_idx(bases: str) -> list:
    if index_intl21(bases) < 100 :
        sub_li = [bases[0] + x + bases[2]  for x in ["A", "C", "G", "T"]]
        # print([index_intl21(x) for x in sub_li])  
        return  [index_intl21(x) for x in sub_li]
    else:
        sub_li = [ x + bases[1] + bases[2]  for x in ["A", "C", "G", "T"]]
        return  [index_intl21(x) for x in sub_li]


    
def obtain_mean(bases: str,list2append: list) -> list:
    
    idx_li = obtain_N_sub_idx(bases)
    
    for i in range(16):
        ele = round(sum([matrix[x][i] for x in idx_li]) / 4 , 1)
        list2append.append(ele)

for col in range(num_columns):

    # AAN
    append_avg(aan_averages)
    # ACN
    append_avg(acn_averages)
    # AGN
    append_avg(agn_averages)
    # ATN
    append_avg(atn_averages)
    
    # ANA
    

    # CAN
    append_avg(can_averages)
    # CCN
    append_avg(ccn_averages)
    # CGN
    append_avg(cgn_averages)
    # CTN
    append_avg(ctn_averages)

    # GAN
    append_avg(gan_averages)
    # GCN
    append_avg(gcn_averages)
    # GGN
    append_avg(ggn_averages)
    # GTN
    append_avg(gtn_averages)

    # TAN
    append_avg(tan_averages)
    # TCN
    append_avg(tcn_averages)
    # TGN
    append_avg(tgn_averages)
    # TTN
    append_avg(ttn_averages)

    
    start = -1
    end = 3
# insert matrix
insert_idx = 4
# AAN
int_loop_energy_dict["21dg"].insert(insert_idx,aan_averages)
insert_idx += 5
# ACN
int_loop_energy_dict["21dg"].insert(insert_idx,acn_averages)
insert_idx += 5
# AGN
int_loop_energy_dict["21dg"].insert(insert_idx,agn_averages)
insert_idx += 5
# ATN
int_loop_energy_dict["21dg"].insert(insert_idx,atn_averages)
insert_idx += 1
# ANA
obtain_mean("ANA", ana_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,ana_averages)
insert_idx += 1
# ANC
obtain_mean("ANC", anc_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,anc_averages)
insert_idx += 1
# ANG
obtain_mean("ANG", ang_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,ang_averages)
insert_idx += 1
# ANT
obtain_mean("ANT", ant_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,ant_averages)
insert_idx += 1
# ANN
obtain_mean("ANN", ann_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,ann_averages)
insert_idx += 5

# CAN
int_loop_energy_dict["21dg"].insert(insert_idx,can_averages)
insert_idx += 5
# CCN
int_loop_energy_dict["21dg"].insert(insert_idx,ccn_averages)
insert_idx += 5
# CGN
int_loop_energy_dict["21dg"].insert(insert_idx,cgn_averages)
insert_idx += 5
# CTN
int_loop_energy_dict["21dg"].insert(insert_idx,ctn_averages)
insert_idx += 1
# CNA
obtain_mean("CNA", cna_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,cna_averages)
insert_idx += 1
# CNC
obtain_mean("CNC", cnc_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,cnc_averages)
insert_idx += 1
# CNG
obtain_mean("CNG", cng_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,cng_averages)
insert_idx += 1
# CNT
obtain_mean("CNT", cnt_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,cnt_averages)
insert_idx += 1
# CNT
obtain_mean("CNN", cnn_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,cnn_averages)
insert_idx += 5

# GAN
int_loop_energy_dict["21dg"].insert(insert_idx,gan_averages)
insert_idx += 5
# GCN
int_loop_energy_dict["21dg"].insert(insert_idx,gcn_averages)
insert_idx += 5
# GGN
int_loop_energy_dict["21dg"].insert(insert_idx,ggn_averages)
insert_idx += 5
# GTN
int_loop_energy_dict["21dg"].insert(insert_idx,gtn_averages)
insert_idx += 1
# GNA
obtain_mean("GNA", gna_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,gna_averages)
insert_idx += 1
# GNC
obtain_mean("GNC", gnc_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,gnc_averages)
insert_idx += 1
# GNG
obtain_mean("GNG", gng_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,gng_averages)
insert_idx += 1
# GNT
obtain_mean("GNT", gnt_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,gnt_averages)
insert_idx += 1
# GNT
obtain_mean("GNN", gnn_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,gnn_averages)
insert_idx += 5


# TAN
int_loop_energy_dict["21dg"].insert(insert_idx,tan_averages)
insert_idx += 5
# TCN
int_loop_energy_dict["21dg"].insert(insert_idx,tcn_averages)
insert_idx += 5
# TGN
int_loop_energy_dict["21dg"].insert(insert_idx,tgn_averages)
insert_idx += 5
# TTN
int_loop_energy_dict["21dg"].insert(insert_idx,ttn_averages)
insert_idx += 1
# TNA
obtain_mean("TNA", tna_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,tna_averages)
insert_idx += 1
# TNC
obtain_mean("TNC", tnc_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,tnc_averages)
insert_idx += 1
# TNG
obtain_mean("TNG", tng_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,tng_averages)
insert_idx += 1
# TNT
obtain_mean("TNT", tnt_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,tnt_averages)
insert_idx += 1
# TNT
obtain_mean("TNN", tnn_averages)
int_loop_energy_dict["21dg"].insert(insert_idx,tnn_averages)
insert_idx += 5



# fill N matrix
n_dict = defaultdict(list)
basen_li = ["A", "C", "G", "T", "N"]
for x in basen_li:
    for y in basen_li:
        idx_li = obtain_N_sub_idx(f"N{x}{y}")
        for i in range(16):
            ele = round(sum([matrix[x][i] for x in idx_li]) / 4 , 1)
            n_dict[f"N{x}{y}"].append(ele)
        int_loop_energy_dict["21dg"].append(n_dict[f"N{x}{y}"])

# print(len(int_loop_energy_dict["21dg"]))   

# format print dict
for i, x in enumerate(int_loop_energy_dict["21dg"]):
    first_base = basen_li[i // 25]
    second_base = basen_li[i % 25 // 5]
    third_base = basen_li[i % 25 % 5]
    print(f"\t\t# {first_base} {second_base} {third_base}")
    
    print("\t\t[")
    print("\t\t\t",end="")
    print(", ".join([str(y) for y in x]))
    print("\t\t],")
    
    


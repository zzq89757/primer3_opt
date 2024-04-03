# Define parameters
S_A_A = 240
S_A_C = 173
S_A_G = 208
S_A_T = 239
S_A_N = 215

S_C_A = 129
S_C_C = 266
S_C_G = 278
S_C_T = 208
S_C_N = 220

S_G_A = 135
S_G_C = 267
S_G_G = 266
S_G_T = 173
S_G_N = 210

S_T_A = 169
S_T_C = 135
S_T_G = 129
S_T_T = 240
S_T_N = 168

S_N_A = 168
S_N_C = 210
S_N_G = 220
S_N_T = 215
S_N_N = 203

H_A_A = 91
H_A_C = 65
H_A_G = 78
H_A_T = 86
H_A_N = 80

H_C_A = 58
H_C_C = 110
H_C_G = 119
H_C_T = 78
H_C_N = 91

H_G_A = 56
H_G_C = 111
H_G_G = 110
H_G_T = 65
H_G_N = 85

H_T_A = 60
H_T_C = 56
H_T_G = 58
H_T_T = 91
H_T_N = 66

H_N_A = 66
H_N_C = 85
H_N_G = 91
H_N_T = 80
H_N_N = 80

# Delta G's of disruption * 1000
G_A_A = 1900
G_A_C = 1300
G_A_G = 1600
G_A_T = 1500
G_A_N = 1575

G_C_A = 1900
G_C_C = 3100
G_C_G = 3600
G_C_T = 1600
G_C_N = 2550

G_G_A = 1600
G_G_C = 3100
G_G_G = 3100
G_G_T = 1300
G_G_N = 2275

G_T_A = 900
G_T_C = 1600
G_T_G = 1900
G_T_T = 1900
G_T_N = 1575

G_N_A = 1575
G_N_C = 2275
G_N_G = 2550
G_N_T = 1575
G_N_N = 1994

# New parameters
DS_A_A = 222
DS_A_C = 224
DS_A_G = 210
DS_A_T = 204
DS_A_N = 224

DS_C_A = 227
DS_C_C = 199
DS_C_G = 272
DS_C_T = 210
DS_C_N = 272

DS_G_A = 222
DS_G_C = 244
DS_G_G = 199
DS_G_T = 224
DS_G_N = 244

DS_T_A = 213
DS_T_C = 222
DS_T_G = 227
DS_T_T = 222
DS_T_N = 227

DS_N_A = 168
DS_N_C = 210
DS_N_G = 220
DS_N_T = 215
DS_N_N = 220

# Parameters for disruption
DH_A_A = 79
DH_A_C = 84
DH_A_G = 78
DH_A_T = 72
DH_A_N = 72

DH_C_A = 85
DH_C_C = 80
DH_C_G = 106
DH_C_T = 78
DH_C_N = 78

DH_G_A = 82
DH_G_C = 98
DH_G_G = 80
DH_G_T = 84
DH_G_N = 80

DH_T_A = 72
DH_T_C = 82
DH_T_G = 85
DH_T_T = 79
DH_T_N = 72

DH_N_A = 72
DH_N_C = 80
DH_N_G = 78
DH_N_T = 72
DH_N_N = 72

# Delta G values of disruption * 1000
DG_A_A = 1000
DG_A_C = 1440
DG_A_G = 1280
DG_A_T = 880
DG_A_N = 880

DG_C_A = 1450
DG_C_C = 1840
DG_C_G = 2170
DG_C_T = 1280
DG_C_N = 1450

DG_G_A = 1300
DG_G_C = 2240
DG_G_G = 1840
DG_G_T = 1440
DG_G_N = 1300

DG_T_A = 580
DG_T_C = 1300
DG_T_G = 1450
DG_T_T = 1000
DG_T_N = 580

DG_N_A = 580
DG_N_C = 1300
DG_N_G = 1280
DG_N_T = 880
DG_N_N = 580

# symmetry correction if seq is symmetrical
# if(sym == 1) {
#       ds+=14;
#     }

# Terminal AT penalty 
# if(strncmp("A", s, 1)==0
#        || strncmp("T", s, 1)==0)  {
#       ds += -41;
#       dh += -23;
#     } else if (strncmp("C", s, 1)==0
#                || strncmp("G", s, 1)==0) {
#       ds += 28;
#       dh += -1;
#     }
#     s+=len;
#     if(strncmp("T", s, 1)==0
#        || strncmp("A", s, 1)==0) {
#       ds += -41;
#       dh += -23;
#     } else if (strncmp("C", s, 1)==0
#                || strncmp("G", s, 1)==0) {
#       ds += 28;
#       dh += -1;
#     }

# base judge
# c = *s; s++ if ('A' == c) { dh += 79; ds += 222; goto A_STATE2; } else if ('T' == c) { dh += 72; ds += 204; goto T_STATE2; } else if ('G' == c) { dh += 78; ds += 210; goto G_STATE2; } else if ('C' == c) { dh += 84; ds += 224; goto C_STATE2; } else if ('N' == c) { dh += 72; ds += 224

# NN
#define STATE2(LAST)     \
  #  CATID2(LAST,_STATE2): \
  #  c = *s; s++;         \
  #  DO_PAIR2(LAST,A)      \
  #  else DO_PAIR2(LAST,T) \
  #  else DO_PAIR2(LAST,G) \
  #  else DO_PAIR2(LAST,C) \
  #  else DO_PAIR2(LAST,N) \
# DO_PAIR2(LAST,THIS)          \
#      if (CATID2(THIS,_CHAR) == c) { \
#       dh += CATID5(DH,_,LAST,_,THIS); \
#       ds += CATID5(DS,_,LAST,_,THIS); \
#       goto CATID2(THIS,_STATE2);

# init value
  # delta_H = dh * -100.0;  /*
  #                          * Nearest-neighbor thermodynamic values for dh
  #                          * are given in 100 cal/mol of interaction.
  #                          */
  # delta_S = ds * -0.1;     /*
  #                           * Nearest-neighbor thermodynamic values for ds
  #                           * are in in .1 cal/K per mol of interaction.
  #                           */
  # ret.Tm=0;  /* Melting temperature */
  
  #  // salt_corrections----------------------------------------------------------------------------------
  #    K_mM = K_mM + divalent_to_monovalent(divalent_conc, dntp_conc);
  #    delta_S = delta_S + 0.368 * (len - 1) * log(K_mM / 1000.0 );
  #    if(sym == 1) { /* primer is symmetrical */
  #       /* Equation A */
  #       ret.Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/1000000000.0)) - T_KELVIN;
  #       ret.Tm -= dmso_conc * dmso_fact;
  #       ret.Tm += (0.453 * ((double) GC_count) / len - 2.88) * formamide_conc;
  #       if (annealing_temp > 0.0) {
  #          ddG = delta_H - (annealing_temp + T_KELVIN) * delta_S;
  #          ka = exp(-ddG / (1.987 * (annealing_temp + T_KELVIN)));
  #          ret.bound = (1 / (1 + sqrt(1/((DNA_nM/1000000000.0) * ka)))) * 100;
  #       }
  #    } else {
  #       /* Equation B */
  #       ret.Tm = delta_H / (delta_S + 1.987 * log(DNA_nM/4000000000.0)) - T_KELVIN;
  #       ret.Tm -= dmso_conc * dmso_fact;
  #       ret.Tm += (0.453 * ((double) GC_count) / len - 2.88) * formamide_conc;
  #       if (annealing_temp > 0.0) {
  #          ddG = delta_H - (annealing_temp + T_KELVIN) * delta_S;
  #          ka = exp(-ddG / (1.987 * (annealing_temp + T_KELVIN)));
  #          ret.bound = (1 / (1 + sqrt(1/((DNA_nM/4000000000.0) * ka)))) * 100;
  #       }
  #    }
  
#   double divalent_to_monovalent(double divalent,
#                               double dntp)
# {
#    if(divalent==0) dntp=0;
#    if(divalent<0 || dntp<0) return OLIGOTM_ERROR;
#    if(divalent<dntp)
#      /* According to theory, melting temperature does not depend on
#         divalent cations */
#      divalent=dntp;
#    return 120 * (sqrt(divalent-dntp));
# }
if __name__ == "__main__":
  print(symmetry("ATTT"))
from ClassesBlackOil import *
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
'--------------------'
dg = 0.84
do = 0.86
Pb = 5000  # psia
# P = 3626  # psia
T = 122  # °F
Tsep = 80  # °F
Psep = 100  # psia
'--------------------'

Matrizonha = []
Linha = []
P_P = []
Rs_Rs = []
Bo_Bo = []
Co_Co = []
uo_uo = []
Z_Z = []
rho_g_rho_g = []
Bg_Bg = []
Cg_Cg = []
ug_ug = []
Rho_oleo_Rho_oleo = []

for P in range(14, 7000, 100):
    PVT = BlackOil()
    PVT.do = do
    PVT.T = converte_T_para_R(T, 'f')
    PVT.API = PVT.fase_oleo_grau_API_com_do__API__()
    PVT.dg = dg
    PVT.Mg = PVT.fase_gas_massa_do_gas__Mg__()
    PVT.T = T
    PVT.P = P
    if P <= Pb:
        PVT.T = converte_T_para_F(T, 'f')
        PVT.Rs = PVT.fase_oleo_razao_de_solubilidade_standing_1947_P_menorIgual_Pb__Rs__()
        PVT.Bo = PVT.fase_oleo_fator_volume_formacao_de_oleo_standing_1947_P_menorIgual_Pb__Bo__()
    else:
        PVT.P = Pb
        PVT.Rs = PVT.fase_oleo_razao_de_solubilidade_standing_1947_P_menorIgual_Pb__Rs__()
        PVT.T = converte_T_para_F(122, 'f')
        PVT.Co = PVT.fase_oleo_ompressibilidade_isotermica_oleo_petrosky_e_farshad_1993_P_maiorIgual_Pb__Co__()
        PVT.T = converte_T_para_F(122, 'f')
        PVT.Bob = PVT.fase_oleo_fator_volume_formacao_de_oleo_standing_1947_P_menorIgual_Pb__Bo__()
        PVT.P = P
        PVT.Pb = Pb
        PVT.Bo = PVT.fase_oleo_fator_volume_formacao_de_oleo_P_maior_Pb__Bo__()
    PVT.Pb = Pb

    "---"
    if P <= Pb:
        PVT.P = P
        PVT.Rho_oleo = PVT.fase_oleo_massa_especifica_oleo__Rho_oleo__()
    else:
        PVT.P = Pb
        PVT.Rho_ob = PVT.fase_oleo_massa_especifica_oleo__Rho_oleo__()
        PVT.Rho_oleo = PVT.fase_oleo_massa_especifica_oleo__Rho_oleo__()
    "---"
    PVT.T = converte_T_para_R(T, 'f')
    PVT.Ppc = PVT.fase_gas_pressao_temperatura_pseudocritica__Ppr__Tpr__Ppc__Tpc__()[2]
    PVT.Tpc = PVT.fase_gas_pressao_temperatura_pseudocritica__Ppr__Tpr__Ppc__Tpc__()[3]
    PVT.Ppr = PVT.fase_gas_pressao_temperatura_pseudocritica__Ppr__Tpr__Ppc__Tpc__()[0]
    PVT.Tpr = PVT.fase_gas_pressao_temperatura_pseudocritica__Ppr__Tpr__Ppc__Tpc__()[1]
    PVT.Z = PVT.fator_z_correlacao_papay()
    PVT.T = T
    PVT.uod = PVT.fase_oleo_viscosidade_do_oleo_morto_beggs_e_robinson_1975__uo__()
    if P <= Pb:
        PVT.P = P
        PVT.uob = PVT.fase_oleo_viscosidade_do_oleo_saturado_beggs_e_robinson_1975_P_menorIgual_Pb__uob__()
        PVT.uo = PVT.uob
        PVT.T = converte_T_para_R(T, 'f')
        PVT.rho_g = PVT.fase_gas_massa_especifica_gas__rho_g__()
        PVT.ug = PVT.fase_gas_viscosidade_do_gas_lee__ug__()
    else:
        PVT.P = Pb
        PVT.uob = PVT.fase_oleo_viscosidade_do_oleo_saturado_beggs_e_robinson_1975_P_menorIgual_Pb__uob__()
        PVT.P = P
        PVT.Pb = Pb
        PVT.uo = PVT.fase_oleo_viscosidade_do_oleo_sub_saturado_beal_standing_1981_P_maiorIgual_Pb__uo__()
    PVT.T = T
    PVT.P = P
    PVT.Bg = PVT.fase_gas_fator_volume_formacao_de_gas__Bg__()
    if PVT.P < PVT.Pb:
        PVT.P = P
        PVT.Pb = Pb
        PVT.Co = PVT.fase_oleo_compressiblidade_isotermica_oleo_P_menor_Pb__Co__()
    else:
        PVT.T = converte_T_para_F(T, 'f')
        PVT.P = Pb
        PVT.Rsb = PVT.fase_oleo_razao_de_solubilidade_standing_1947_P_menorIgual_Pb__Rs__()  # razão de solubilidade
        PVT.Co = PVT.fase_oleo_ompressibilidade_isotermica_oleo_petrosky_e_farshad_1993_P_maiorIgual_Pb__Co__()
    PVT.Cg = PVT.fase_gas_compressibilidade_isotermica_do_gas__Cg__()
    Linha.append(P)
    listaTabela = [PVT.Pb, PVT.Rs, PVT.Bo, PVT.Co, PVT.uo, PVT.Rho_oleo, PVT.Z, PVT.rho_g, PVT.Bg, PVT.Cg, PVT.ug]
    P_P.append(P)
    Rs_Rs.append(PVT.Rs)
    Bo_Bo.append(PVT.Bo)
    Co_Co.append(PVT.Co)
    uo_uo.append(PVT.uo)
    Z_Z.append(PVT.Z)
    rho_g_rho_g.append(PVT.rho_g)
    Bg_Bg.append(PVT.Bg)
    Cg_Cg.append(PVT.Cg)
    ug_ug.append(PVT.ug)
    Rho_oleo_Rho_oleo.append(PVT.Rho_oleo)
    Matrizonha.append(listaTabela)

coluna = 'Pb[Psi] Rs[SCF/STB] Bo[bbl/STB] Co[1/Psi] uo[cP] rho_oleo[lb/ft³] Z rho_gas[lb/ft³] Bg[m³/m³std] Cg[1/Pa] ug[cP]'.split()
linha = Linha
DADOS = np.random.randint(0, 1, len(linha)*len(coluna)).reshape(len(linha), len(coluna))
tabela = pd.DataFrame(data=DADOS, index=linha, columns=coluna)

for l in range(len(tabela)):
    for c in range(0, 11):
        tabela.iat[l, c] = Matrizonha[l][c]
print(tabela)

plt.plot(P_P, Rs_Rs, 'b--o')
plt.title('Rs')
plt.xlabel('Psi')
plt.ylabel('SCF/STB')
plt.show()

plt.plot(P_P, Bo_Bo, 'b--o')
plt.title('Bo')
plt.xlabel('Psi')
plt.ylabel('bbl/STB')
plt.show()

plt.plot(P_P, Co_Co, 'b--o')
plt.title('Co')
plt.xlabel('Psi')
plt.ylabel('1/Psi')
plt.show()

plt.plot(P_P, uo_uo, 'b--o')
plt.title('uo')
plt.xlabel('Psi')
plt.ylabel('cP')
plt.show()

plt.plot(P_P, Rho_oleo_Rho_oleo, 'b--o')
plt.title('rho_oleo')
plt.xlabel('Psi')
plt.ylabel('lb/ft³')
plt.show()

plt.plot(P_P, rho_g_rho_g, 'b--o')
plt.title('rho_gas')
plt.xlabel('Psi')
plt.ylabel('lb/ft³')
plt.show()

plt.plot(P_P, Bg_Bg, 'b--o')
plt.title('Bg')
plt.xlabel('Psi')
plt.ylabel('m³/m³ std')
plt.show()

plt.plot(P_P, Cg_Cg, 'b--o')
plt.title('Cg')
plt.xlabel('Psi')
plt.ylabel('1/Pa')
plt.show()

plt.plot(P_P, ug_ug, 'b--o')
plt.title('ug')
plt.xlabel('Psi')
plt.ylabel('cP')
plt.show()

tabela.to_excel('Tabela_PVT_BlackOil.xlsx', sheet_name='Página-1', index=True)

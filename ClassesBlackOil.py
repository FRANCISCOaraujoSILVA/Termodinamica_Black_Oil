"""
Quando usar a modelagem Black-Oil
    ➢As condições em que se deseja avaliar o fluido está longe do ponto crítico;
    ➢ A Temperatura e/ou a Composição do fluido podem ser consideradas constantes;
    ➢ O fluido apresenta uma única fase;
     Se 𝐵𝑂𝑆 < 2 (na pressão de saturação - lembrar do gráfico): Abordagem Black-Oil;

Tabela PVT Black-Oil para uma temperatura de reservatório constante de modo a contemplar as seguintes propriedades PVT.

Para fase óleo:
    - Pressão de saturação;
    - razão de solubilidade gás-óleo;
    - fator volume-formação do óleo; Bo
    - compressibilidade do óleo; Co
    - viscosidade do óleo.
Para a fase gás:
    - Fator Z;
    - Massa específica do gás;
    - fator volume-formação do gás; bg
    - compressibilidade do gás; cg
    - viscosidade do gás. Ug


> A Tabela PVT deve ser calculada a partir das correlações implementadas numericamente de modo a desenvolver uma
ferramenta de cálculo de fácil manuseio.
> Elaborar três casos distintos e avaliar as propriedades PVT em função da pressão de modo que apresente tanto a região
de óleo sub-saturado como de óleo saturado.
> Cada grupo deve elaborar uma apresentação de forma didática expondo a formulação das correlações, a metodologia de
resolução e desenvolvimento da ferramenta, os resultados e a conclusão. Cada grupo deve entregar a ferramenta de cálculo
desenvolvida e a apresentação elaborada
"""

import math as M
import numpy as np


class FatorZ:
    def __init__(self, Ppr=0, Tpr=0, zc=0, x0=0):
        """
        :param Ppr: Pressão pseudoreduzida, adimensional
        :param Tpr: Temperatura pseudoreduzida, adimensional
        :param zc: z crítico (metano), adimensional
        :param x0: Valor que zera função objetivo
        """
        self.Ppr = Ppr
        self.Tpr = Tpr
        self.zc = zc
        self.x0 = x0

    def fator_z_correlacao_de_brill_e_beggs(self):
        """
        :param Ppr: Pressão pseudoreduzida, adimensional
        :param Tpr: Temperatura pseudoreduzida, adimensional
        :return: Fator de compressibilidade do gás
        """
        Ppr = self.Ppr
        Tpr = self.Tpr

        A = 1.39 * (Tpr - 0.92) ** 0.5 - 0.36 * Tpr - 0.101
        B = (0.62 - 0.23 * Tpr) * Ppr + ((0.066 / (Tpr - 0.86)) - 0.037) * Ppr ** 2 + (
                   0.32 / 10 ** (9 * (Tpr - 1))) * Ppr ** 6
        C = 0.132 - 0.32 * M.log10(Tpr)
        D = 10 ** (0.3106 - 0.49 * Tpr + 0.1824 * Tpr ** 2)
        Z = A + ((1 - A) / M.exp(B)) + C * Ppr ** D
        return Z

    def fator_z_correlacao_papay(self):  # Essa correlação é simples mas tem suas limitações
        """
        :param Ppr: Pressão pseudoreduzida, adimensional
        :param Tpr: Temperatura pseudoreduzida, adimensional
        :return: Fator de compressibilidade do gás
        """
        Ppr = self.Ppr
        Tpr = self.Tpr

        Z = 1 - ((3.53 * Ppr) / (10 ** (0.9813 * Tpr))) + ((0.274 * Ppr ** 2) / (10 ** (0.8157 * Tpr)))
        return Z

    def fator_z_correlacao_de_hall_yarborough(self):
        """
        :param Ppr: Pressão pseudoreduzida, adimensional
        :param Tpr: Temperatura pseudoreduzida, adimensional
        :return: Fator de compressibilidade do gás
        """
        Ppr = self.Ppr
        Tpr = self.Tpr

        t = 1 / Tpr
        X1 = 0.06125 * t * M.exp(-1.2 * (1 - t)**2)
        X2 = 14.76 * t - 9.76 * t**2 + 4.58 * t**3
        X3 = 90.7 * t - 242.2 * t**2 + 42.4 * t**3
        X4 = 2.18 + (2.82 * t)

        Pert = 10 ** -6
        Parad = 10 ** -11
        maxit = 5000
        iter = 0
        x0 = (X1 * Ppr)/FatorZ.fator_z_correlacao_de_brill_e_beggs(self)
        # Aqui, podemos usar a Correlação de Papay também!

        """
        O melhor chute: basta fazer o passo contrário do "Z = (X1 * Ppr) / x0" e isolar x0.
        O melhor valor de z pode ser encontrado pela correlação de Brill e Beggs (se ajustou melhor).
        Com isso, com esse valor de x0 podemos entrar no laço e encontrar o melhor valor de Y, para enfim, encontrar
        o valor de Z.
        """

        F = lambda Y: - X1 * Ppr + ((Y + Y ** 2 + Y ** 3 - Y ** 4) / (1 - Y) ** 3) - X2 * Y ** 2 + X3 * Y ** X4
        while True:
            iter += 1
            xold = x0
            x0 = xold - ((Pert * xold * F(xold)) / (F(xold + Pert * xold) - F(x0)))
            Erro = ((xold - x0) / xold) * 100
            if Erro <= Parad or iter >= maxit:
                break
        Z = (X1 * Ppr) / x0
        return f'FATOR Z PELA CORRELAÇÃO DE HALL-YARBOROUGH>> {Z}'

    def fator_z_correlacao_dranchukabukassem(self):
        """
        :param Ppr: Pressão pseudoreduzida, adimensional
        :param Tpr: Temperatura pseudoreduzida, adimensional
        :param zc: z crítico (metano), adimensional
        :param x0: Valor que zera função objetivo
        :return: Fator de compressibilidade do gás
        """
        Ppr = self.Ppr
        Tpr = self.Tpr
        zc = self.zc
        x0 = self.x0

        A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11 = 0.3265, -1.0700, -0.5339, 0.01569, -0.05165, 0.5475, -0.7361, \
                                                       0.1844, 0.1056, 0.6134, 0.7210
        # Função Objetivo
        F = lambda z: 1 + (A1 + A2 / Tpr + A3 / Tpr ** 3 + A4 / Tpr ** 4 + A5 / Tpr ** 5) * ((zc * Ppr) / (z * Tpr)) + (
                    A6 + A7 / Tpr + A8 / Tpr ** 2) * \
                    ((zc * Ppr) / (z * Tpr)) ** 2 - A9 * (A7 / Tpr + A8 / Tpr ** 2) * ((zc * Ppr) / (z * Tpr)) ** 5 + \
                    A10 * (1 + A11 * ((zc * Ppr) / (z * Tpr)) ** 2) * (
                                (((zc * Ppr) / (z * Tpr)) ** 2) * M.exp(-A11 * ((zc * Ppr) / (z * Tpr)) ** 2)) / (
                                Tpr ** 3) - z

        Pert = 10 ** -6
        Parad = 10 ** -11
        maxit = 50000
        iter = 0

        while True:
            iter += 1
            xold = x0
            x0 = xold - ((Pert * xold * F(xold)) / (F(xold + Pert * xold) - F(x0)))
            Erro = ((xold - x0) / xold) * 100
            if Erro <= Parad or iter >= maxit:
                break
        z = x0
        return f'FATOR Z PELA CORRELAÇÃO DE DRANCHUK & ABU-KASSEM>> {z}'


class BlackOil(FatorZ):
    def __init__(self, Rho_oleo=0, P=0, dgn=0, API=0, do=0, Rs=0, Rsb=0, dg=0, T=0, Tsep=0, Psep=0, Pb=0, Bo=0, Bg=0, Bob=0, Co=0
                 , uob=0, uod=0, Correl_Bo=0, Correl_Rs=0, rho_o_sc=0, rho_g_sc=0, Rho_ob=0,  n=0.172, Z=0,
                 Tpr=0, Tpc=0, Ppr=0, Ppc=0, Mar=28.96, Mg=0, dgas=0, dar=1.225, Psc=14.7, Tsc=60, Yn2=0, Yco2=0,
                 Yh2s=0, R=10.73, rho_g=0, ug=0, Cg=0):


        self.Rho_oleo = Rho_oleo
        self.P = P
        self.dgn = dgn
        self.API = API
        self.do = do
        self.Rs = Rs
        self.Rsb = Rsb
        self.dg = dg
        self.T = T
        self.Tsep = Tsep
        self.Psep = Psep
        self.Pb = Pb
        self.Bo = Bo
        self.Bg = Bg
        self.Bob = Bob
        self.Co = Co
        self.uob = uob
        self.uod = uod
        self.Correl_Bo = Correl_Bo  # Escolher um fator volume formação de óleo para derivar em relação a P
        self.Correl_Rs = Correl_Rs  # Escolher um método de razão de solubilidade do óleo para derivar em relação a P
        self.rho_o_sc = rho_o_sc
        self.rho_g_sc = rho_g_sc
        self.Rho_ob = Rho_ob
        self.n = n
        self.Z = Z
        self.Tpr = Tpr
        self.Tpc = Tpc
        self.Ppr = Ppr
        self.Ppc = Ppc
        self.Mar = Mar
        self.Mg = Mg
        self.dgas = dgas
        self.dar = dar
        self.Psc = Psc
        self.Tsc = Tsc
        self.Yn2 = Yn2
        self.Yco2 = Yco2
        self.Yh2s = Yh2s
        self.R = R
        self.rho_g = rho_g
        self.ug = ug
        self.Cg = Cg

    def fase_oleo_densidade_relativa_do_oleo_com_API__do__(self):
        """
        :param API: Grau API, adimensional
        :return: Densidade relativa do óleo, adimensional
        """
        API = self.API

        do = 141.5 / (API + 131.5)
        return do

    def fase_oleo_grau_API_com_do__API__(self):
        """
        :param do: Densidade relatia do óleo, adimensional
        :return:  Grau API, adimensional
        """
        do = self.do
        API = (141.5/do) - 131.5
        return API

    def pressao_de_bolha_Standing_1947__Pb__(self):
        """
        :param API: Grau API, adimensional
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param dg: Densidade relativa do óleo, adimensional
        :param T: Temperatura, °F
        :return: Pressõa de bolha, psia
        """
        API = self.API
        Rs = self.Rs
        dg = self.dg
        T = self.T

        a = 0.00091 * T - 0.0125 * API
        Pb = 18.2 * (((Rs / dg) ** 0.83) * 10 ** a - 1.4)
        return Pb

    def pressao_de_bolha_Vasquez_e_Beggs_1980__Pb__(self):
        """
        Nota: dgn, gravidade específica do gás normalizada
        :param API: Grau API, adimensional
        :param Tsep: Temperatura no separador, °F
        :param Psep: Pressõa no separador, psia
        :param dg: Densidade relativa do gás, adimensional
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param T: Temperatura, °R
        :return: Pressõa do ponto de bolha, psia
        """
        API = self.API
        Tsep = self.Tsep
        Psep = self.Psep
        dg = self.dg
        Rs = self.Rs
        T = self.T

        if API <= 30:
            C1 = 27.624
            C2 = 0.914328
            C3 = 11.172
        else:
            C1 = 56.18
            C2 = 0.84246
            C3 = 10.393
        dgn = dg * (1 + 5.912 * 10 ** -5 * API * Tsep * np.log10(Psep / 114.7))
        a = -C3 * API / T  # T em °R
        Pb = ((C1 * Rs / dgn) * 10 ** a) ** C2
        return Pb

    def pressao_de_bolha_Glaso_1980__Pb__(self):
        """
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param dg: Densidade relativa do gás, adimensional
        :param n: Para Black_Oil (n=0.172), para óleo_volátil (n=0.1302)
        :param API: Grau API, adimensional
        :param T: Temperatura, °F
        :return: Pressõa do ponto de bolha, psia
        """

        Rs = self.Rs
        dg = self.dg
        API = self.API
        T = self.T
        n = 0.172

        A = ((Rs / dg) ** 0.816) * (T ** n) / API ** 0.989
        Pb = 10 ** (1.7669 + 1.7447 * np.log10(A) - 0.30208 * (np.log10(A)) ** 2)  # Pb em psia
        return Pb

    def pressao_de_bolha_petrosky_e_farshad_1993__Pb__(self):
        """
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param API: Grau API, adimensional
        :param dg: Densidade relativa do gás, adimensional
        :param T: Temperatura, °F
        :return: Pressõa do ponto de bolha, psia
        """
        Rs = self.Rs
        API = self.API
        dg = self.dg
        T = self.T

        a = 7.916 * 10 ** -4 * API ** 1.5410 - 4.561 * 10 ** -5 * T ** 1.3911
        Pb = ((112.727 * Rs ** 0.577421) / (dg ** 0.8439 * 10 ** a)) - 1391.051  # Pb em psia
        return Pb

    def fase_oleo_razao_de_solubilidade_standing_1947_P_menorIgual_Pb__Rs__(self):
        """
        se P > Pb:
        Rs = Rsb  # Rs calculado em P = Pb  (acho que Rsb é a razão de solubilidade na pressõa de bolha, SCF/STB)
        Essa parte eu não entendi.
        Uma dica, para esse caso, é construir o bloco if no main, para avaliar se P<=Pb


        :param P: Pressão, psia
        :param dg: Densidade relativa do gás, adimensional
        :param API: Grau API, adimensional
        :param T: Temperatura, °F
        :return: Razão de solubilidade Gás-Óleo, SCF/STB
        """
        P = self.P
        dg = self.dg
        API = self.API
        T = self.T

        Rs = dg * (((P / 18.2) + 1.4) * 10 ** (0.0125 * API - 0.00091 * T)) ** (1 / 0.83)
        return Rs

    def fase_oleo_razao_de_solubilidade_vasquez_e_beggs_1980_P_menorIgual_Pb__Rs__(self):
        """
        :param API: Grau API, adimensional
        :param P: Pressão, psia
        :param dgn: Gravidade específica do gás normalizada
        :param T: Temperatura, °R
        :return: Razão de solubilidade Gás-Óleo, SCF/STB
        """
        API = self.API
        P = self.P
        dgn = self.dgn
        T = self.T

        if API <= 30:
            C1 = 0.0362
            C2 = 1.0937
            C3 = 25.7240
        else:
            C1 = 0.0178
            C2 = 1.1870
            C3 = 23.931
        Rs = (C1 * dgn * P ** C2) * np.exp(C3 * (API / T))
        return Rs

    def fase_oleo_razao_de_solubilidade_glaso_1980_P_menorIgual_Pb__Rs__(self):
        """
        :param P: Pressão, psia
        :param API: Grau API, adimensional
        :param dg: Densidade relativa do gás, adimensional
        :param T: Temperatura, °F
        :return: Razão de solubilidade Gás-Óleo, SCF/STB
        """
        P = self.P
        API = self.API
        dg = self.dg
        T = self.T

        a = 2.8869 - (14.1811 - 3.3093 * np.log10(P)) ** 0.5
        Rs = dg * (((API ** 0.989) / (T ** 0.1722)) * 10 ** a) ** 1.2255
        return Rs

    def fase_oleo_razao_de_solubilidade_Petrosky_1993_P_menorIgual_Pb__Rs__(self):
        """
        :param P: Pressão, psia
        :param API: Grau API, adimensional
        :param dg: Densidade relativa do gás, adimensional
        :param T: Temperatura, °F
        :return: Razão de solubilidade Gás-Óleo, SCF/STB
        """
        P = self.P
        API = self.API
        dg = self.dg
        T = self.T

        a = (7.916 * 10 ** -4) * (API ** 1.541) - (4.561 * 10 ** -5) * T ** 1.3911
        Rs = (((P / 112.727) + 12.34) * dg ** 0.8439 * 10 ** a) ** 1.73184
        return Rs

    def fase_oleo__compressibilidade_isotermica_oleo_standing_1974_P_maiorIgual_Pb__Co__(self):
        """
        :param P: Pressão, psia
        :param Pb: Pressõa de bolha, psia
        :param Rho_ob: Massa específica do óleo na pressõa de bolha, lb/ft³
        :return: Compressibilidade isotérmica do óleo, 1/psia
        """
        P = self.P
        Pb = self.Pb
        Rho_ob = self.Rho_ob

        Co = 10 ** -6 * np.exp((Rho_ob + 0.004347 * (P - Pb) - 79.1) / (0.0007141 * (P - Pb) - 12.938))  # Co em 1/psia
        return Co

    def fase_oleo_compressibilidade_isotermica_oleo_vasques_e_beggs_1980_P_maiorIgual_Pb__Co__(self):
        """
        :param P: Pressão, psia
        :param Rs: Razão de solubilidade na pressão de bolha, SCF/STB
        :param API: Grau API, adimensional
        :param dgn: Gravidade específica do gás normalizada
        :param T: Temperatura, °F
        :return: Compressibilidade isotérmica do óleo, 1/psia
        """
        P = self.P
        Rs = self.Rs
        API = self.API
        dgn = self.dgn
        T = self.T
        Co = (-1433 + 5 * Rs + 17.2 * T - 1180 * dgn + 12.61 * API) / (10 ** 5 * P)  # Co em 1/psia
        return Co

    def fase_oleo_ompressibilidade_isotermica_oleo_petrosky_e_farshad_1993_P_maiorIgual_Pb__Co__(self):
        """
        :param P: Pressão, psia
        :param Rsb: Razão de solubilidade na pressão de bolha, SCF/STB
        :param API: Grau API, adimensional
        :param dg: Densidade relativa do gás, adimensional
        :param T: Temperatura, °F
        :return: Compressibilidade isotérmica do óleo, 1/psia
        """
        P = self.P
        Rs = self.Rs
        API = self.API
        dg = self.dg
        T = self.T

        Co = (1.705 * 10 ** -7) * (Rs ** 0.69357) * (dg ** 0.1885) * (API ** 0.3272) * (T ** 0.6729) * (P ** -0.5906)
        return Co

    def fase_oleo_compressiblidade_isotermica_oleo_P_menor_Pb__Co__(self):
        """
        Nota: Se dBo/dP < Bg * dRs/dP -> Co > 0.
        Nota: Devo conferir essas novas unidades
        :param Correl_Bo: Escolher um fator volume formação de óleo para derivar em relação a P
        :param Correl_Rs: Escolher uma expressõa de razão de solubilidade para derivar em relação a P
        :param Bo: Fator volume formação de óleo, não falou unidade
        :param Bg: Fator volume formação de gás, não falou unidade
        :return: Compressibilidade isotérmica do óleo, não falou unidade
        """
        # Correl_Bo = self.Correl_Bo
        # Correl_Rs = self.Correl_Rs
        Bo = self.Bo
        Bg = self.Bg
        Bob = self.Bob
        Co = self.Co
        P = self.P
        Pb = self. Pb
        dg = self.dg
        API = self.API
        T = self.T

        # dBo_dP = Correl_Bo
        # dRs_dP = Correl_Rs
        dBo_dP = -Bob * Co * np.exp(-Co * (P - Pb))  # Essa derivada veio do caso quando P>Pb, da função acima
        dRs_dP = dg / 0.83 * (((P/18.2 + 1.4) * 10**(0.0125 * API - 0.00091 * T))**(1/0.83 - 1)) *\
                 10**(0.0125 * API - 0.00091 * T) * (1/18.2)
        # derivada de Rs pela correlação de standing, quando P<=Pb

        Co = (-1 / Bo) * dBo_dP + (Bg / Bo) * dRs_dP
        return Co

    def fase_oleo_fator_volume_formacao_de_oleo_P_maior_Pb__Bo__(self):
        """
        Nota: pg 84, o que seria Bob, seria o fator volume formação de óleo na pressõa de bolha?
        :param P: Pressão, não falou unidade
        :param Pb: Pressão de bolha, não falou unidade
        :param Bob: Fator volume formação calculado na pressão de bolha
        :param Co: Compressibilidade isotérmica do óleo, não falou unidade
        :return:Fator volume formação do óleo, não falou unidade
        """
        P = self.P
        Pb = self.Pb
        Bob = self.Bob
        Co = self.Co

        Bo = Bob * np.exp(-Co * (P - Pb))
        return Bo

    def fase_oleo_fator_volume_formacao_de_oleo_standing_1947_P_menorIgual_Pb__Bo__(self):
        """
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param dg: Densidade relativa do gás, adimensional
        :param do: Densidade relativa do óleo, adimensional
        :param T: Temperatura, °F
        :return:Fator volume formação de óleo, bbl/STB
        """
        Rs = self.Rs
        dg = self.dg
        do = self.do
        T = self.T

        Bo = 0.9759 + 0.00012 * (Rs * (dg / do) ** 0.5 + 1.25 * T) ** 1.2  # Bo em bbl/STB
        return Bo

    def fase_oleo_ator_volume_formacao_de_oleo_vasquez_e_beggs_1980_P_menor_Pb__Bo__(self):
        """
        :param API: Grau API, adimensional
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param dgn: Gravidade específica do gás normalizada
        :param T: Temperatura, °R
        :return: Fator volume formação de óleo, bbl/STB
        """
        API = self.API
        Rs = self.Rs
        dgn = self.dgn
        T = self.T

        if API <= 30:
            C1 = 4.677 * 10 ** -4
            C2 = 1.751 * 10 ** -5
            C3 = -1.811 * 10 ** -8
        else:
            C1 = 4.670 * 10 ** -4
            C2 = 1.100 * 10 ** -5
            C3 = 1.337 * 10 ** -9
        Bo = 1 + C1 * Rs + (T - 520) * (API / dgn) * (C2 + C3 * Rs)
        return Bo

    def fase_oleo_fator_volume_formacao_de_oleo_glaso_1980_P_menor_Pb__Bo__Bob__(self):
        """
        Nota: Usando 45 amostras. Bob é um parâmetro de correção
        :param Rs: Razão de solubilidade Gás-Óleo na pressão de bolha, SCF/STB
        :param do: Densidade relativa do óleo, adimensional
        :param dg: Densidade relativa do gás, adimensional
        :param T: Temperatura, °F
        :return: Fator volume formação de óleo, bbl/STB
        """
        Rs = self.Rs
        do = self.do
        dg = self.dg
        T = self.T

        Bob = Rs * (dg / do) ** 0.526 + 0.986 * T  # Bob é apenas um parâmetro e correlação
        a = -6.58511 + 2.91329 * np.log10(Bob) - 0.27683 * (np.log10(Bob)) ** 2
        Bo = 1 + 10 ** a
        return Bo, Bob

    def fase_oleo_fator_volume_formacao_de_oleo_na_pressao_de_bolha__Bob(self):
        """
        Nota: Usando 45 amostras. Bob é um parâmetro de correção
        :param Rs: Razão de solubilidade na pressõa de bolhaGás-Óleo, SCF/STB
        :param do: Densidade relativa do óleo, adimensional
        :param dg: Densidade relativa do gás, adimensional
        :param T: Temperatura, °F
        :return: Fator volume formação de óleo na pressão de bolha, bbl/STB
        """
        Rs = self.Rs
        do = self.do
        dg = self.dg
        T = self.T

        Bob = Rs * (dg / do) ** 0.526 + 0.986 * T  # Bob é apenas um parâmetro e correlação
        return Bob

    def fase_oleo_fator_volume_formacao_de_oleo_petrosky_e_farshad_1993_P_menor_Pb__Bo__(self):
        """
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param dg: Densidade relativa do gás, adimensional
        :param do: Densidade relativa do óleo, adimensional
        :param T: Temperatura, °F
        :return: Fator volume formação de óleo, bbl/STB
        """
        Rs = self.Rs
        dg = self.dg
        do = self.do
        T = self.T

        Bo = 1.0113 + 7.2046 * 10 ** -5 * (
                    Rs ** 0.3738 * (dg ** 0.2914) / (do ** 0.6265) + 0.24626 * T ** 0.5371) ** 3.0936
        return Bo

    def fase_oleo_massa_especifica_oleo__Rho_oleo__(self):
        """
        Nota: na página 96 não tinha alguma unidades, e na 97 nõa tenho certeza.
        [rho_o_sc=0.1, rho_g_sc=0.1], O que seria: rho_o_sc e rho_g_sc?
        :param P: Pressão, não tinha unidade
        :param Pb: Pressão de bolha, não tinha unidade
        :param Co: Compressibilidade isotérmica do óleo, não tinha unidade
        :param Rho_ob: Será será a densidade do gás na pressõa de bolha???????????
        :param Bo: Fator volume formação de óleo, bbl/STB
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param do: Densidade relativa do óleo, adiminensional
        :param dg: Densidade relativa do gás, adimensional
        :return: Massa específica do óleo, lb/ft³
        """
        rho_o_sc = self.rho_o_sc
        rho_g_sc = self.rho_g_sc
        P = self.P
        Pb = self.Pb
        Co = self.Co
        Rho_ob = self.Rho_ob
        Bo = self.Bo
        Rs = self.Rs
        do = self.do
        dg = self.dg

        if P > Pb:
            Rho_oleo = Rho_ob * np.exp(Co * (P - Pb))
        elif P <= Pb:
            Rho_oleo = (62.4 * do + 0.0136 * Rs * dg) / Bo  # lb/ft³
        else:  # Para P<Pb
            Rho_oleo = (rho_o_sc + Rs * rho_g_sc) / Bo
        return Rho_oleo

    def fase_oleo_viscosidade_do_oleo_morto_beal_standing_1981__uo__(self):
        """
        :param API: Grau API, adimensional
        :param T: Temperatura, °R
        :return: Viscosidade do óleo morto (uod), cP
        """
        API = self.API
        T = self.T

        A = 10 ** (0.43 + (8.33 / API))
        uod = 0.32 + ((1.8 * 10 ** 7) / (API ** 4.53)) * (360 / (T - 260)) ** A
        return uod

    def fase_oleo_viscosidade_do_oleo_morto_beggs_e_robinson_1975__uo__(self):
        """
        :param API: Grau API, adimensional
        :param T: Temperatura, °F
        :return: Viscosidade do óleo morto (uod), cP
        """
        API = self.API
        T = self.T

        A = 10 ** (3.0324 - (0.02023 * API))
        uod = (10 ** (A * (T ** -1.163))) - 1
        return uod

    def fase_oleo_viscosidade_do_oleo_morto_bergman_2004__uod__(self):
        """
        Nota: O último logaritmo da fórmula é o logaritmo natural em python.
        :param API: Grau API, adimensional
        :param T: Temperatura, °F
        :return: Viscosidade do óleo morto (uod), cP
        """
        API = self.API
        T = self.T

        X = np.exp(22.33 - 0.194 * API + 0.00033 * API ** 2 - (3.2 - 0.0185 * API) * np.log(T + 310))
        uod = np.exp(X) - 1
        return uod

    def fase_oleo_viscosidade_do_oleo_saturado_standing_1981_P_menorIgual_Pb__uob__(self):
        """
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param uod: Viscosidade do óleo morto, cP
        :return: Viscosiddde do óleo saturado, cP
        """
        Rs = self.Rs
        uod = self.uod
        a = 10 ** ((-7.4 * 10 ** -4 * Rs) + (2.2 * 10 ** -7 * Rs ** 2))
        b = (0.68 / (10 ** ((8.62 * 10 ** -5) * Rs))) + (0.25 / (10 ** ((1.1 * 10 ** -3) * Rs))) + (
                    0.062 / (10 ** ((3.74 * 10 ** -3) * Rs)))
        uob = a * uod ** b
        return uob

    def fase_oleo_viscosidade_do_oleo_saturado_beggs_e_robinson_1975_P_menorIgual_Pb__uob__(self):
        """
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param uod: Viscosidade do óleo morto, cP
        :return: Viscosiddde do óleo saturado, cP
        """
        Rs = self.Rs
        uod = self.uod

        a = 10.715 * (Rs + 100) ** -0.515
        b = 5.44 * (Rs + 150) ** -0.338
        uob = a * uod ** b
        return uob

    def fase_oleo_viscosidade_do_oleo_saturado_bergman_1975_P_menorIgual_Pb__uob__(self):
        """
        :param Rs: Razão de solubilidade Gás-Óleo, SCF/STB
        :param uod: iscosidade do óleo morto, cP
        :return: Viscosiddde do óleo saturado, cP
        """
        Rs = self.Rs
        uod = self.uod

        a = np.exp(4.768 - 0.8359 * np.log(Rs + 300))
        b = 0.555 + (133.5/(Rs + 300))
        uob = a * uod ** b
        return uob

    def fase_oleo_viscosidade_do_oleo_sub_saturado_beal_standing_1981_P_maiorIgual_Pb__uo__(self):
        """
        :param P: Pressão, psi
        :param Pb: Pressõa de bolha, psi
        :param uob: Viscosidade do óleo saturado na pressão de bolha, cP
        :return: Viscosidade do óleo sub-saturado, cP
        """
        P = self.P
        Pb = self.Pb
        uob = self.uob

        uo = uob + (0.001 * (P - Pb)) * (0.024 * uob ** 1.6 + 0.038 * uob ** 0.56)
        return uo

    def fase_oleo_viscosidade_do_oleo_sub_saturado_beggs_e_robinson_1975_P_maiorIgual_Pb__uo__(self):
        """
        :param P: Pressão, psi
        :param Pb: Pressõa de bolha, psi
        :param uob: Viscosidade do óleo saturado na pressão de bolha, cP
        :return: Viscosidade do óleo sub-saturado, cP
        """
        P = self.P
        Pb = self.Pb
        uob = self.uob

        m = 2.6 * (P ** 1.187) * np.exp(-11.513 - 8.98 * (10 ** -5) * P)
        uo = uob * (P / Pb) ** m
        return uo

    def fase_oleo_viscosidade_do_oleo_sub_saturado_bergman_2004_P_maiorIgual_Pb__uo__(self):
        """
        :param P: Pressão, psi
        :param Pb: Pressõa de bolha, psi
        :param uob: Viscosidade do óleo saturado na pressão de bolha, cP
        :return: Viscosidade do óleo sub-saturado, cP
        """
        P = self.P
        Pb = self.Pb
        uob = self.uob

        alpha = 6.5698 * (10 ** -7) * (np.log(uob)) ** 2 - 1.48211 * (10 ** -5) * (np.log(uob)) + 2.27877 * 10 ** -4
        betha = 2.24623 * (10 ** -2) * (np.log(uob)) + 0.873204
        uo = uob * np.exp(alpha * (P - Pb) ** betha)
        return uo

    def fase_gas_massa_especifica_gas__rho_g__(self):
        """
        :param P: Pressão, Pa
        :param T: Temperatura, K
        :param Mg: Massa molecular do gás, g/mol
        :param R: Constante universal dos gases, 8.314 (Pa⋅m³)/(mol⋅K)54
        :param Z: Fator de compressibilidade isotérmico do gás, adimensional
        :return: rho_g
        """
        P = self.P
        Mg = self.Mg
        Z = self.Z
        R = self.R
        T = self.T

        rho_g = (P*Mg)/(Z*R*T)
        return rho_g

    def fase_gas_massa_do_gas__Mg__(self):
        """
        Nota: Mar = 28.96 lb/lbmol
        Mar e Mg precisam ter a mesma unidade
        :param dg: Densidade relativa do gás, adimensional
        :param Mar: Mar do ar
        :return: Massa do gás, adimensional
        """
        dg = self.dg
        Mar = self.Mar

        Mg = Mar * dg
        return Mg

    def fase_gas_densidade_relativa_do_gas__dg__(self):
        """
        Nota: Pressão na condição padrão: Pa, Temperatura na condição padrão: K
        Nota: Em 101,325 KPa e 15°C, dar = 1.225 kg/m³
        :param dgas: Massa específica de um gás na condição padrão (kg/m³)
        :param dar: Massa específica do ar na condição padrão (kg/m³)
        :return: dg (densidade relativa do gás), kg/m³,
        Ppc (pressõa pseudocrítica) - adimensional,
        Tpc (temperatura pseudocrítica) - adimensional
        """
        dgas = self.dgas
        dar = self.dar

        dg = dgas / dar
        return dg

    def fase_gas_pressao_temperatura_pseudocritica__Ppr__Tpr__Ppc__Tpc__(self):
        dg = self.dg
        P = self.P
        T = self.T
        if dg < 0.75:  # gás seco
            Ppc = (677 + 15 * dg - 37.5 * dg ** 2)
            Tpc = (168 + 325 * dg - 12.5 * dg ** 2)
        else:  # gás úmido
            Ppc = (706 - 51.7 * dg - 11.1 * dg ** 2)
            Tpc = (187 + 330 * dg - 71.5 * dg ** 2)

        Ppr = P/Ppc
        Tpr = T/Tpc

        return Ppr, Tpr, Ppc, Tpc

    def fase_gas_fator_volume_formacao_de_gas__Bg__(self):
        """
        :param P: Pressão, mesma unidade de Psc
        :param Z: Fator de compressibilidade (unidimensional)
        :param T: Temperatura, mesma unidade de Tsc
        :param Psc: 14.7 psia (pressõa na condição de superfície)
        :param Tsc: 60 °F (temperatura na condição de superfície)
        :return: Bg, Fator volume formação do gás
        """
        P = self.P
        Z = self.Z
        T = self.T
        Psc = self.Psc
        Tsc = self.Tsc

        Bg = (Psc / Tsc) * ((Z * T) / P)
        return Bg

    def fase_gas_compressibilidade_isotermica_do_gas__Cg__(self):
        """
        Nota: dZ_dPpr foi obtido a partir da derivada manual da Correlação de Papay
        :param Z: Fator de compressiblidade do gás, adimensional
        :param Tpr: Temperatura pseudoreduzida
        :param Ppr: Pressõa pseudoreduzida
        :param Ppc: Pressão pseudocrítica
        :return: Compressibilidde isotérmica do gás
        """
        Z = self.Z
        Tpr = self.Tpr
        Ppr = self.Ppr
        Ppc = self.Ppc

        dZ_dPpr = - (3.53/(10**(0.9813*Tpr))) + ((2*0.274 * Ppr) / (10**(0.8157 * Tpr)))
        Cpr = 1/Ppr - (1/Z) * dZ_dPpr
        Cg = Cpr/Ppc
        return Cg

    def fase_gas_viscosidade_do_gas_dempsey_1965__ug__(self):
        """
        Nota: UG, viscosidade do gás não corrigida, cP. Na página 48, será que a correlação está incompleta?
        :param T: Temperatura, °F
        :param dg: Densidade relativa do gás, adimensional
        :param Yn2: Fração molar do componente
        :param Yco2: Fração molar do componente
        :param Yh2s: Fração molar do componente
        :return: Viscosidade do gás na condição de pressão atmosférica, cP
        """
        dg = self.dg
        T = self.T
        Yn2 = self.Yn2
        Yco2 = self.Yco2
        Yh2s = self.Yh2s

        UG = (1.709 * 10 ** -5 - 2.062 * 10 ** -6 * dg) * T + 8.188 * 10 ** -3 - 6.15 * 10 ** -3 * np.log10(dg)
        delta_ug_n2 = Yn2 * (8.43 * 10 ** -3 * np.log10(dg) + 9.59 * 10 ** -3)
        delta_ug_co2 = Yco2 * (9.08 * 10 ** -3 * np.log10(dg) + 6.24 * 10 ** -3)
        delta_ug_h2s = Yh2s * (8.49 * 10 ** -3 * np.log10(dg) + 3.73 * 10 ** -3)
        ug = UG + delta_ug_n2 + delta_ug_co2 + delta_ug_h2s
        return ug

    def fase_gas_viscosidade_do_gas_lee__ug__(self):
        """
        :param Mg: Peso molecular do gás, lbm/(lb mol)
        :param densidade_gas: densidade do gás, lb/ft³
        :param T: Tempeeratura, °R
        :return: Viscosidade do gás, cP
        """
        Mg = self.Mg
        rho_g = self.rho_g
        T = self.T

        xv = 3.448 + (986.4 / T) + 0.01009 * Mg
        yv = 2.4 - 0.2 * xv
        kv = ((9.379 + 0.0160 * Mg) * T ** 1.5) / (209.2 + 19.26 * Mg + T)
        ug = (10 ** -4) * kv * np.exp(xv * (rho_g / 62.4) ** yv)
        return ug

    def fase_gas_viscosidade_do_gas_sutton_2007__ug__(self):
        """
        :param Mg: Peso molecular, g/(g mol)
        :param dgas: Densidade do gás, g/cm³
        :param Tpr: Temperatura pseudoreduzida, adimensional
        :param Tpc: Temperatura pseudocrítica, adimensional
        :param Ppc: Pressão pseudocrítica, adimensional
        :param T: Temperatura, °R
        :return: Viscosidade do gás, cP
        """
        Mg = self.Mg
        dgas = self.dgas   # Acho que aqui é rho_g
        Tpr = self.Tpr
        Tpc = self.Tpc
        Ppc = self.Ppc
        T = self.T

        X = 3.47 + (1588 / T) + 0.0009 * Mg
        Y = 1.66378 - 0.04679 * X
        A = X * dgas ** Y
        K = (0.807 * Tpr ** 0.618 - 0.357 * np.exp(-0.449 * Tpr) + 0.34 * np.exp(-4.058 * Tpr) + 0.018) / \
            (0.9490 * (Tpc / (Mg ** 3 * Ppc ** 4)) ** (1 / 6))
        ug = 10 ** -4 * K * np.exp(A)
        return ug


"""------------------------------------------------------------------------------------------------------------------"""
"Def's para converter"


def converte_T_para_R(*args):  # inserir a temperatura com a unidade dentro de uma string. Converte para Rankine
    T = list(args)

    if T[1].upper() == "C":
        T = float(T[0]) * 9 / 5 + 491.67
    elif T[1].upper() == "R":
        T = float(T[0])
    elif T[1].upper() == "F":
        T = (float(T[0]) + 459.67)
    elif T[1].upper() == "K":
        T = float(T[0]) * 1.8
    return T


def converte_T_para_F(*args):  # inserir a temperatura com a unidade dentro de uma string. Converter para Fahrenheit
    T = list(args)

    if T[1].upper() == "C":
        T = float(T[0]) * 1.8 + 32
    elif T[1].upper() == "R":
        T = float(T[0]) - 459.67
    elif T[1].upper() == "F":
        T = float(T[0])
    elif T[1].upper() == "K":
        T = (float(T[0]) - 273.15) * 1.8 - 32
    return T


def converte_T_para_K(*args):  # inserir a temperatura com a unidade dentro de uma string. Converter para Kelvin
    T = list(args)

    if T[1].upper() == "C":
        T = float(T[0]) + 273.15
    elif T[1].upper() == "R":
        T = float(T[0]) * 0.556
    elif T[1].upper() == "F":
        T = (float(T[0]) - 32)/1.8 + 273.15
    elif T[1].upper() == "K":
        T = float(T[0])
    return T


def converte_P_para_Psi(*args):  # inserir a temperatura com a unidade dentro de uma string. Converter para Psi
    P = list(args)

    if P[1] == "BAR":
        P = float(P[0]) * 14.514
    elif P[1].upper() == "PA":
        P = float(P[0]) / 6895
    elif P[1].upper() == "ATM":
        P = float(P[0]) * 14.696
    elif P[1].upper() == "TORR":
        P = float(P[0]) / 51.715
    elif P[1].upper() == "MMHG":
        P = float(P[0]) / 51.715
    elif P[1].upper() == "KGF/CM2":
        P = float(P[0]) * 14.223
    elif P[1].upper() == "KGF/IN2":
        P = float(P[0]) * 2.205
    elif P[1].upper() == "PSI":
        P = float(P[0])
    return P


def converte_P_para_Pa(*args):
    P = list(args)

    if P[1] == "BAR":
        P = float(P[0]) * 100000
    elif P[1].upper() == "PA":
        P = float(P[0])
    elif P[1].upper() == "ATM":
        P = float(P[0]) * 101325
    elif P[1].upper() == "TORR":
        P = float(P[0]) * 133.322
    elif P[1].upper() == "MMHG":
        P = float(P[0]) * 133.32
    elif P[1].upper() == "KGF/CM2":
        P = float(P[0]) * 98066.5
    elif P[1].upper() == "KGF/IN2":
        P = float(P[0]) * 6894759.09
    elif P[1].upper() == "PSI":
        P = float(P[0]) * 6894.76
    return P

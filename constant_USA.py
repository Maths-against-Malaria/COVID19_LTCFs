import numpy as np

# Number of Erlang stages NE, NP, NI, NL
NH = (16, 16, 16, 16)

# Durations at stages DE,DP,DI,DL. Dtest is the waiting time before tests results (days)
DH = (3.7, 1, 5, 5)
Dtest = [2]  # (1/2, 1, 2, 3, 4)

# Testing rate per Erlang stage (assumes that we have the same number of Erlang stages
DXi = [5]  # (1, 2, 5, 7, 14)  # Time spent between 2 consecutive tests (days)
test_rate = [1/i for i in DXi]

# drifting rates (epsilon, phi, gamma, delta)
Rate = [i / j for i, j in zip(NH, DH)]
epsi, phi, gamma, delta = Rate

# daily incoming infections in the general and staff population
LExt = 50

# Population (General, Staff, Risk group)
N = 331e6
NPop = (329.177e6, 423e3, 1.4e6)

N_init_inf = 100

# Initial number of individuals in the different compartments
y0 = [0 for i in np.arange(5 * (NH[0] + NH[1] + NH[2] + NH[3] + 2) + 3)]
y0[0], y0[1], y0[2] = (329.177e6 - N_init_inf, 423e3, 1.4e6)  # Susceptible individuals in SGe, SSt, SRi
y0[NH[0] + NH[1] + 3] = N_init_inf  # Initial fully infectious individuals in the general population

# Sensitivity of the test sH
sH = np.array([[0,    0,     0,    0],              # No test
               [0,    0.10,  0.65, 0.35],           # Poor test
               [0.03, 0.30,  0.75, 0.50],           # Intermediate test
               [0.15, 0.60,  0.80, 0.60],           # Good test
               [0.25, 0.75,  0.90, 0.65],           # Very good test
               [0.30, 0.80,  0.95, 0.85]])          # Excellent test

# Contagiousness at the infectious stages P, I, L
cH = (0.5, 1, 0.5)

# Probability  fsick of showing symptoms in the Ge, and St sub-populations
fsick = 0.50

# Probability fiso of being isolated in the Ge, and St sub-populations
fiso = 0.48

# Probability fDead of dying from the COVID-19 in the Ge, and St sub-populations
fdead = 0.04

# Probability of showing symptoms and die in the general and staff populations
fside = fsick * fdead

# Probability of being sick and isolated
fsiso = fsick * fiso

# Probability  fsick of showing symptoms in the Ri sub-populations
fsick_Ri = 0.70

# Probability fDead of dying from the Covid-19 in the Ri St sub-populations
fdead_Ri = 0.07

# Probability of showing symptoms and die in the risk group
fside_ri = fsick_Ri * fdead_Ri

# Maximum carrying capacity Qmax of quarantine ward (ICU capacity) per 10000 inhabitants
Qmax = 30 * N / 10000

# Simulation time span
sim_time = 800

# sustainability period of quarantine measures (till the end of the simulation)
tiso = (20, sim_time)

# Sustainability period for the general distancing measures (first lockdown)
tdista = 35  # March    15th, 2020
tdistb = 75  # May      29th
tdistc = 185  # August   02nd
tdistd = 200  # October  01st
tdiste = 250  # November 05th
tdistf = 260  # November 24th
tdistg = 309  # December 01st
tdisth = 315  # December 10th
tdisti = 325  # December 20th
tdistj = 350  # January  08th, 2021
tdistk = 500  # April    15th

# Additional test control for Ge per disease stage (P, I, L)
fpos_mat = np.array([[0, 0, 0],
                     [0, 0, 0]])

# Effectiveness of home isolation in the general population, and the LTCF employees
phome = 0.75

# Parameters of the basic reproduction number R0 (_R0, a, tR0max, seasonality)
ParamRO = [[3.2, 0.35, 354, "yes"]]

# Contact rates at stages P, I, L
denBeta = sum([i * j for i, j in zip(cH, DH[1:])])

# Mapping the ranges of the population list
# Stage E
EGeL = 3
EGeR = EGeL + NH[0]

EStnL = NH[0] + NH[1] + NH[2] + NH[3] + 3
EStnR = EStnL + NH[0]

EStsL = 2 * EStnL - 3
EStsR = EStsL + NH[0]

EStpL = 3 * EStnL - 6
EStpR = EStpL + NH[0]

ERiL = 4 * EStnL - 9
ERiR = ERiL + NH[0]

# Stage P
PGeL = EGeR
PGeR = PGeL + NH[1]

PStnL = EStnR
PStnR = PStnL + NH[1]

PStsL = EStsR
PStsR = PStsL + NH[1]

PStpL = EStpR
PStpR = PStpL + NH[1]

PRiL = ERiR
PRiR = PRiL + NH[1]

# Stage I
IGeL = PGeR
IGeR = IGeL + NH[2]

IStnL = PStnR
IStnR = IStnL + NH[2]

IStsL = PStsR
IStsR = IStsL + NH[2]

IStpL = PStpR
IStpR = IStpL + NH[2]

IRiL = PRiR
IRiR = IRiL + NH[2]

# Stage L
LGeL = IGeR
LGeR = LGeL + NH[3]

LStnL = IStnR
LStnR = LStnL + NH[3]

LStsL = IStsR
LStsR = LStsL + NH[3]

LStpL = IStpR
LStpR = LStpL + NH[3]

LRiL = IRiR
LRiR = LRiL + NH[3]

# Recovered
RecGe = LRiR  # 5 * NH[0] + 5 * NH[1] + 5 * NH[2] + 5 * NH[3] + 3
RecStn = RecGe + 1
RecSts = RecStn + 1
RecStp = RecSts + 1
RecRi = RecStp + 1

# Dead
DeGe = RecRi + 1
DeStn = DeGe + 1
DeSts = DeStn + 1
DeStp = DeSts + 1
DeRi = DeStp + 1

# Population intervals
HLeft = [EGeL, EStnL, EStsL, EStpL, ERiL, PGeL, PStnL, PStsL, PStpL, PRiL, IGeL, IStnL, IStsL, IStpL, IRiL, LGeL, LStnL,
         LStsL, LStpL, LRiL]
HRight = [EGeR, EStnR, EStsR, EStpR, ERiR, PGeR, PStnR, PStsR, PStpR, PRiR, IGeR, IStnR, IStsR, IStpR, IRiR, LGeR,
          LStnR, LStsR, LStpR, LRiR]

# Total contacts in the General population, Staff population, and risk group per day.
nPop = (60, 60, 30)

xm = 0.9992
ym = 0.0006
um = ym*nPop[0]*NPop[0]/(nPop[1]*NPop[1])
vm = 0.20

pm = (1 - xm - ym)*nPop[0]*NPop[0]/(nPop[2]*NPop[2])
qm = (1 - um - vm)*nPop[1]*NPop[1]/(nPop[2]*NPop[2])

pcontreduc_WT_mat = np.array([[0, 0.60, 0.25, 0.50, 0.40, 0.60, 0.45, 0.55, 0.70, 0.55, 0.70, 0],      # pGe
                              [0, 0.60, 0.25, 0.50, 0.40, 0.60, 0.45, 0.55, 0.70, 0.55, 0.70, 0],      # pGeSt
                              [0, 0.60, 0.45, 0.50, 0.40, 0.55, 0.50, 0.60, 0.65, 0.60, 0.65, 0],      # pRiGe
                              [0, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55],   # pStSt
                              [0, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40],   # pStRi
                              [0, 0.20, 0.15, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10]])  # pRiRi

pcontreduc_NT_mat = np.array([[0, 0.60, 0.25, 0.50, 0.40, 0.60, 0.45, 0.55, 0.70, 0.55, 0.70, 0],      # pGe
                              [0, 0.60, 0.25, 0.50, 0.40, 0.60, 0.45, 0.55, 0.70, 0.55, 0.70, 0],      # pGeSt
                              [0, 0.60, 0.45, 0.50, 0.40, 0.55, 0.50, 0.60, 0.65, 0.60, 0.65, 0],      # pRiGe
                              [0, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55, 0.55],   # pStSt
                              [0, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 0.40],   # pStRi
                              [0, 0.20, 0.15, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.20, 0.10]])  # pRiRi
# Unadjusted mixing matrix with contact reduction (with testing)
Xfinal_init_NT = []
Xfinal_init_WT = []

xinit = np.array(
    [[xm * nPop[0] * N ** 2 / NPop[0],   ym * nPop[0] * N ** 2 / NPop[1],  (1 - xm - ym) * nPop[0] * N ** 2 / NPop[2]],
     [um * nPop[1] * N ** 2 / NPop[0],   vm * nPop[1] * N ** 2 / NPop[1],  (1 - um - vm) * nPop[1] * N ** 2 / NPop[2]],
     [pm * nPop[2] * N ** 2 / NPop[0],   qm * nPop[2] * N ** 2 / NPop[1],  (1 - pm - qm) * nPop[2] * N ** 2 / NPop[2]]])

for i in range(12):
    pcontreduc_loop = list(pcontreduc_WT_mat[:, i])

    M_mat = np.array([[1-pcontreduc_loop[0], 1-pcontreduc_loop[1],  1-pcontreduc_loop[2]],
                      [1-pcontreduc_loop[1], 1-pcontreduc_loop[3],  1-pcontreduc_loop[4]],
                      [1-pcontreduc_loop[2], 1-pcontreduc_loop[4],  1-pcontreduc_loop[5]]])
    Xfinal_init_WT.append(np.multiply(xinit, M_mat))

# Unadjusted contact matrix with no test
for i in range(12):
    pcontreduc_loop = list(pcontreduc_NT_mat[:, i])

    M_mat = np.array([[1-pcontreduc_loop[0], 1-pcontreduc_loop[1],  1-pcontreduc_loop[2]],
                      [1-pcontreduc_loop[1], 1-pcontreduc_loop[3],  1-pcontreduc_loop[4]],
                      [1-pcontreduc_loop[2], 1-pcontreduc_loop[4],  1-pcontreduc_loop[5]]])
    Xfinal_init_NT.append(np.multiply(xinit, M_mat))

# Number of infected compartments
n_inf_comp = 3 * sum(NH)

F = np.zeros((n_inf_comp, n_inf_comp))
V = np.zeros((n_inf_comp, n_inf_comp))
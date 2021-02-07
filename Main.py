# Import the libraries and required functions
import pandas as pd
from scipy.integrate import solve_ivp
from Model import *

Func = FuncClass()
Mod = ModelFunc()

colonne = ["population", "susceptible", "infected", "recovered", "dead", "test", "testing_rate", "waiting_time",
           "seasonal_fluctuation", "simulation_step"]
df = pd.DataFrame(columns=colonne)

aa = len(ParamRO)
bb = len(sH)
cc = len(test_rate)
dd = len(Dtest)

for i in range(aa):  # with and without seasonality
    SeFl = ParamRO[i]
    BRep = SeFl[0:3]
    season = SeFl[3]

    for j in range(bb):  # Test Quality
        TestQual = sH[j, ]
        fh = []
        if j == 0:      # No test
            FPos = fpos_mat[0]
            XF_Init = Xfinal_init_NT
        else:           # Testing control
            FPos = fpos_mat[1]
            XF_Init = Xfinal_init_WT

        for k in range(cc):  # Testing rate
            xi = np.repeat(test_rate[k], 4)
            fh = [TestQ * TestR for TestQ, TestR in zip(TestQual, xi)]
            f_e, f_p, f_i, f_l = fh
            for cnt in range(dd):  # Waiting time

                df2 = pd.DataFrame(columns=colonne)  # Empty data frame

                alpha = 1 / Dtest[cnt]

                # Solving the system
                soln = solve_ivp(lambda t, y: Mod.f(t, y, BRep, f_e, f_p, f_i, f_l, alpha, FPos, XF_Init), [0, sim_time]
                                 , y0, method="RK45", dense_output=True)

                # Building the data frame

                # Column: Susceptibles
                SGe = soln.y[0]
                SSt = soln.y[1]
                SRi = soln.y[2]

                l = len(SGe)

                SStpos = [0 for i in range(l)]

                Susc = SGe.tolist() + SSt.tolist() + SStpos + SRi.tolist() 

                # Column: Population
                Ge = ["general" for i in range(l)]
                St = ["staff" for i in range(l)]
                Stpos = ["staff_removed" for i in range(l)]
                Ri = ["risk" for i in range(l)]

                popul = Ge + St + Stpos + Ri 

                # Column: Infected

                EGe = [0 for i in range(l)]
                PGe = [0 for i in range(l)]
                IGe = [0 for i in range(l)]
                LGe = [0 for i in range(l)]

                for i in range(3, EGeR):
                    EGe = [m + n for m, n in zip(EGe, soln.y[i])]
                for i in range(PGeL, PGeR):
                    PGe = [m + n for m, n in zip(PGe, soln.y[i])]
                for i in range(IGeL, IGeR):
                    IGe = [m + n for m, n in zip(IGe, soln.y[i])]
                for i in range(LGeL, LGeR):
                    LGe = [m + n for m, n in zip(LGe, soln.y[i])]

                out11 = [m + n for m, n in zip(EGe, PGe)]
                out12 = [m + n for m, n in zip(IGe, LGe)]
                out1 = [m + n for m, n in zip(out11, out12)]

                # Staff population

                ESt = [0 for i in range(l)]
                PSt = [0 for i in range(l)]
                ISt = [0 for i in range(l)]
                LSt = [0 for i in range(l)]

                for i in range(EStnL, EStnR):
                    ESt = [m + n for m, n in zip(ESt, soln.y[i])]
                for i in range(PStnL, PStnR):
                    PSt = [m + n for m, n in zip(PSt, soln.y[i])]
                for i in range(IStnL, IStnR):
                    ISt = [m + n for m, n in zip(ISt, soln.y[i])]
                for i in range(LStnL, LStnR):
                    LSt = [m + n for m, n in zip(LSt, soln.y[i])]

                for i in range(EStsL, EStsR):
                    ESt = [m + n for m, n in zip(ESt, soln.y[i])]
                for i in range(PStsL, PStsR):
                    PSt = [m + n for m, n in zip(PSt, soln.y[i])]
                for i in range(IStsL, IStsR):
                    ISt = [m + n for m, n in zip(ISt, soln.y[i])]
                for i in range(LStsL, LStsR):
                    LSt = [m + n for m, n in zip(LSt, soln.y[i])]

                for i in range(EStpL, EStpR):
                    ESt = [m + n for m, n in zip(ESt, soln.y[i])]
                for i in range(PStpL, PStpR):
                    PSt = [m + n for m, n in zip(PSt, soln.y[i])]
                for i in range(IStpL, IStpR):
                    ISt = [m + n for m, n in zip(ISt, soln.y[i])]
                for i in range(LStpL, LStpR):
                    LSt = [m + n for m, n in zip(LSt, soln.y[i])]

                out21 = [m + n for m, n in zip(ESt, PSt)]
                out22 = [m + n for m, n in zip(ISt, LSt)]
                out2 = [m + n for m, n in zip(out21, out22)]

                # Removed staff (St,+)
                EStpos = [0 for i in range(l)]
                PStpos = [0 for i in range(l)]
                IStpos = [0 for i in range(l)]
                LStpos = [0 for i in range(l)]

                for i in range(EStpL, EStpR):
                    EStpos = [m + n for m, n in zip(EStpos, soln.y[i])]
                for i in range(PStpL, PStpR):
                    PStpos = [m + n for m, n in zip(PStpos, soln.y[i])]
                for i in range(IStpL, IStpR):
                    IStpos = [m + n for m, n in zip(IStpos, soln.y[i])]
                for i in range(LStpL, LStpR):
                    LStpos = [m + n for m, n in zip(LStpos, soln.y[i])]

                out31 = [m + n for m, n in zip(EStpos, PStpos)]
                out32 = [m + n for m, n in zip(IStpos, LStpos)]
                out3 = [m + n for m, n in zip(out31, out32)]

                # Risk group population
                ERi = [0 for i in range(l)]
                PRi = [0 for i in range(l)]
                IRi = [0 for i in range(l)]
                LRi = [0 for i in range(l)]

                for i in range(ERiL, ERiR):
                    ERi = [m + n for m, n in zip(ERi, soln.y[i])]
                for i in range(PRiL, PRiR):
                    PRi = [m + n for m, n in zip(PRi, soln.y[i])]
                for i in range(IRiL, IRiR):
                    IRi = [m + n for m, n in zip(IRi, soln.y[i])]
                for i in range(LRiL, LRiR):
                    LRi = [m + n for m, n in zip(LRi, soln.y[i])]

                out41 = [m + n for m, n in zip(ERi, PRi)]
                out42 = [m + n for m, n in zip(IRi, LRi)]
                out4 = [m + n for m, n in zip(out41, out42)]

                Infec = out1 + out2 + out3 + out4

                # Recovered
                RGe = soln.y[RecGe]

                RSt = [0 for i in range(l)]
                for i in range(RecStn, RecStp + 1):
                    RSt = [m + n for m, n in zip(RSt, soln.y[i])]

                RStpos = soln.y[RecStp]

                RRi = soln.y[RecRi]

                Recov = RGe.tolist() + RSt + RStpos.tolist() + RRi.tolist()

                # Death toll
                DGe = soln.y[DeGe]

                DSt = [0 for i in range(l)]
                for i in range(DeStn, DeStp + 1):
                    DSt = [m + n for m, n in zip(DSt, soln.y[i])]

                DStpos = soln.y[DeStp]

                DRi = soln.y[DeRi]

                Dead = DGe.tolist() + DSt + DStpos.tolist() + DRi.tolist()

                # Column: Test quality
                if j == 0:
                    TQ = ["NT" for ii in range(4 * l)]
                elif j == 1:
                    TQ = ["Poor" for ii in range(4 * l)]
                elif j == 2:
                    TQ = ["Intermediate" for ii in range(4 * l)]
                elif j == 3:
                    TQ = ["Good" for ii in range(4 * l)]
                elif j == 4:
                    TQ = ["Very_good" for ii in range(4 * l)]
                else:
                    TQ = ["Excellent" for ii in range(4 * l)]

                # Column: Test rate
                if k == 0:
                    TR = ["everyday" for ii in range(4 * l)]
                elif k == 1:
                    TR = ["every_two_days" for ii in range(4 * l)]
                elif k == 2:
                    TR = ["every_five_days" for ii in range(4 * l)]
                elif k == 3:
                    TR = ["once_per_week" for ii in range(4 * l)]
                else:
                    TR = ["once_every_2_weeks" for ii in range(4 * l)]

                #  Column: Waiting Time
                if cnt == 0:
                    WT = ["half" for ii in range(4 * l)]
                elif cnt == 1:
                    WT = ["one" for ii in range(4 * l)]
                elif cnt == 2:
                    WT = ["two" for ii in range(4 * l)]
                elif cnt == 3:
                    WT = ["three" for ii in range(4 * l)]
                else:
                    WT = ["four" for ii in range(4 * l)]

                # Column: Seasonal fluctuation

                SF = [season for ii in range(4 * l)]

                # Column: simulation time
                sim_step = 4 * soln.t.tolist()

                # saving the simulated data
                df2 = pd.DataFrame({"population": popul, "susceptible": Susc, "infected": Infec, "recovered": Recov,
                                    "dead": Dead, "test": TQ, "testing_rate": TR, "waiting_time": WT,
                                    "seasonal_fluctuation": SF,
                                    "simulation_step": sim_step})

                # saving the final dataset
                df = df.append(df2, ignore_index=True)

# Saving the data from the simulation
df.to_csv("covid_19_risk_group_simulation_LTCF_QT_NCR.csv")

print("Simulation finished successfully. Please check your data")

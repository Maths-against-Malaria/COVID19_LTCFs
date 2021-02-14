from Functions import *

func = FuncClass()

# Beta parameter at time t = 0
beta = func.betah((cH, func.r0(0, ParamRO[0][0:3]) / denBeta))

# The adjusted mixing matrix
coef_NGM = func.contmat(beta/func.r0(0, ParamRO[0][0:3]), Xfinal_init_NT[0])

class ModelFunc:
    """This class contains the function that builds the differential equations of the model"""

    def f(self, t, y, brep, fe, fp, fi, fl, alpha, fpos, Xfinal_init):
        """This function computes the differential equations"""

        # Initialize the outputs
        out = [0 for i in range(5 * (sum(NH) + 2) + 3)]

        # Total population per sub-population and disease stage
        h_sum = [func.hsum(y, i, j) for i, j in zip(HLeft, HRight)]

        # Total number of individuals that can be isolated in quarantine wards
        aq = fsiso * (func.hsum(h_sum, 10, 13) + func.hsum(h_sum, 15, 18)) + sum([h_sum[i] for i in np.arange(3, 20, 5)]
                                                                                 )
        # Time information
        tps = (t, tiso[0], tiso[1])

        # Adjusted mixing matrix
        if tdista <= t < tdistb:
            mixmat = np.multiply(Xfinal_init[1], coef_NGM)
        elif tdistb <= t < tdistc:
            mixmat = np.multiply(Xfinal_init[2], coef_NGM)
        elif tdistc <= t < tdistd:
            mixmat = np.multiply(Xfinal_init[3], coef_NGM)
        elif tdistd <= t < tdiste:
            mixmat = np.multiply(Xfinal_init[4], coef_NGM)
        elif tdiste <= t < tdistf:
            mixmat = np.multiply(Xfinal_init[5], coef_NGM)
        elif t >= tdistf:  # tdistf <= t < tdistg:  #
            mixmat = np.multiply(Xfinal_init[6], coef_NGM)
        # Uncomment the lines below to implement the extra contact reductions of the USA model
        #elif tdistg <= t < tdisth:
        #    mixmat = np.multiply(Xfinal_init[7], coef_NGM)
        #elif tdisth <= t < tdisti:
        #    mixmat = np.multiply(Xfinal_init[8], coef_NGM)
        #elif tdisti <= t < tdistj:
        #   mixmat = np.multiply(Xfinal_init[9], coef_NGM)
        #elif tdistj <= t < tdistk:
        #    mixmat = np.multiply(Xfinal_init[10], coef_NGM)
        #elif t >= tdistk:
        #    mixmat = np.multiply(Xfinal_init[11], coef_NGM)
        else:
            mixmat = np.multiply(Xfinal_init[0], coef_NGM)

        # Effective infective population at stage H
        peff_stp = func.heff(phome,
                             (h_sum[8], func.hiso(tps, (aq, Qmax, h_sum[8])), func.hhome(tps, (aq, Qmax, h_sum[8]))))

        ieff_ge = func.heff(phome,
                            (h_sum[10], func.hiso(tps, (aq, Qmax, fsiso * h_sum[10])),
                             func.hhome(tps, (aq, Qmax, fsiso * h_sum[10]))))
        ieff_geri = func.heff(1,
                              (h_sum[10], func.hiso(tps, (aq, Qmax, fsiso * h_sum[10])),
                               func.hhome(tps, (aq, Qmax, fsiso * h_sum[10]))))

        ieff_stn = func.heff(phome,
                             (h_sum[11], func.hiso(tps, (aq, Qmax, fsiso * h_sum[11])),
                              func.hhome(tps, (aq, Qmax, fsiso * h_sum[11]))))
        ieff_stnri = func.heff(1,
                               (h_sum[11], func.hiso(tps, (aq, Qmax, fsiso * h_sum[11])),
                                func.hhome(tps, (aq, Qmax, fsiso * h_sum[11]))))

        ieff_sts = func.heff(phome,
                             (h_sum[12], func.hiso(tps, (aq, Qmax, fsiso * h_sum[12])),
                              func.hhome(tps, (aq, Qmax, fsiso * h_sum[12]))))
        ieff_stsri = func.heff(1,
                               (h_sum[12], func.hiso(tps, (aq, Qmax, fsiso * h_sum[12])),
                                func.hhome(tps, (aq, Qmax, fsiso * h_sum[12]))))

        ieff_stp = func.heff(phome, (
            h_sum[13], func.hiso(tps, (aq, Qmax, h_sum[13])), func.hhome(tps, (aq, Qmax, fsiso * h_sum[13]))))

        ieff_ri = func.heffri((h_sum[14], func.hisori(tps, fsick_Ri * h_sum[14])))

        leff_ge = func.heff(phome,
                            (h_sum[15], func.hiso(tps, (aq, Qmax, fsiso * h_sum[15])),
                             func.hhome(tps, (aq, Qmax, fsiso * h_sum[15]))))
        leff_geri = func.heff(1,
                              (h_sum[15], func.hiso(tps, (aq, Qmax, fsiso * h_sum[15])),
                               func.hhome(tps, (aq, Qmax, fsiso * h_sum[15]))))

        leff_stn = func.heff(phome,
                             (h_sum[16], func.hiso(tps, (aq, Qmax, fsiso * h_sum[11])),
                              func.hhome(tps, (aq, Qmax, fsiso * h_sum[16]))))
        leff_stnri = func.heff(1,
                               (h_sum[16], func.hiso(tps, (aq, Qmax, fsiso * h_sum[11])),
                                func.hhome(tps, (aq, Qmax, fsiso * h_sum[16]))))

        leff_sts = func.heff(phome,
                             (h_sum[17], func.hiso(tps, (aq, Qmax, fsiso * h_sum[17])),
                              func.hhome(tps, (aq, Qmax, fsiso * h_sum[17]))))
        leff_stsri = func.heff(1,
                               (h_sum[17], func.hiso(tps, (aq, Qmax, fsiso * h_sum[17])),
                                func.hhome(tps, (aq, Qmax, fsiso * h_sum[17]))))



        leff_stp = func.heff(phome,
                             (h_sum[18], func.hiso(tps, (aq, Qmax, h_sum[18])),
                              func.hhome(tps, (aq, Qmax, fsiso * h_sum[18]))))

        leff_ri = func.heffri((h_sum[19], func.hisori(tps, fsick_Ri * h_sum[19])))


        # Beta parameter at time t
        beta_h = func.betah((cH, func.r0(t, brep) / denBeta))

        # Force of infection in the sub-populations Ge, St, and Ri
        l_ge = (LExt + (beta_h[0] * sum(
            [i * j for i, j in zip(mixmat[0, ], [h_sum[5], sum([h_sum[6], h_sum[7], peff_stp]), h_sum[9]])])
                        + beta_h[1] * sum(
                    [i * j for i, j in zip(mixmat[0, ], [ieff_ge, sum([ieff_stn, ieff_sts, ieff_stp]), ieff_ri])])
                        + beta_h[2] * sum(
                    [i * j for i, j in zip(mixmat[0, ], [leff_ge, sum([leff_stn, leff_sts, leff_stp]), leff_ri])]))) / N

        l_st = (LExt + (beta_h[0] * sum(
            [i * j for i, j in zip(mixmat[1, ], [h_sum[5], sum([h_sum[6], h_sum[7], peff_stp]), h_sum[9]])])
                        + beta_h[1] * sum(
                    [i * j for i, j in zip(mixmat[1, ], [ieff_ge, sum([ieff_stn, ieff_sts, ieff_stp]), ieff_ri])])
                        + beta_h[2] * sum(
                    [i * j for i, j in zip(mixmat[1, ], [leff_ge, sum([leff_stn, leff_sts, leff_stp]), leff_ri])]))) / N

        l_ri = (mixmat[2, 0] * (beta_h[0] * (1 - func.in_interval_dist(t, fpos[0])) * h_sum[5]
                                + beta_h[1] * (1 - func.in_interval_dist(t, fpos[1])) * ieff_geri
                                + beta_h[2] * (1 - func.in_interval_dist(t, fpos[2])) * leff_geri)
                + beta_h[0] * sum([i * j for i, j in zip(mixmat[2, 1:], [sum([h_sum[6], h_sum[7], peff_stp]), h_sum[9]])
                                   ])
                + beta_h[1] * sum([i * j for i, j in zip(mixmat[2, 1:], [sum([ieff_stnri, ieff_stsri]), ieff_ri])])
                + beta_h[2] * sum([i * j for i, j in zip(mixmat[2, 1:], [sum([leff_stnri, leff_stsri]), leff_ri])])) / N

        ###################################### Susceptible population ###############################################

        # General population
        out[0] = -l_ge * y[0]

        # LTCF Employees
        out[1] = -l_st * y[1]

        # LTCF Resident risk group
        out[2] = -l_ri * y[2]

        ###################################### population dynamic at the early infected stage E #####################

        # General population
        out[3] = l_ge * y[0] - epsi * y[3]
        for i in range(4, EGeR):
            out[i] = epsi * (y[i - 1] - y[i])

            # LTCF Employees
            # St,-
        out[EStnL] = l_st * y[1] - (epsi + fe) * y[EStnL]
        for i in range(EStnL + 1, EStnR):
            out[i] = epsi * (y[i - 1] - y[i]) - fe * y[i]

            # St,*
        out[EStsL] = fe * y[EStnL] - (epsi + alpha)*y[EStsL]
        for i, k in zip(range(EStsL + 1, EStsR), range(EStnL + 1, EStnR)):
            out[i] = fe * y[k] + epsi * (y[i - 1] - y[i]) - alpha * y[i]

            # St,+
        out[EStpL] = alpha * y[EStsL] - epsi * y[EStpL]
        for i, k in zip(range(EStpL + 1, EStpR), range(EStsL + 1, EStsR)):
            out[i] = alpha * y[k] + epsi * (y[i - 1] - y[i])

            # LTCF Resident risk group
        out[ERiL] = l_ri * y[2] - epsi * y[ERiL]
        for i in range(ERiL + 1, ERiR):
            out[i] = epsi * (y[i - 1] - y[i])

        ###################################### Population dynamic at the prodromal stage P #####################

        # General population
        out[PGeL] = epsi * y[EGeR - 1] - phi * y[PGeL]
        for i in range(PGeL + 1, PGeR):
            out[i] = phi * (y[i - 1] - y[i])

            # LTCF Employees
            # St,-
        out[PStnL] = epsi * y[EStnR - 1] - (phi + fp) * y[PStnL]
        for i in range(PStnL + 1, PStnR):
            out[i] = phi * (y[i - 1] - y[i]) - fp * y[i]

            # St,*
        out[PStsL] = epsi * y[EStsR - 1] + fp * y[PStnL] - (phi + alpha) * y[PStsL]
        for i, k in zip(range(PStsL + 1, PStsR), range(PStnL + 1, PStnR)):
            out[i] = fp * y[k] + phi * (y[i - 1] - y[i]) - alpha * y[i]

            # St,+
        out[PStpL] = epsi * y[EStpR - 1] + alpha * y[PStsL] - phi * y[PStpL]
        for i, k in zip(range(PStpL + 1, PStpR), range(PStsL + 1, PStsR)):
            out[i] = alpha * y[k] + phi * (y[i - 1] - y[i])

            # LTCF Resident risk group
        out[PRiL] = epsi * y[ERiR - 1] - phi * y[PRiL]
        for i in range(PRiL + 1, PRiR):
            out[i] = phi * (y[i - 1] - y[i])

        ############################# population dynamic at the fully infectious stage I ########################

        # General population
        out[IGeL] = phi * y[PGeR - 1] - gamma * y[IGeL]
        for i in range(IGeL + 1, IGeR):
            out[i] = gamma * (y[i - 1] - y[i])

            # LTCF Employees
            # St,-
        out[IStnL] = phi * y[PStnR - 1] - (gamma + fi) * y[IStnL]
        for i in range(IStnL + 1, IStnR):
            out[i] = gamma * (y[i - 1] - y[i]) - fi * y[i]

            # St,*
        out[IStsL] = phi * y[PStsR - 1] + fi * y[IStnL] - (gamma + alpha) * y[IStsL]
        for i, k in zip(range(IStsL + 1, IStsR), range(IStnL + 1, IStnR)):
            out[i] = fi * y[k] + gamma * (y[i - 1] - y[i]) - alpha * y[i]

            # St,+
        out[IStpL] = phi * y[PStpR - 1] + alpha * y[IStsL] - gamma * y[IStpL]
        for i, k in zip(range(IStpL + 1, IStpR), range(IStsL + 1, IStsR)):
            out[i] = alpha * y[k] + gamma * (y[i - 1] - y[i])

            # LTCF Resident risk group
        out[IRiL] = phi * y[PRiR - 1] - gamma * y[IRiL]
        for i in range(IRiL + 1, IRiR):
            out[i] = gamma * (y[i - 1] - y[i])

        ############################# population dynamic at the late infectious stage L ########################

        # General population
        out[LGeL] = gamma * y[IGeR - 1] - delta * y[LGeL]
        for i in range(LGeL + 1, LGeR):
            out[i] = delta * (y[i - 1] - y[i])

            # LTCF Employees
            # St,-
        out[LStnL] = gamma * y[IStnR - 1] - (delta + fl) * y[LStnL]
        for i in range(LStnL + 1, LStnR):
            out[i] = delta * (y[i - 1] - y[i]) - fl * y[i]

            # St,*
        out[LStsL] = gamma * y[IStsR - 1] + fl * y[LStnL] - (delta + alpha) * y[LStsL]
        for i, k in zip(range(LStsL + 1, LStsR), range(LStnL + 1, LStnR)):
            out[i] = fl * y[k] + delta * (y[i - 1] - y[i]) - alpha * y[i]

            # St,+
        out[LStpL] = gamma * y[IStpR - 1] + alpha * y[LStsL] - delta * y[LStpL]
        for i, k in zip(range(LStpL + 1, LStpR), range(LStsL + 1, LStsR)):
            out[i] = alpha * y[k] + delta * (y[i - 1] - y[i])

            # LTCF Resident risk group
        out[LRiL] = gamma * y[IRiR - 1] - delta * y[LRiL]
        for i in range(LRiL + 1, LRiR):
            out[i] = delta * (y[i - 1] - y[i])

        ############################# population dynamic at the recovered stage R ########################

        # General population
        out[RecGe] = delta * (1 - fsick * fdead) * y[LGeR - 1]

        # LTCF Employees
        # St,-
        out[RecStn] = delta * (1 - fsick * fdead) * y[LStnR - 1]

        # St,*
        out[RecSts] = delta * (1 - fsick * fdead) * y[LStsR - 1]

        # St,+
        out[RecStp] = delta * (1 - fsick * fdead) * y[LStpR - 1]

        # LTCF Resident risk group
        out[RecRi] = delta * (1 - fsick_Ri * fdead_Ri) * y[LRiR - 1]

        ############################# Death toll D ########################

        # General population
        out[DeGe] = delta * fsick * fdead * y[LGeR - 1]

        # LTCF Employees
        # St,-
        out[DeStn] = delta * fsick * fdead * y[LStnR - 1]

        # St,*
        out[DeSts] = delta * fsick * fdead * y[LStsR - 1]

        # St,+
        out[DeStp] = delta * fsick * fdead * y[LStpR - 1]

        # LTCF Resident risk group
        out[DeRi] = delta * fsick_Ri * fdead_Ri * y[LRiR - 1]

        return out

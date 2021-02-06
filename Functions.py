#from constant_Germany import *
from constant_USA import *    # Uncomment this line and comment the line above if you are simulating USA

class FuncClass:
    """This class contains all the functions we need to define the model"""

    def pcont(self, t, inp_var):
        """Returns the parameter for general contact reduction in Staff and general
           population. The parameter is active during the period from tdist1 to tdist2.
           Parameters:
               t = time list containing;
                t0 = time t;
                t1 = day when the distancing measures start;
                t2 = day when the distancing measures end;
               inp_var = tuple of general contact reduction in the sub-populations.
        """
        if tdista <= t < tdistb:
            return inp_var[0]
        elif tdistb <= t < tdistc:
            return inp_var[1]
        elif tdistc <= t < tdistd:
            return inp_var[2]
        elif tdistd <= t < tdiste:
            return inp_var[3]
        elif tdiste <= t < tdistf:
            return inp_var[4]
        elif t >= tdistf:
            return inp_var[5]
        else:
            return 0, 0, 0, 0, 0, 0

    def pcontrige(self, t, inp_var):
        """Returns the parameter for general contact reduction between the resident risk group and
           the general population. The parameter is active during the period from tdist1 to tdist2, and from
           tdist2 to tiso2.
           Parameters:
               t = time list containing;
                t0 = time t;
                t1 = day when the distancing measures start;
                t2 = day when the distancing measures end;
                t3 = day when the isolation measures start;
                t4 = day when the isolation measures end;
               inp_var = probability for Ge to be tested positive at  the disease stages P, I, L.
               pcontrigeh = contact reduction between risk group and general population at stage H.
        """
        if tdista <= t < tdistb:
            return inp_var[0]
        elif tdistb <= t < tdistc:
            return inp_var[1]
        elif tdistc <= t < tdistd:
            return inp_var[2]
        elif tdistd <= t < tdiste:
            return inp_var[3]
        elif tdiste <= t < tdistf:
            return inp_var[4]
        elif t >= tdistf:
            return inp_var[5]
        else:
            return 0

        #  if tdista <= t < tdistb:
         #   return inp_var[0]
        # elif tdistb <= t < tdistc:
          #  return inp_var[1]
       # elif tdistc <= t < tdistd:
         #    return inp_var[2]
        # else:
          #  return inp_var[3]

    def in_interval_dist(self, t, inp_var):
        """Returns the input parameter if time frame is within the period when distancing measure apply.
        The parameter is active during the period from tdist1 to tdist2.
           Parameters:
                t = time list containing;
                t0 = time t;
                t1 = day when the distancing measures start;
                t2 = day when the distancing measures end;
               inp_var = probability for Ge to be tested positive at  the disease stages P, I, L.
        """

        # is_in_interval_dist = t[2] >= t[0] >= t[1]
        #print(t)
        if tdista <= t <= tdistd:
            return inp_var
        else:
            return 0
        # if is_in_interval_dist:
        #    return inp_var
        # else:
        #    return 0

    def r0(self, t, inp_var):
        """This function returns the value of the seasonal basic reproduction number at time t.
        parameters:
            t = time t;
            inp_var = tuple containing:
                R0 = average basic reproductive number;
                a = amplitude of the seasonal effect (0<a<1);
                tR0Max = day in the year when the R0 is at max.
        """
        return inp_var[0] * (1 + inp_var[1] * np.cos(2 * np.pi * (t - inp_var[2]) / 365))

    def betah(self, inp_var):
        """This function returns the value of the rate of infective contact $\beta_H$ at stage H (P, I, L) and time t
        parameters:
            inp_var = tuple containing:
                cH = list of contagiousness level of stage H (P, I, L);
                R0(t)/denBeta.
        """
        out = [i * inp_var[1] for i in inp_var[0]]
        return out

    def hsum(self, pop, inp_var1, inp_var2):
        """This function returns the total number of individuals in one of the three sub-populations (Ge, St, Ri)
        who are at a particular disease stage H (P, I, L).
        parameters:
            pop = list containing the number of individuals for each sub-population in all the compartments of the model;
            inp_var1, inp_var2 = the range HL, HR in which the values are located."""

        return sum([pop[i] for i in range(inp_var1, inp_var2)])

    def hiso(self, t, inp_var):
        """This function returns the total number of individuals from the general population and the staff
        who are isolated in the quarantine ward at time t.
        parameters:
            t = tuple containing:
                t0 = current time;
                t1 = day when isolation measures start;
                t2 = day when isolation measures end.
            inp_var = tuple containing:
                aq = total number of individuals who can be isolated at time t0;
                Qmax = quarantine maximum carrying capacity;
                popIso = total number of individuals that are showing symptoms and are isolated at stage H
                (fsick*fiso*Hsum for Ge, St-, St* and Hsum for St+).
        """

        is_in_interval_iso_with_space = (t[2] >= t[0] >= t[1] and inp_var[0] <= inp_var[1])
        is_in_interval_iso_no_space = (t[2] >= t[0] >= t[1] and inp_var[0] > inp_var[1])

        if is_in_interval_iso_with_space:
            return inp_var[2]
        elif is_in_interval_iso_no_space:
            return inp_var[2] * inp_var[1] / inp_var[0]
        else:
            return 0

    def hisori(self, t, inp_var):
        """This function returns the total number of individuals from the resident risk group who are isolated
        in the quarantine ward at time t.
        parameters:
            t = tuple containing:
                t0 = current time;
                t1 = day when isolation measures start;
                t2 = day when isolation measures end.
            inp_var = total number of individuals that are showing symptoms at stage H (fsick*Hsum).
        """

        is_in_interval_iso = (t[2] >= t[0] >= t[1])

        if is_in_interval_iso:
            return inp_var
        else:
            return 0

    def hhome(self, t, inp_var):
        """This function returns the total number of individuals from the general population and the staff
        who are isolated at home at time t.
        parameters:
            t = tuple containing:
                t0 = current time;
                t1 = day when isolation measures start;
                t2 = day when isolation measures end.
            inp_var = tuple containing:
                aq = total number of individuals who can be isolated at time t0;
                Qmax = quarantine maximum carrying capacity;
                popIso = total number of individuals that are showing symptoms and are isolated at stage H
                (fsick*fiso*Hsum for Ge, St-, St* and Hsum for St+).
        """

        is_in_interval_iso_no_space = (t[2] >= t[0] >= t[1] and inp_var[0] > inp_var[1])

        if is_in_interval_iso_no_space:
            return inp_var[2]*(1 - inp_var[1] / inp_var[0])
        else:
            return 0

    def heff(self, cont_red_home, inp_var):
        """This function returns the total number of individuals from the general and staff population who can effectively
         infect susceptible individuals in the General and staff populations.
         parameter:
            cont_red_home = parameter defining the contact reduction in home isolation;
            inp_var = a tuple containing hsum, hiso, hhome for a given sub-population at a particular stage H (P, I, L).
        """
        return inp_var[0] - inp_var[1] - cont_red_home * inp_var[2]

    def heffri(self, inp_var):
        """This function returns the total number of individuals from the general and staff population who can effectively
         infect susceptible individuals in the risk group.
         parameter:
            inp_var = a tuple containing hsum, hisori for a given sub-population at a particular stage H (P, I, L);
        """
        return inp_var[0] - inp_var[1]

    def contmat(self, inp_var, xinit):
        """This function finds the value of R0 in the system at the initial time. The value obtained  is therefore divided
        by the value of R0 we want in our population. The coefficient that is obtained is used to multiply the mixing matrix
        parameter:
            inp_var =
        """

        # Initial mixing matrix
        xGG = xinit[0][0]
        xGS = xinit[0][1]
        xGR = xinit[0][2]
        xSG = xinit[1][0]
        xSS = xinit[1][1]
        xSR = xinit[1][2]
        xRG = xinit[2][0]
        xRS = xinit[2][1]
        xRR = xinit[2][2]

        betaGe = [i * NPop[0] / N for i in inp_var]
        # E1(Ge)
        for i in np.arange(3 * NH[0], 3 * NH[0] + NH[1]):  # Derivative with respect to Pk(Ge)
            F[0, i] = betaGe[0] * xGG

        for i in np.arange(3 * NH[0] + NH[1], 3 * NH[0] + 2 * NH[1]):  # Derivative with respect to Pk(St,-)
            F[0, i] = betaGe[0] * xGS

        for i in np.arange(3 * NH[0] + 2 * NH[1], 3 * NH[0] + 3 * NH[1]):  # Derivative with respect to Pk(Ri)
            F[0, i] = betaGe[0] * xGR

        for i in np.arange(3 * NH[0] + 3 * NH[1], 3 * NH[0] + 3 * NH[1] + NH[2]):  # Derivative with respect to Ik(Ge)
            F[0, i] = betaGe[1] * xGG

        for i in np.arange(3 * NH[0] + 3 * NH[1] + NH[2],
                           3 * NH[0] + 3 * NH[1] + 2 * NH[2]):  # Derivative with respect to Ik(St,-)
            F[0, i] = betaGe[1] * xGS

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 2 * NH[2],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2]):  # Derivative with respect to Ik(Ri)
            F[0, i] = betaGe[1] * xGR

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3]):  # Derivative with respect to Lk(Ge)
            F[0, i] = betaGe[2] * xGG

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3]):  # Derivative with respect to Lk(St,-)
            F[0, i] = betaGe[2] * xGS

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 3 * NH[3]):  # Derivative with respect to Lk(Ri)
            F[0, i] = betaGe[2] * xGR

        betaSt = [i * NPop[1] / N for i in inp_var]

        # E1(st,-)
        for i in np.arange(3 * NH[0], 3 * NH[0] + NH[1]):  # Derivative with respect to Pk(Ge)
            F[NH[0], i] = betaSt[0] * xSG

        for i in np.arange(3 * NH[0] + NH[1], 3 * NH[0] + 2 * NH[1]):  # Derivative with respect to Pk(St,-)
            F[NH[0], i] = betaSt[0] * xSS

        for i in np.arange(3 * NH[0] + 2 * NH[1], 3 * NH[0] + 3 * NH[1]):  # Derivative with respect to Pk(Ri)
            F[NH[0], i] = betaSt[0] * xSR

        for i in np.arange(3 * NH[0] + 3 * NH[1], 3 * NH[0] + 3 * NH[1] + NH[2]):  # Derivative with respect to Ik(Ge)
            F[NH[0], i] = betaSt[1] * xSG

        for i in np.arange(3 * NH[0] + 3 * NH[1] + NH[2],
                           3 * NH[0] + 3 * NH[1] + 2 * NH[2]):  # Derivative with respect to Ik(St,-)
            F[NH[0], i] = betaSt[1] * xSS

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 2 * NH[2],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2]):  # Derivative with respect to Ik(Ri)
            F[NH[0], i] = betaSt[1] * xSR

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3]):  # Derivative with respect to Lk(Ge)
            F[NH[0], i] = betaSt[2] * xSG

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3]):  # Derivative with respect to Lk(St,-)
            F[NH[0], i] = betaSt[2] * xSS

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 3 * NH[3]):  # Derivative with respect to Lk(Ri)
            F[NH[0], i] = betaSt[2] * xSR

        betaRi = [i * NPop[2] / N for i in inp_var]

        # E1(Ri)
        for i in np.arange(3 * NH[0], 3 * NH[0] + NH[1]):  # Derivative with respect to Pk(Ge)
            F[2 * NH[0], i] = betaRi[0] * xRG

        for i in np.arange(3 * NH[0] + NH[1], 3 * NH[0] + 2 * NH[1]):  # Derivative with respect to Pk(St,-)
            F[2 * NH[0], i] = betaRi[0] * xRS

        for i in np.arange(3 * NH[0] + 2 * NH[1], 3 * NH[0] + 3 * NH[1]):  # Derivative with respect to Pk(Ri)
            F[2 * NH[0], i] = betaRi[0] * xRR

        for i in np.arange(3 * NH[0] + 3 * NH[1], 3 * NH[0] + 3 * NH[1] + NH[2]):  # Derivative with respect to Ik(Ge)
            F[2 * NH[0], i] = betaRi[1] * xRG

        for i in np.arange(3 * NH[0] + 3 * NH[1] + NH[2],
                           3 * NH[0] + 3 * NH[1] + 2 * NH[2]):  # Derivative with respect to Ik(St,-)
            F[2 * NH[0], i] = betaRi[1] * xRS

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 2 * NH[2],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2]):  # Derivative with respect to Ik(Ri)
            F[2 * NH[0], i] = betaRi[1] * xRR

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3]):  # Derivative with respect to Lk(Ge)
            F[2 * NH[0], i] = betaRi[2] * xRG

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3]):  # Derivative with respect to Lk(St,-)
            F[2 * NH[0], i] = betaRi[2] * xRS

        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3],
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 3 * NH[3]):  # Derivative with respect to Lk(Ri)
            F[2 * NH[0], i] = betaRi[2] * xRR

        # Building the drift matrix V of the next generation matrix

        # E stage
        # Ge
        V[0, 0] = Rate[0]  # Derivative with respect to E1(Ge)
        for i in np.arange(1, NH[0]):  # Derivative with respect to Ek(Ge)
            V[i, i - 1] = -Rate[0]
            V[i, i] = Rate[0]

        # St,-
        V[NH[0], NH[0]] = Rate[0]  # Derivative with respect to E1(St,-)
        for i in np.arange(NH[0] + 1, 2 * NH[0]):  # Derivative with respect to Ek(St,-)
            V[i, i - 1] = -Rate[0]
            V[i, i] = Rate[0]

        # Ri
        V[2 * NH[0], 2 * NH[0]] = Rate[0]  # Derivative with respect to E1(Ri)
        for i in np.arange(2 * NH[0] + 1, 3 * NH[0]):  # Derivative with respect to Ek(Ri)
            V[i, i - 1] = -Rate[0]
            V[i, i] = Rate[0]

        # P stage
        # Ge
        V[3 * NH[0], NH[0] - 1] = -Rate[0]  # Derivative with respect to EnE(Ge)
        V[3 * NH[0], 3 * NH[0]] = Rate[1]  # Derivative with respect to P1(Ge)
        for i in np.arange(3 * NH[0] + 1, 3 * NH[0] + NH[1]):  # Derivative with respect to Pk(Ge)
            V[i, i - 1] = -Rate[1]
            V[i, i] = Rate[1]

        # St,-
        V[3 * NH[0] + NH[1], 2 * NH[0] - 1] = -Rate[0]  # Derivative with respect to EnE(St,-)
        V[3 * NH[0] + NH[1], 3 * NH[0] + NH[1]] = Rate[1]  # Derivative with respect to P1(st,-)
        for i in np.arange(3 * NH[0] + NH[1] + 1, 3 * NH[0] + 2 * NH[1]):  # Derivative with respect to Pk(st,-)
            V[i, i - 1] = -Rate[1]
            V[i, i] = Rate[1]

        # Ri
        V[3 * NH[0] + 2 * NH[1], 3 * NH[0] - 1] = -Rate[0]  # Derivative with respect to EnE(Ri)
        V[3 * NH[0] + 2 * NH[1], 3 * NH[0] + 2 * NH[1]] = Rate[1]  # Derivative with respect to P1(Ri)
        for i in np.arange(3 * NH[0] + 2 * NH[1] + 1, 3 * NH[0] + 3 * NH[1]):  # Derivative with respect to Pk(Ri)
            V[i, i - 1] = -Rate[1]
            V[i, i] = Rate[1]

        # I stage
        # Ge
        V[3 * NH[0] + 3 * NH[1], 3 * NH[0] + NH[1] - 1] = -Rate[1]  # Derivative with respect to PnP(Ge)
        V[3 * NH[0] + 3 * NH[1], 3 * NH[0] + 3 * NH[1]] = Rate[2]  # Derivative with respect to I1(Ge)
        for i in np.arange(3 * NH[0] + 3 * NH[1] + 1,
                           3 * NH[0] + 3 * NH[1] + NH[2]):  # Derivative with respect to Ik(Ge)
            V[i, i - 1] = -Rate[2]
            V[i, i] = Rate[2]

        # St,-
        V[3 * NH[0] + 3 * NH[1] + NH[2], 3 * NH[0] + 2 * NH[1] - 1] = -Rate[1]  # Derivative with respect to PnP(St,-)
        V[3 * NH[0] + 3 * NH[1] + NH[2], 3 * NH[0] + 3 * NH[1] + NH[2]] = Rate[2]  # Derivative with respect to I1(St,-)
        for i in np.arange(3 * NH[0] + 3 * NH[1] + NH[2] + 1,
                           3 * NH[0] + 3 * NH[1] + 2 * NH[2]):  # Derivative with respect to Ik(St,-)
            V[i, i - 1] = -Rate[2]
            V[i, i] = Rate[2]

        # Ri
        V[3 * NH[0] + 3 * NH[1] + 2 * NH[2], 3 * NH[0] + 3 * NH[1] - 1] = -Rate[1]  # Derivative with respect to PnP(Ri)
        V[3 * NH[0] + 3 * NH[1] + 2 * NH[2], 3 * NH[0] + 3 * NH[1] + 2 * NH[2]] = Rate[
            2]  # Derivative with respect to I1(Ri)
        for i in np.arange(3 * NH[0] + 3 * NH[1] + 2 * NH[2] + 1,
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2]):  # Derivative with respect to Ik(Ri)
            V[i, i - 1] = -Rate[2]
            V[i, i] = Rate[2]

        # L stage
        # Ge
        V[3 * NH[0] + 3 * NH[1] + 3 * NH[2], 3 * NH[0] + 3 * NH[1] + NH[2] - 1] = -Rate[
            2]  # Derivative with respect to InI(Ge)
        V[3 * NH[0] + 3 * NH[1] + 3 * NH[2], 3 * NH[0] + 3 * NH[1] + 3 * NH[2]] = Rate[
            3]  # Derivative with respect to L1(Ge)
        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 1,
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3]):  # Derivative with respect to Lk(Ge)
            V[i, i - 1] = -Rate[3]
            V[i, i] = Rate[3]

        # St,-
        V[3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3], 3 * NH[0] + 3 * NH[1] + 2 * NH[2] - 1] = -Rate[
            2]  # Derivative with respect to InI(St,-)
        V[3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3], 3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3]] = Rate[
            3]  # Derivative with respect to L1(St,-)
        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + NH[3] + 1,
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3]):  # Derivative with respect to Lk(St,-)
            V[i, i - 1] = -Rate[3]
            V[i, i] = Rate[3]

        # Ri
        V[3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3], 3 * NH[0] + 3 * NH[1] + 3 * NH[2] - 1] = -Rate[
            2]  # Derivative with respect to InI(Ri)
        V[3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3], 3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3]] = Rate[
            3]  # Derivative with respect to L1(Ri)
        for i in np.arange(3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 2 * NH[3] + 1,
                           3 * NH[0] + 3 * NH[1] + 3 * NH[2] + 3 * NH[3]):  # Derivative with respect to Lk(Ri)
            V[i, i - 1] = -Rate[3]
            V[i, i] = Rate[3]

        # Next Generation Matrix M
        IvV = np.linalg.inv(V)  # Inverse of the drift matrix
        M = np.matmul(F, IvV)  # Next generation matrix

        # Eigenvalues
        EigVal, EigVec = np.linalg.eig(M)  # EigVal = Eigenvalues of the NGM, EigVec = Eigenvectors of the NGM
        Max_Eig = max(np.absolute(EigVal))  # Max_Eig = Max spectral radius

        # Coefficient
        coef = 1 / Max_Eig  # Coef used as a multiplication parameter for the mixing matrix to assure that the
        # Max_Eig = R0

        return coef
        #return [coef * xGG, coef * xGS, coef * xGR, coef * xSG, coef * xSS, coef * xSR, coef * xRG, coef * xRS,
        #        coef * xRR]

import numpy as np


class Field:

    def __init__(self, pulsar, anomalies):
        self.a = anomalies
        self.p = pulsar

    def r1(self, theta, phi, n):
        return self.a.r[n] * (np.sin(self.a.theta[n]) * np.sin(theta) * np.cos(phi - self.a.phi[n]) + np.cos(self.a.theta[n]) * np.cos(theta))

    def r2(self, theta, phi, n):
        return self.a.r[n] * (np.sin(self.a.theta[n]) * np.cos(theta) * np.cos(phi - self.a.phi[n]) - np.cos(self.a.theta[n]) * np.sin(theta))

    def r3(self, phi, n):
        return - self.a.r[n] * np.sin(self.a.theta[n]) * np.sin(phi - self.a.phi[n])

    def m1(self, theta, phi, n):
        return self.a.m[n] * (np.sin(self.a.theta_m[n]) * np.sin(theta) * np.cos(phi - self.a.phi_m[n]) + np.cos(self.a.theta_m[n]) * np.cos(theta))

    def m2(self, theta, phi, n):
        # changed sin(theta) -> cos(theta)
        return self.a.m[n] * (np.sin(self.a.theta_m[n]) * np.cos(theta) * np.cos(phi - self.a.phi_m[n]) - np.cos(self.a.theta_m[n]) * np.sin(theta))

    def m3(self, phi,n):
        return - self.a.m[n] * np.sin(self.a.theta_m[n]) * np.sin(phi - self.a.phi_m[n])

    # r Is The Local Name Of One Of The Values Of Matrix Rho The Same Is Used Also Bellow
    def T(self, r, theta, phi, n):
        return self.m1(theta, phi, n) * r - (self.m1(theta,phi,n) * self.r1(theta,phi,n) +
             self.m2(theta,phi,n) * self.r2(theta, phi, n) + self.m3(phi, n) * self.r3(phi, n))

    def D(self, r, theta, phi, n):
        return self.a.r[n] ** 2 + r ** 2 - 2. * self.a.r[n] * r * (np.sin(self.a.theta[n]) * np.sin(theta) * np.cos(phi - self.a.phi[n]) + np.cos(self.a.theta[n]) * np.cos(theta))

    def Bb1(self, r, theta, phi, n):
        # This Local Variable Is Used To Calculate D^-2,5
        DL = 1. / (self.D(r, theta, phi, n) ** 2.5)
        return - DL * (3. * self.T(r, theta, phi, n) * self.r1(theta, phi, n) - 3. * self.T(r, theta, phi, n) * r + self.D(r, theta, phi, n) * self.m1(theta, phi, n))

    def Bb2(self, r, theta, phi, n):
        # This Local Variable Is Used To Calculate D^-2.5
        DL = 1. / (self.D(r, theta, phi, n) ** 2.5)
        return - DL * (3. * self.T(r, theta, phi, n) * self.r2(theta, phi, n) + self.D(r, theta, phi, n) * self.m2(theta, phi, n))

    def Bb3(self, r, theta, phi, n):
        # This Local Variable Is Used To Calculate D^-2.5
        DL = 1. / (self.D(r, theta, phi, n) ** 2.5)
        return - DL * (3. * self.T(r, theta, phi, n) * self.r3(phi, n) + self.D(r, theta, phi, n) * self.m3(phi, n))

    def B1(self, r, theta, phi):
        bl = 0.
        for n in range(self.a.num):
            bl += self.Bb1(r, theta, phi, n)
        return bl

    def B2(self, r, theta, phi):
        bl = 0.
        for n in range(self.a.num):
            bl += self.Bb2(r, theta, phi, n)
        return bl

    def B3(self, r, theta, phi):
        bl = 0.0
        for n in range (self.a.num):
            bl += self.Bb3(r, theta, phi, n)
        return bl

    def H1(self, r, theta):
        return 2. * np.cos(theta) / (r ** 3)

    def H2(self, r, theta):
        return np.sin(theta) / (r ** 3)

    # Here are functions  for rotation axes  and angle business
    def O1(self, theta, phi):
        return np.sin(self.p.alpha) * np.sin(theta) * np.sin(phi) + np.cos(self.p.alpha) * np.cos(theta)

    def O2(self, theta, phi):
        return np.sin(self.p.alpha) * np.cos(theta) * np.sin(phi) - np.cos(self.p.alpha) * np.sin(theta)

    def O3(self, theta, phi):
        return np.sin(self.p.alpha) * np.cos(phi)

    def BFull(self, r, theta, phi):
        """
        Full magnetic field: sum of global and local dipole fields
        :param r:
        :param theta:
        :param phi:
        :return:
        """
        return np.sqrt((self.B1(r, theta, phi) + self.H1(r, theta)) ** 2 + (self.B2(r, theta, phi) + self.H2(r, theta)) ** 2 + self.B3(r, theta, phi) ** 2)

    def HFull(self, r, theta):
        return np.sqrt(self.H1(r, theta) ** 2 + self.H2(r, theta) ** 2)

    def fsigma(self, r, theta, phi):
        return (self.O1(theta, phi) * (self.B1(r, theta, phi) + self.H1(r, theta)) +
            self.O2(theta, phi) * (self.B2(r, theta, phi) + self.H2(r, theta)) + self.O3(theta, phi) * self.B3(r, theta, phi)) / self.BFull(r, theta, phi)

    def derivs (self, x, r):
        der_ = []
        der_.append((self.B2(r, x[0], x[1]) + self.H2(r, x[0])) / (r * (self.B1(r, x[0], x[1]) + self.H1(r, x[0]))))
        der_.append(self.B3(r, x[0], x[1]) / (r * np.sin(x[0]) * (self.B1(r, x[0], x[1]) + self.H1(r, x[0]))))
        return der_

    def dtheta_dphi(self, th_ph, r):
        """ ODEs to calculate line """
        dtheta_dr = (self.B2(r, th_ph[0], th_ph[1]) + self.H2(r, th_ph[0])) / (r * (self.B1(r, th_ph[0], th_ph[1]) + self.H1(r, th_ph[0])))
        dphi_dr = self.B3(r, th_ph[0], th_ph[1]) / (r * np.sin(th_ph[0]) * (self.B1(r, th_ph[0], th_ph[1]) + self.H1(r, th_ph[0])))
        return np.array([dtheta_dr, dphi_dr])


def main():
    print('Bye')

if __name__ == "__main__":
    main()

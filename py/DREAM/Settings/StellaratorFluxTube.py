import numpy as np

from DREAM.DREAMException import DREAMException

from .StellaratorGeometry import StellaratorGeometryPackage


def _evaluate_fourier(coeff_cos, coeff_sin, m, xn, theta, phi):
    phase = np.outer(m, theta) - np.outer(xn, phi)
    return np.sum(coeff_cos[:, None] * np.cos(phase) + coeff_sin[:, None] * np.sin(phase), axis=0)


def _evaluate_fourier_dtheta(coeff_cos, coeff_sin, m, xn, theta, phi):
    phase = np.outer(m, theta) - np.outer(xn, phi)
    return np.sum(-m[:, None] * coeff_cos[:, None] * np.sin(phase) + m[:, None] * coeff_sin[:, None] * np.cos(phase), axis=0)


def _evaluate_fourier_dphi(coeff_cos, coeff_sin, m, xn, theta, phi):
    phase = np.outer(m, theta) - np.outer(xn, phi)
    return np.sum(xn[:, None] * coeff_cos[:, None] * np.sin(phase) - xn[:, None] * coeff_sin[:, None] * np.cos(phase), axis=0)


class FluxTubeEvaluator:
    def __init__(self, geometry_package):
        if isinstance(geometry_package, StellaratorGeometryPackage):
            self.package = geometry_package
        else:
            self.package = StellaratorGeometryPackage.read(geometry_package)

        if self.package.boozer is None or "vmec" not in self.package.boozer:
            raise DREAMException(
                "FluxTubeEvaluator: This geometry package does not contain the spectral/Boozer-capable block required for field-line evaluation."
            )

        self.vmec = self.package.boozer["vmec"]
        self.s_grid = np.asarray(self.package.boozer.get("s", self.package.profiles["s"]), dtype=np.float64)
        self.iota = np.asarray(self.package.boozer.get("iota", self.package.profiles["iota"]), dtype=np.float64)
        self.m = np.asarray(self.vmec["xm"], dtype=np.float64)
        self.xn = np.asarray(self.vmec["xn"], dtype=np.float64)
        self.m_nyq = np.asarray(self.vmec["xm_nyq"], dtype=np.float64)
        self.xn_nyq = np.asarray(self.vmec["xn_nyq"], dtype=np.float64)

    def _surface_index(self, s=None, rho=None):
        if s is None and rho is None:
            raise DREAMException("FluxTubeEvaluator: Either 's' or 'rho' must be specified.")

        if rho is not None:
            s = (float(rho) / self.package.a) ** 2

        return int(np.argmin(np.abs(self.s_grid - float(s))))

    def _fieldline_rhs(self, surface_index, theta, phi):
        bsupu = _evaluate_fourier(
            np.asarray(self.vmec["bsupumnc"][surface_index], dtype=np.float64),
            np.asarray(self.vmec["bsupumns"][surface_index], dtype=np.float64),
            self.m_nyq,
            self.xn_nyq,
            np.asarray([theta]),
            np.asarray([phi]),
        )[0]
        bsupv = _evaluate_fourier(
            np.asarray(self.vmec["bsupvmnc"][surface_index], dtype=np.float64),
            np.asarray(self.vmec["bsupvmns"][surface_index], dtype=np.float64),
            self.m_nyq,
            self.xn_nyq,
            np.asarray([theta]),
            np.asarray([phi]),
        )[0]

        if abs(bsupv) < 1e-14:
            raise DREAMException("FluxTubeEvaluator: Encountered |B^phi| ~ 0 while tracing the field line.")

        return bsupu / bsupv

    def evaluate(self, *, s=None, rho=None, alpha=0.0, nturns=1, npoints=256):
        surface_index = self._surface_index(s=s, rho=rho)
        phi_end = float(nturns) * 2.0 * np.pi / float(self.package.nfp)
        phi = np.linspace(0.0, phi_end, int(npoints), endpoint=True)
        theta = np.zeros_like(phi)
        theta[0] = float(alpha)

        for i in range(1, phi.size):
            h = phi[i] - phi[i - 1]
            th = theta[i - 1]
            ph = phi[i - 1]

            k1 = self._fieldline_rhs(surface_index, th, ph)
            k2 = self._fieldline_rhs(surface_index, th + 0.5 * h * k1, ph + 0.5 * h)
            k3 = self._fieldline_rhs(surface_index, th + 0.5 * h * k2, ph + 0.5 * h)
            k4 = self._fieldline_rhs(surface_index, th + h * k3, ph + h)
            theta[i] = th + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

        rmnc = np.asarray(self.vmec["rmnc"][surface_index], dtype=np.float64)
        rmns = np.asarray(self.vmec["rmns"][surface_index], dtype=np.float64)
        zmns = np.asarray(self.vmec["zmns"][surface_index], dtype=np.float64)
        zmnc = np.asarray(self.vmec["zmnc"][surface_index], dtype=np.float64)
        lmns = np.asarray(self.vmec["lmns"][surface_index], dtype=np.float64)
        lmnc = np.asarray(self.vmec["lmnc"][surface_index], dtype=np.float64)

        bmnc = np.asarray(self.vmec["bmnc"][surface_index], dtype=np.float64)
        bmns = np.asarray(self.vmec["bmns"][surface_index], dtype=np.float64)
        gmnc = np.asarray(self.vmec["gmnc"][surface_index], dtype=np.float64)
        gmns = np.asarray(self.vmec["gmns"][surface_index], dtype=np.float64)

        bsupumnc = np.asarray(self.vmec["bsupumnc"][surface_index], dtype=np.float64)
        bsupumns = np.asarray(self.vmec["bsupumns"][surface_index], dtype=np.float64)
        bsupvmnc = np.asarray(self.vmec["bsupvmnc"][surface_index], dtype=np.float64)
        bsupvmns = np.asarray(self.vmec["bsupvmns"][surface_index], dtype=np.float64)

        R = _evaluate_fourier(rmnc, rmns, self.m, self.xn, theta, phi)
        Z = _evaluate_fourier(zmnc, zmns, self.m, self.xn, theta, phi)
        B = _evaluate_fourier(bmnc, bmns, self.m_nyq, self.xn_nyq, theta, phi)
        sqrtg = _evaluate_fourier(gmnc, gmns, self.m_nyq, self.xn_nyq, theta, phi)
        bsupu = _evaluate_fourier(bsupumnc, bsupumns, self.m_nyq, self.xn_nyq, theta, phi)
        bsupv = _evaluate_fourier(bsupvmnc, bsupvmns, self.m_nyq, self.xn_nyq, theta, phi)
        lambda_theta = _evaluate_fourier_dtheta(lmnc, lmns, self.m, self.xn, theta, phi)
        lambda_phi = _evaluate_fourier_dphi(lmnc, lmns, self.m, self.xn, theta, phi)
        R_theta = _evaluate_fourier_dtheta(rmnc, rmns, self.m, self.xn, theta, phi)
        Z_theta = _evaluate_fourier_dtheta(zmnc, zmns, self.m, self.xn, theta, phi)
        R_phi = _evaluate_fourier_dphi(rmnc, rmns, self.m, self.xn, theta, phi)
        Z_phi = _evaluate_fourier_dphi(zmnc, zmns, self.m, self.xn, theta, phi)

        g_tt = R_theta**2 + Z_theta**2
        g_tp = R_theta * R_phi + Z_theta * Z_phi

        return {
            "surface_index": surface_index,
            "s": float(self.s_grid[surface_index]),
            "rho": float(self.package.grid["rho"][min(surface_index, self.package.grid["rho"].size - 1)]),
            "theta": theta,
            "zeta": phi,
            "R": R,
            "Z": Z,
            "|B|": B,
            "sqrt(g)": sqrtg,
            "bsupu": bsupu,
            "bsupv": bsupv,
            "g_tt": g_tt,
            "g_tp": g_tp,
            "lambda_t": lambda_theta,
            "lambda_p": lambda_phi,
            "iota": float(self.iota[min(surface_index, self.iota.size - 1)]),
        }

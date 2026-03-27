# Implementation of stellarator magnetic equilibrium data.

import pathlib
import warnings

import matplotlib.pyplot as plt
import numpy as np

from DREAM.DREAMException import DREAMException

from .NumericalMagneticField import NumericalMagneticField
from .StellaratorFluxTube import FluxTubeEvaluator
from .StellaratorGeometry import StellaratorGeometryPackage
from .StellaratorGeometryProviders import DescProvider, PackageProvider, VmecJaxProvider, load_geometry_package


PROVIDER_DESC = "desc"
PROVIDER_VMEC_JAX = "vmec_jax"
PROVIDER_PACKAGE = "package"


class StellaratorMagneticField(NumericalMagneticField):
    def __init__(
        self,
        source,
        nr,
        ntheta=129,
        nphi=129,
        cache_filename=None,
        provider=None,
        write_cache=False,
        with_boozer=False,
    ):
        self.filename = str(source) if source is not None else None
        self.nr = nr
        self.ntheta = ntheta
        self.nphi = nphi
        self.cache_filename = cache_filename
        self.provider_name = provider
        self.write_cache = bool(write_cache)
        self.with_boozer = bool(with_boozer)

        self.package = None
        self.provider = None

        self.rho = None
        self.theta = None
        self.phi = None
        self.f_passing = None
        self.B_min = None
        self.B_max = None
        self.G = None
        self.I = None
        self.iota = None
        self.psi_T = None
        self.B = None
        self.BdotGradPhi = None
        self.Jacobian = None
        self.g_tt = None
        self.g_tp = None
        self.lambda_t = None
        self.lambda_p = None
        self.R = None
        self.Z = None
        self.R0 = None
        self.a = None
        self.nfp = None

        if self.provider_name is None:
            if self.filename and StellaratorGeometryPackage.is_geometry_package(self.filename):
                self.provider_name = PROVIDER_PACKAGE
            else:
                self.provider_name = PROVIDER_DESC

        if self.cache_filename is not None and load_geometry_package(self.cache_filename) is not None:
            self.package = load_geometry_package(self.cache_filename)
        elif self.cache_filename is not None:
            cache_path = pathlib.Path(self.cache_filename if str(self.cache_filename).endswith(".h5") else f"{self.cache_filename}.h5")
            if cache_path.is_file() and not StellaratorGeometryPackage.is_geometry_package(cache_path):
                self.package = StellaratorGeometryPackage.from_legacy_cache(
                    cache_path,
                    source_kind="legacy_cache",
                    provider="legacy_cache",
                    source_path=self.filename,
                )
        elif self.provider_name == PROVIDER_PACKAGE:
            self.package = StellaratorGeometryPackage.read(self.filename)

        if self.package is not None:
            self._hydrate_from_package(self.package)
        else:
            self.provider = self._make_provider(self.provider_name)

        super().__init__(self.a if self.a is not None else 0.0, self.R0 if self.R0 is not None else 0.0)

    def _make_provider(self, provider_name):
        provider_name = str(provider_name).lower()
        if provider_name == PROVIDER_DESC:
            return DescProvider(self.filename)
        if provider_name == PROVIDER_VMEC_JAX:
            return VmecJaxProvider(self.filename)
        if provider_name == PROVIDER_PACKAGE:
            return PackageProvider(self.filename)

        raise DREAMException(f"StellaratorMagneticField: Unrecognized stellarator geometry provider '{provider_name}'.")

    def _hydrate_from_package(self, package):
        self.package = package
        kernel_data = package.to_kernel_data()
        self.rho = kernel_data["rho"]
        self.theta = kernel_data["theta"]
        self.phi = kernel_data["phi"]
        self.f_passing = kernel_data["f_passing"]
        self.B_min = kernel_data["B_min"]
        self.B_max = kernel_data["B_max"]
        self.G = kernel_data["G"]
        self.I = kernel_data["I"]
        self.iota = kernel_data["iota"]
        self.psi_T = kernel_data["psi_T"]
        self.B = kernel_data["B"]
        self.BdotGradPhi = kernel_data["BdotGradPhi"]
        self.Jacobian = kernel_data["Jacobian"]
        self.g_tt = kernel_data["g_tt"]
        self.g_tp = kernel_data["g_tp"]
        self.lambda_t = kernel_data["lambda_t"]
        self.lambda_p = kernel_data["lambda_p"]
        self.R = kernel_data["R"]
        self.Z = kernel_data["Z"]
        self.R0 = float(kernel_data["R0"])
        self.a = float(kernel_data["a"])
        self.nfp = int(kernel_data["nfp"])
        self.provider_name = package.metadata.get("provider", self.provider_name)

    def load(self, savedata=None, savefilename=None):
        if self.package is None:
            self.package = self.provider.build_package(
                nr=self.nr,
                ntheta=self.ntheta,
                nphi=self.nphi,
                with_boozer=self.with_boozer,
            )
            self._hydrate_from_package(self.package)

        if savedata is None:
            savedata = self.write_cache
        if savefilename is None:
            savefilename = self.cache_filename

        if savedata and savefilename is not None:
            self.package.write(savefilename)

        return self.package

    def getMinorRadius(self):
        return float(self.a)

    def getMajorRadius(self):
        return float(self.R0)

    def getFluxTubeEvaluator(self):
        if self.package is None:
            self.load()
        return FluxTubeEvaluator(self.package)

    def visualize(self, phi, nrho=None, ntheta=None, ax=None, show=None):
        red = (249 / 255, 65 / 255, 68 / 255)
        black = (87 / 255, 117 / 255, 144 / 255)
        gray = (120 / 255, 120 / 255, 120 / 255)

        genax = ax is None
        if genax:
            ax = plt.axes()
            if show is None:
                show = True

        if np.isscalar(phi):
            phi = [phi]

        if nrho is None:
            nrho = self.nr - 1
        if ntheta is None:
            ntheta = self.ntheta

        if self.package is not None:
            R = self.package.reshape_sampled("R")
            Z = self.package.reshape_sampled("Z")
            for iphi in phi:
                idx = int(np.argmin(np.abs(self.phi - iphi)))
                ax.plot(R[idx, :, :].T[:, :-1], Z[idx, :, :].T[:, :-1], color=gray, linewidth=0.5)
                ax.plot(R[idx, :, :].T[:, -1], Z[idx, :, :].T[:, -1], color=red, linewidth=2)
                ax.plot(R[idx, 0, 0], Z[idx, 0, 0], "s", color=red)
        else:
            warnings.warn("StellaratorMagneticField: Geometry was not loaded before visualization.", RuntimeWarning)

        ax.plot(self.R0, 0, "X", color=black)
        ax.set_xlabel("$R$ (m)")
        ax.set_ylabel("$Z$ (m)")
        ax.axis("equal")

        if show:
            plt.show()

        return ax

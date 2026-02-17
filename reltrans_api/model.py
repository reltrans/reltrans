from dataclasses import dataclass
from typing import Optional, Any
import numpy as np

from .wrapper import Reltrans, DCP_Parameters


@dataclass
class ReltransModel:
    """
    High-level Python interface for reltrans DCP model.
    Provides a user-friendly API for spectrum computation.
    """

    path: Optional[str] = None
    engine: Optional[Any] = None

    def __post_init__(self):
        if self.engine is None:
            self.engine = Reltrans(self.path)

    def compute_spectrum(
        self,
        energy: np.ndarray,
        **params: Any,
    ) -> np.ndarray:
        """
        Compute DCP spectrum.
        """

        dcp_params = DCP_Parameters(**params)
        return self.engine.dcp(energy, dcp_params)


import numpy as np
from reltrans_api import ReltransModel


class FakeEngine:
    def dcp(self, energy, params):
        return np.full(len(energy) - 1, 7.0, dtype=np.float32)


def test_api_without_native():
    energy = np.linspace(1, 10, 11)

    model = ReltransModel(engine=FakeEngine())
    output = model.compute_spectrum(energy)

    assert output.shape == (10,)
    assert (output == 7.0).all()


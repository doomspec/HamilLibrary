from typing import Tuple


class MeasurementMethod:
    def get_groups(self, H, options=None) -> Tuple[list, list]:
        pass

class QubitMeasurementMethod(MeasurementMethod):
    def get_groups(self, qubit_H, options=None):
        pass

class LargeQubitMeasurementMethod(MeasurementMethod):
    def get_groups(self, qubit_H, options=None):
        pass
"""
GCS module for mitochondrial 1C metabolism model.
Full reaction representation using fixed literature flux.
"""

from ModelParameters import *
class GlycineCleavageSystem:

    def __init__(self):
        self.ran = False  # Flag to indicate if the model has been run

    def mito_gcs_step(self, pool, step_timespan):
        """

        GLY + THF + NADP+ → CH2THF + CO2 + NH3 + NADPH

        Uses fixed flux (literature ~0.3 mM/h = 300 µM/h)
        """

        print(" mito_gcs")
        print(" gly:", pool[GLY],
              " thf:", pool[THF],
              " nadpp:", pool[NADPP],
              " ch2thf:", pool[CH2THF],
              " co2:", pool[CO2])

        # Fixed flux (µM/h)
        V_GCS = 300.0

        dt_hr = step_timespan
        delta = V_GCS * dt_hr

        # Substrate consumption
        pool[GLY] -= delta
        pool[THF] -= delta
        pool[NADPP] -= delta

        # Product formation
        pool[CH2THF] += delta
        pool[CO2] += delta
        pool[NADPH] += delta
        pool[H_PLUS] += delta

        print(" gcs_flux:", V_GCS)
        print(" gly:", pool[GLY],
              " thf:", pool[THF],
              " ch2thf:", pool[CH2THF],
              " nadph:", pool[NADPH],
              " co2:", pool[CO2])

        return V_GCS

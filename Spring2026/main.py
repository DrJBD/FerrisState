"""Main entry point for the conversion models"""
import sys
from pprint import pprint
from ModelParameters import *
from SerToFormateMmModel import SerToFormateMmModel
from ReportingTable import report

def main():
    """
    Main entry point for the mitochondrial serine-to-formate model.
    """

    experiment_runtime_in_sec = (
        EXPERIMENT_RUNTIME_IN_HOURS * 3600
    )

    model = SerToFormateMmModel(
        MITO_POOL.copy(),
        experiment_runtime_in_sec,
        EXPERIMENT_TIMESTEPS
    )

    result_pool = model.run()

    report.pprint()
    report.save_plots()

    print(
        'Formate produced from',
        INITIAL_SER_CONCENTRATION,
        'µM SER:',
        result_pool[FORMATE],
        'µM'
    )


if __name__ == "__main__":
    main()


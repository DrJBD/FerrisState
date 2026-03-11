"""Main entry point for the conversion models"""
import sys
from pprint import pprint
from ModelParameters import *
from SerToFormateMmModel import SerToFormateMmModel
from ReportingTable import report

def main():
    """
    This is the main function that will be executed when the script is run.
    """
    model = SerToFormateMmModel(MITO_POOL.copy(), EXPERIMENT_RUNTIME_IN_HOURS, EXPERIMENT_TIMESTEPS)
    result_pool = model.run()

    report.pprint()
    report.save_plots()
    print('Formate produced from', INITIAL_SER_CONCENTRATION, 'M SER:', result_pool[FORMATE], 'M')

# Execute the main function
main()

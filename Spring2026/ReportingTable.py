"""
Store the intermediate values of the simulation for end-of-run report generation.
"""
import csv
import os
import pandas
import matplotlib.pyplot as plt
from ModelParameters import *
from pprint import pprint

class ReportingTable:

    def __init__(self):
        self.log = {}


    def add_step(self, step, values):
        """
        Add key/values as a step in the log.  The values will merge into an existing entry if the step
        has already been created.  New values take priority.
        """
        to_add = values.copy()
        if step in self.log:
            to_add = to_add | self.log[step]
        self.log[step] = to_add


    def pprint(self):
        """
        Pretty print the self.log variable and anything else we need.
        """
        print(f"ReportingTable: KM_SER_MAX: {KM_SER_MAX}")
        print(f"ReportingTable: ALPHA_KM_SER: {ALPHA_KM_SER}")
        print(f"ReportingTable: K05_THF: {K05_THF}")
        print(f"ReportingTable: EXPERIMENT_RUNTIME_IN_HOURS: {EXPERIMENT_RUNTIME_IN_HOURS}")
        print(f"ReportingTable: EXPERIMENT_TIMESTEPS: {EXPERIMENT_TIMESTEPS}")
        print(f"ReportingTable: INITIAL_SER_CONCENTRATION: {INITIAL_SER_CONCENTRATION}")

        if not self.log:
            print('ReportingTable: No steps were added to the log')
        else:
            pprint(f'ReportingTable: Min step {min(self.log)}.  Max step {max(self.log)}')
            # pprint(self.log)


    def save_csv(self, filename):
        # Get all the value names stored even if they are not in every step.  Will enter blanks in the CSV
        columns = set([v for values in self.log.values() for v in values])

        with open(filename, 'w', newline='') as f:
            fieldnames = ['step'] + sorted(columns)
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()

            for step, concentrations in self.log.items():
                row_cells = dict.fromkeys(fieldnames, '')
                row_cells['step'] = step
                for key in concentrations:
                    row_cells[key] = concentrations[key]
                writer.writerow(row_cells)


    def save_plots(self):
        # Load the CSV
        filename = 'Spring2026_report.csv'
        self.save_csv(filename)
        df = pandas.read_csv(filename)

        # First column is the x-axis (step)
        x_col = df.columns[0]
        data_cols = df.columns[1:]

        # Output directory for plots
        output_dir = "plots"
        os.makedirs(output_dir, exist_ok=True)

        for col in data_cols:
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(df[x_col], df[col], linewidth=1.5)
            ax.set_xlabel(x_col.capitalize(), fontsize=12)
            ax.set_ylabel(col, fontsize=12)
            ax.set_title(col, fontsize=14)
            ax.grid(True, linestyle="--", alpha=0.5)
            plt.tight_layout()
            # Sanitize filename
            safe_name = col.replace("/", "_").replace("+", "plus").replace(" ", "_")
            fig.savefig(os.path.join(output_dir, f"{safe_name}.png"), dpi=150)
            plt.close(fig)
            print(f"Saved plot for: {col}")

        print(f"\nDone! {len(data_cols)} plots saved to '{output_dir}/'")

report = ReportingTable()
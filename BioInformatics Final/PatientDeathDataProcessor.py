import csv
import numpy as np


class PatientDeathDataProcessor:
    """
    This class processes clinical data to create a standardized mapping of patient IDs to their
    normalized days-to-death values. It handles loading the data, mapping it, and normalizing it.
    """

    def __init__(self, cancer_type):
        """
        Initializes the processor with the specified cancer type. Sets up containers for patient IDs,
        raw days-to-death data, and the final standardized death data.
        """
        self.cancer_type = cancer_type
        self.patient_ids = []
        self.days_to_death_data = []
        self.standardized_death_data = {}

    def _load_data(self):
        """
        Loads patient IDs and days-to-death data from a clinical data file specific to the cancer type.
        The data file is expected to be in tab-separated values (TSV) format.
        """
        line_num = 0
        file_name = f"{self.cancer_type}.clin.merged.picked.txt"
        with open(file_name) as tsv:
            for line in csv.reader(tsv, dialect="excel-tab"):
                if line_num == 0:
                    # Extract patient identifiers from the first row
                    self.patient_ids = line[1:]
                elif line_num > 1 and line[0] == "days_to_death":
                    # Extract days to death from the relevant row
                    self.days_to_death_data = line[1:]
                line_num += 1

    def _map_patients_to_death(self):
        """
        Maps patient codes (derived from patient IDs) to their corresponding days-to-death values.
        Filters out patients with no recorded death (where days_to_death is "NA").
        """
        patient_codes = [patient_id.split("-")[2] for patient_id in self.patient_ids]
        self.standardized_death_data = {
            patient_codes[i]: float(self.days_to_death_data[i]) for i in range(len(self.days_to_death_data)) if
            self.days_to_death_data[i] != "NA"
        }

    def _normalize_death_data(self):
        """
        Normalizes the days-to-death data by calculating the mean and standard deviation.
        Then, it adjusts each patient's days-to-death value to a z-score, which represents how many
        standard deviations the value is from the mean.
        """
        death_values = list(self.standardized_death_data.values())
        mean_death = np.mean(death_values)
        std_dev_death = np.std(death_values)

        print("AVERAGE DEATH ", mean_death)
        print("STD DEATH", std_dev_death)

        # Normalize the days to death values to z-scores
        self.standardized_death_data = {
            patient: (self.standardized_death_data[patient] - mean_death) / std_dev_death for patient in self.standardized_death_data.keys()
        }

    def get_standardized_patients_to_death(self):
        """
        Orchestrates the data processing by calling the load, map, and normalize methods in sequence.
        Returns the final dictionary mapping patient IDs to their standardized days-to-death values.
        """
        self._load_data()
        self._map_patients_to_death()
        self._normalize_death_data()
        return self.standardized_death_data

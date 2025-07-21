from pathlib import Path
from libra_toolbox.neutron_detection.activation_foils.calibration import (
    CheckSource,
    co60,
    cs137,
    mn54,
    na22,
    ActivationFoil,
    nb93_n2n,
    zr90_n2n,
)
import libra_toolbox.neutron_detection.activation_foils.compass as compass
from libra_toolbox.neutron_detection.activation_foils.compass import (
    Measurement,
    CheckSourceMeasurement,
    SampleMeasurement,
)
from libra_toolbox.tritium.model import ureg
from datetime import date, datetime
import json
from zoneinfo import ZoneInfo
from download_raw_foil_data import download_and_extract_foil_data
import copy

#####################################################################
##################### CHANGE THIS FOR EVERY RUN #####################
#####################################################################
measurement_directory = "BABY_1L_Run6_250530"

# Path to save the extracted files
output_path = Path("../../data/neutron_detection/")
activation_foil_path = output_path / "activation_foils"  

measurement_directory_path = activation_foil_path / measurement_directory


################ Check Source Calibration Information ###################
uCi_to_Bq = 3.7e4

#####################################################################
##################### CHANGE THIS FOR EVERY RUN #####################
#####################################################################

co60_checksource = CheckSource(
    co60, activity_date=date(2014, 3, 19), activity=0.872 * uCi_to_Bq
)
cs137_checksource = CheckSource(
    cs137, activity_date=date(2023, 9, 29), activity=9.38 * uCi_to_Bq
)
mn54_checksource = CheckSource(
    mn54, activity_date=date(2016, 5, 2), activity=6.27 * uCi_to_Bq
)
na22_checksource = CheckSource(
    na22, activity_date=date(2023, 9, 29), activity=9.98 * uCi_to_Bq
)
check_source_dict = {
    "Co60 Count 1": {
        "directory": measurement_directory_path
        / "Co60_0_872uCi_2014Mar19_count1/UNFILTERED",
        "check_source": co60_checksource,
    },
    "Co60 Count 2": {
        "directory": measurement_directory_path
        / "Co60_0_872uCi_2014Mar19_count2/UNFILTERED",
        "check_source": co60_checksource,
    },
    "Cs137 Count 1": {
        "directory": measurement_directory_path
        / "Cs137_9_38uCi_2023Sep29_count1/UNFILTERED",
        "check_source": cs137_checksource,
    },
    "Cs137 Count 2": {
        "directory": measurement_directory_path
        / "Cs137_9_38uCi_2023Sep29_count2/UNFILTERED",
        "check_source": cs137_checksource,
    },
    "Mn54 Count 1": {
        "directory": measurement_directory_path
        / "Mn54_6_27uCi_2016May2_count1/UNFILTERED",
        "check_source": mn54_checksource,
    },
    "Mn54 Count 2": {
        "directory": measurement_directory_path 
        / "Mn54_6_27uCi_2016May2_count2/UNFILTERED",
        "check_source": mn54_checksource,
    },
    "Na22 Count 1": {
        "directory": measurement_directory_path 
        / "Na22_9_98uCi_2023Sep29_count1/UNFILTERED",
        "check_source": na22_checksource,
    },
    "Na22 Count 2": {
        "directory": measurement_directory_path 
        / "Na22_9_98uCi_2023Sep29_count2/UNFILTERED",
        "check_source": na22_checksource,
    },
}

background_dir = measurement_directory_path / "Background_20250602_count1/UNFILTERED"

################ Foil Information ###################

def get_distance_to_source_from_dict(foil_dict: dict):
    distance_to_source_dict = foil_dict["distance_to_source"]
    # unit from string with pint
    unit = ureg.parse_units(distance_to_source_dict["unit"])
    return (distance_to_source_dict["value"] * unit).to(ureg.cm).magnitude
    

def get_mass_from_dict(foil_dict: dict):
    foil_mass = foil_dict["mass"]["value"]
    # unit from string with pint
    unit = ureg.parse_units(foil_dict["mass"]["unit"])
    return (foil_mass * unit).to(ureg.g).magnitude
    

def get_thickness_from_dict(foil_dict: dict):
    foil_thickness = foil_dict["thickness"]["value"]
    # unit from string with pint
    unit = ureg.parse_units(foil_dict["thickness"]["unit"])
    return (foil_thickness * unit).to(ureg.cm).magnitude


def get_foil(foil_element_symbol, foil_designator=None):
    """Get information about a specific foil from the general data file.
    Args:
        foil_element_symbol (str): The chemical symbol of the foil element (e.g., "Zr" for Zirconium).
        foil_designator (str, optional): The designator of the foil (e.g., "Nb Packet #1")
    Returns:
        ActivationFoil: An ActivationFoil object containing the foil's properties.
        distance_to_source (float): The distance from the foil to the neutron source in cm.
    """

    with open("../../data/general.json", "r") as f:
        general_data = json.load(f)
        foils = general_data["neutron_detection"]["foils"]["materials"]
        foil_of_specified_element_count = 0
        for foil in foils:
            if foil["material"] == foil_element_symbol:
                foil_of_specified_element_count += 1
                # if no foil_designator is provided, or if it matches the foil's designator
                if foil_designator is None or foil["designator"] == foil_designator:
                    # Get distance to generator
                    distance_to_source = get_distance_to_source_from_dict(foil)

                    # Get mass
                    foil_mass = get_mass_from_dict(foil)

                    # get foil thickness
                    foil_thickness = get_thickness_from_dict(foil)

                    # Get foil name
                    foil_name = foil["designator"]
                    if foil_name is None:
                        foil_name = foil_element_symbol

    if foil_of_specified_element_count == 0:
        raise ValueError(
            f"No foils found for element {foil_element_symbol} with designator {foil_designator}"
        )
    elif foil_of_specified_element_count > 1:
        print(
            f"Warning: Multiple foils found for element {foil_element_symbol} with designator {foil_designator}. Using the last one found."
        )

    if foil_element_symbol == "Zr":
        foil_density = 6.505
        foil_mass_attenuation_coefficient = 0.08590  # cm^2/g at 1 MeV
        foil_reaction = zr90_n2n
    elif foil_element_symbol == "Nb":
        foil_density = 8.582
        foil_mass_attenuation_coefficient = 0.06120  # cm^2/g at 1 MeV
        foil_reaction = nb93_n2n
    else:
        raise ValueError(f"Unsupported foil element symbol: {foil_element_symbol}")
    foil = ActivationFoil(
        reaction=foil_reaction,
        mass=foil_mass,
        name=foil_name,
        density=foil_density,
        thickness=foil_thickness,  # in cm
    )
    foil.mass_attenuation_coefficient = foil_mass_attenuation_coefficient
    print(f"Read in properties of {foil.name} foil")
    return foil, distance_to_source

#####################################################################
##################### CHANGE THIS FOR EVERY RUN #####################
#####################################################################

niobium_foil, nb_distance_to_source = get_foil("Nb")
zirconium_foil, zr_distance_to_source = get_foil("Zr")

foil_source_dict = {
    niobium_foil.name: {
        "measurement_paths": {
            1: measurement_directory_path
                / "Niobium3_20250601_1358_count1/UNFILTERED",
            2: measurement_directory_path
                / "Niobium3_20250602_1123_count2/UNFILTERED"
        },
        "foil": niobium_foil,
        "distance_to_source": nb_distance_to_source
    },
    zirconium_foil.name: {
        "measurement_paths": {
            1: measurement_directory_path
                / "Zirconium1_20250530_1512_count1/UNFILTERED",
            2: measurement_directory_path
                / "Zirconium1_20250531_2115_count2/UNFILTERED"
        },
        "foil": zirconium_foil,
        "distance_to_source": zr_distance_to_source
    }
}


def get_data(download_from_raw=False, url=None,
             check_source_dict=check_source_dict,
             background_dir=background_dir,
             foil_source_dict=foil_source_dict):

    if download_from_raw:
        if url is None:
            raise ValueError("Web URL not provided for downloading data.")
        download_and_extract_foil_data(url, activation_foil_path)
        # Process data
        check_source_measurements, background_meas = read_checksources_from_directory(
                                        check_source_dict, background_dir
                                        )
        foil_measurements = read_foil_measurements_from_dir(foil_source_dict)

        # save spectra to h5 for future, faster use
        print('Saving spectra to h5...')
        save_measurements(check_source_measurements,
                          background_meas,
                          foil_measurements)
    else:
        measurements = Measurement.from_h5(activation_foil_path / "activation_data.h5")
        foil_measurements = copy.deepcopy(foil_source_dict)
        check_source_measurements = {}
        # Get list of foil measurement names
        foil_measurement_names = []
        for foil_name in foil_source_dict.keys():
            for count_num in foil_source_dict[foil_name]["measurement_paths"]:
                foil_measurement_names.append(f"{foil_name} Count {count_num}")

            # Add empty measurements dictionary to foil_source_dict copy
            foil_measurements[foil_name]["measurements"] = {}
            
        for measurement in measurements:
            print(f"Processing {measurement.name} from h5 file...")
            # check if measurement is a check source measurement
            if measurement.name in check_source_dict.keys():
                # May want to change CheckSourceMeasurement in libra-toolbox to make this more seemless
                check_source_meas = CheckSourceMeasurement(measurement.name)
                check_source_meas.__dict__.update(measurement.__dict__)
                check_source_meas.check_source = check_source_dict[measurement.name]["check_source"]
                check_source_measurements[measurement.name] = check_source_meas
            elif measurement.name == "Background":
                background_meas = measurement
            elif measurement.name in  foil_measurement_names:
                # Extract foil name and count number from measurement name
                split_name = measurement.name.split(' ')
                count_num = int(split_name[-1])
                foil_name = " ".join(split_name[:-2])

                foil_meas = SampleMeasurement(measurement)
                foil_meas.__dict__.update(measurement.__dict__)
                foil_meas.foil = foil_source_dict[foil_name]["foil"]
                foil_measurements[foil_name]["measurements"][count_num] = foil_meas
            else:
                print(f"Extra measurement included in h5 file: {measurement.name}")
        
    return check_source_measurements, background_meas, foil_measurements


def save_measurements(check_source_measurements,
                      background_meas,
                      foil_measurements):
    
    measurements = []
    for name, measurement in check_source_measurements.items():
        measurements.append(measurement)
    measurements.append(background_meas)
    for foil_name in foil_measurements.keys():
        for count_num in foil_measurements[foil_name]["measurements"].keys():
            measurements.append(foil_measurements[foil_name]["measurements"][count_num])
    
    for i,measurement in enumerate(measurements):
        if i==0:
            mode = 'w'
        else:
            mode = 'a'
        measurement.to_h5(
            filename= activation_foil_path / "activation_data.h5",
            mode=mode,
            spectrum_only=True
        )


def read_checksources_from_directory(
    check_source_measurements: dict, background_dir: Path
):

    measurements = {}

    for name, values in check_source_measurements.items():
        print(f"Processing {name}...")
        meas = CheckSourceMeasurement.from_directory(values["directory"], name=name)
        meas.check_source = values["check_source"]
        measurements[name] = meas

    print(f"Processing background...")
    background_meas = Measurement.from_directory(
        background_dir,
        name="Background",
        info_file_optional=True,
    )
    return measurements, background_meas


def read_foil_measurements_from_dir(
    foil_measurements: dict
):

    for foil_name in foil_measurements.keys():
        foil_measurements[foil_name]["measurements"] = {}
        foil = foil_measurements[foil_name]["foil"]
        for count_num, measurement_path in foil_measurements[foil_name]["measurement_paths"].items():
            measurement_name = f"{foil_name} Count {count_num}"
            print(f"Processing {measurement_name}...")
            measurement = SampleMeasurement.from_directory(
                source_dir=measurement_path,
                name=measurement_name
            )
            measurement.foil = foil
            foil_measurements[foil_name]["measurements"][count_num] = measurement

    return foil_measurements


# Get the irradiation schedule

with open("../../data/general.json", "r") as f:
    general_data = json.load(f)
irradiations = []
for generator in general_data["generators"]:
    if generator["enabled"] is False:
        continue
    for i, irradiation_period in enumerate(generator["periods"]):
        if i == 0:
            overall_start_time = datetime.strptime(
                irradiation_period["start"], "%m/%d/%Y %H:%M"
            )
        start_time = datetime.strptime(irradiation_period["start"], "%m/%d/%Y %H:%M")
        end_time = datetime.strptime(irradiation_period["end"], "%m/%d/%Y %H:%M")
        irradiations.append(
            {
                "t_on": (start_time - overall_start_time).total_seconds(),
                "t_off": (end_time - overall_start_time).total_seconds(),
            }
        )
time_generator_off = end_time
time_generator_off = time_generator_off.replace(tzinfo=ZoneInfo("America/New_York"))



def calculate_neutron_rate_from_foil(foil_measurements, 
                                     foil_name,
                                     background_meas,
                                     calibration_coeffs,
                                     efficiency_coeffs,
                                     search_width=330,
                                     irradiations=irradiations,
                                     time_generator_off=time_generator_off):
    neutron_rates = {}
    neutron_rate_errs = {}

    for count_num, measurement in foil_measurements[foil_name]["measurements"].items():

        neutron_rates[f"Count {count_num}"] = {}
        neutron_rate_errs[f"Count {count_num}"] = {}

        for detector in measurement.detectors:
            ch = detector.channel_nb

            gamma_emitted, gamma_emitted_err = measurement.get_gamma_emitted(
                background_measurement=background_meas,
                calibration_coeffs=calibration_coeffs[ch],
                efficiency_coeffs=efficiency_coeffs[ch],
                channel_nb=ch,
                search_width=search_width)
            
            neutron_rate = measurement.get_neutron_rate(
                channel_nb=ch,
                photon_counts=gamma_emitted,
                irradiations=irradiations,
                distance=foil_measurements[foil_name]["distance_to_source"],
                time_generator_off=time_generator_off,
                branching_ratio=foil_measurements[foil_name]["foil"].reaction.product.intensity
            )

            neutron_rate_err = measurement.get_neutron_rate(
                channel_nb=ch,
                photon_counts=gamma_emitted_err,
                irradiations=irradiations,
                distance=foil_measurements[foil_name]["distance_to_source"],
                time_generator_off=time_generator_off,
                branching_ratio=foil_measurements[foil_name]["foil"].reaction.product.intensity
            )
            neutron_rates[f"Count {count_num}"][ch] = neutron_rate
            neutron_rate_errs[f"Count {count_num}"][ch] = neutron_rate_err

    return neutron_rates, neutron_rate_errs

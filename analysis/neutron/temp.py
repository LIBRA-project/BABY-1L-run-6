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
from datetime import date
import json


def get_foil(foil_element_symbol, foil_designator=None):
    """Get information about a specific foil from the general data file.
    Args:
        foil_element_symbol (str): The chemical symbol of the foil element (e.g., "Zr" for Zirconium).
        foil_designator (str, optional): The designator of the foil (e.g., "Nb Packet #1")
    Returns:
        ActivationFoil: An ActivationFoil object containing the foil's properties.
        distance_to_source (float): The distance from the foil to the neutron source in cm.
    """
    inches_to_cm = 2.54

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
                    distance_to_source_dict = foil["distance_to_source"]
                    if distance_to_source_dict["unit"] == "cm":
                        distance_to_source = distance_to_source_dict["value"]
                    else:
                        raise ValueError(
                            f"Unsupported unit for distance to source: {distance_to_source_dict['unit']}"
                        )

                    # Get mass
                    foil_mass = foil["mass"]["value"]
                    if foil["mass"]["unit"] != "g":
                        raise ValueError(
                            f"Unsupported unit for mass: {foil['mass']['unit']}"
                        )

                    # get foil thickness
                    foil_thickness = foil["thickness"]["value"]
                    if (
                        foil["thickness"]["unit"] == "inch"
                        or foil["thickness"]["unit"] == "in"
                    ):
                        foil_thickness *= inches_to_cm
                    elif foil["thickness"]["unit"] == "cm":
                        pass  # already in cm
                    else:
                        raise ValueError(
                            f"Unsupported unit for thickness: {foil['thickness']['unit']}"
                        )

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


# Path to save the extracted files
output_dir = Path("../../data/neutron_detection/")
output_file = output_dir / "foil_data.zip"
extracted_directory = output_dir / "activation_foils"

measurement_directory = extracted_directory / "BABY_1L_Run6_250530"

uCi_to_Bq = 3.7e4

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

check_source_measurements = {
    "Co60_1": {
        "directory": measurement_directory
        / "Co60_0_872uCi_2014Mar19_count1/UNFILTERED",
        "check_source": co60_checksource,
    },
    "Co60_2": {
        "directory": measurement_directory
        / "Co60_0_872uCi_2014Mar19_count2/UNFILTERED",
        "check_source": co60_checksource,
    },
    "Cs137_1": {
        "directory": measurement_directory
        / "Cs137_9_38uCi_2023Sep29_count1/UNFILTERED",
        "check_source": cs137_checksource,
    },
    "Cs137_2": {
        "directory": measurement_directory
        / "Cs137_9_38uCi_2023Sep29_count2/UNFILTERED",
        "check_source": cs137_checksource,
    },
    "Mn54_1": {
        "directory": measurement_directory / "Mn54_6_27uCi_2016May2_count1/UNFILTERED",
        "check_source": mn54_checksource,
    },
    "Mn54_2": {
        "directory": measurement_directory / "Mn54_6_27uCi_2016May2_count2/UNFILTERED",
        "check_source": mn54_checksource,
    },
    "Na22_1": {
        "directory": measurement_directory / "Na22_9_98uCi_2023Sep29_count1/UNFILTERED",
        "check_source": na22_checksource,
    },
    "Na22_2": {
        "directory": measurement_directory / "Na22_9_98uCi_2023Sep29_count2/UNFILTERED",
        "check_source": na22_checksource,
    },
}

background_dir = measurement_directory / "Background_20250602_count1/UNFILTERED"

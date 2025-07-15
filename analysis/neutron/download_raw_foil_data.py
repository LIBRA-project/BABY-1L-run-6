from pathlib import Path
import zipfile
import requests

from temp import (
    extracted_directory,
    output_file,
    check_source_measurements,
    background_dir,
)

from libra_toolbox.neutron_detection.activation_foils.compass import (
    Measurement,
    CheckSourceMeasurement,
    SampleMeasurement,
)


def download_and_extract_foil_data(extracted_directory: Path, output_file: Path):

    if extracted_directory.exists():
        print(f"Directory already exists: {extracted_directory}")
    else:
        # URL of the file
        url = "https://zenodo.org/records/15794193/files/BABY_1L_Run6_250530.zip?download=1"

        # Download the file
        print(f"Downloading data from {url}...")
        response = requests.get(url)
        if response.status_code == 200:
            print("Download successful!")
            # Save the file to the specified directory
            with open(output_file, "wb") as f:
                f.write(response.content)
            print(f"File saved to: {output_file}")
        else:
            print(f"Failed to download file. HTTP Status Code: {response.status_code}")

        # Extract the zip file

        # Ensure the extraction directory exists
        extracted_directory.mkdir(parents=True, exist_ok=True)

        # Unzip the file
        with zipfile.ZipFile(output_file, "r") as zip_ref:
            zip_ref.extractall(extracted_directory)
        print(f"Files extracted to: {extracted_directory}")

        # Delete the zip file after extraction
        output_file.unlink(missing_ok=True)


def read_measurement_from_directory(
    check_source_measurements: dict, background_dir: Path
):

    all_measurements = {}

    for name, values in check_source_measurements.items():
        print(f"Processing {name}...")
        meas = CheckSourceMeasurement.from_directory(values["directory"], name=name)
        meas.check_source = values["check_source"]
        all_measurements[name] = meas

    print(f"Processing background...")
    background_meas = Measurement.from_directory(
        background_dir,
        name="Background",
        info_file_optional=True,
    )
    return all_measurements, background_meas


download_and_extract_foil_data(extracted_directory, output_file)


check_source_measurements, background_meas = read_measurement_from_directory(
    check_source_measurements, background_dir
)

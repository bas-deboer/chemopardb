import os
import requests
import time

def handle_protoss_request(pdb_code, output_directory):
    # URLs and headers setup
    post_url = "https://proteins.plus/api/protoss_rest"
    headers = {
        "Accept": "application/json",
        "Content-Type": "application/json"
    }
    post_data = {
        "protoss": {"pdbCode": pdb_code}
    }

    # Initial POST request to check or start the job
    post_response = requests.post(post_url, json=post_data, headers=headers)
    post_result = post_response.json()

    # Check for the location URL in the response
    if 'location' in post_result:
        location_url = post_result['location']
        while True:
            # GET request to fetch the results
            get_response = requests.get(location_url)
            if get_response.status_code == 200:
                get_result = get_response.json()
                # Download files if URLs are provided in the results
                if 'protein' in get_result:
                    download_file(get_result['protein'], os.path.join(output_directory, f"{pdb_code}_protein.pdb"))
                if 'ligands' in get_result:
                    download_file(get_result['ligands'], os.path.join(output_directory, f"{pdb_code}_ligands.sdf"))
                return "Files downloaded successfully."
            elif get_response.status_code == 202:
                get_result = get_response.json()
                if 'message' in get_result:
                    print(get_result['message'])
                # Wait before retrying
                time.sleep(10)
            else:
                return f"Failed to retrieve job results: Status code {get_response.status_code}"
    else:
        return "Job location URL not found in the initial response."

def download_file(url, file_path):
    # Function to download and save a file from a given URL
    response = requests.get(url)
    if response.status_code == 200:
        with open(file_path, 'wb') as file:
            file.write(response.content)
        print(f"{file_path} saved successfully.")
    else:
        print(f"Failed to download {file_path}. Status code: {response.status_code}")

# Example usage
# pdb_code = "YOUR_PDB_CODE"
# output_directory = "YOUR_OUTPUT_DIRECTORY"
# handle_protoss_request(pdb_code, output_directory)

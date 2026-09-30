import urllib.request
import os
import time

# Your updated list of PDB IDs (Duplicates removed and formatted)
pdb_ids = [
    "1AQK", "1AD0", "1AD9", "1B2W", "1BEY", "1BFO", "1BVL", "1DEE", "1DFB",
    "1DN0", "1DQL", "1FGV", "1FL5", "1FN4", "1FVC", "1FVD", "1FVE", "1HEZ",
    "1HZH", "1IGM", "1IGT", "1IT9", "1JPT", "1L7I", "1MCO", "1QLR", "1R24",
    "1R70", "1RHH", "1RZ7", "1S3K", "1T04", "1T3F", "1U6A", "2A9M", "2AGJ",
    "2D7T", "2EH7", "2F5A", "2FB4", "2FGW", "2G75", "2GCY", "2HFF", "2IG2",
    "2JIX", "2QTJ", "2V7N", "2VXV", "2XA8", "2XZA", "2XZC", "3AAZ", "3CHN",
    "3CM9", "3D69", "3EO9", "3EOT", "3F12", "3FZU", "3G6A", "3HC3", "4C2I",
    "4DN3", "4F57", "4FQC", "4FZ8", "4GFZ", "4HH9", "4JY6", "4KTD", "4LLV",
    "4NUG", "4OSU", "4R7N", "4UOK", "4XML", "5AWN", "5BZD", "5CEX", "5CGY",
    "5DK3", "5DR5", "5F7E", "5GS1", "5H32", "5JQD", "5K8A", "5KG9", "1N8Z",
    "5WT9", "4G5Z", "7K8Q", "3IDG", "5GGQ", "8FAB", "6QN8", "6NCP", "6MTS",
    "6OBD", "4HCR", "6CR1", "7D85", "6R8X", "7MXL", "6XY2", "8JEL", "1CZ8",
    "4XXD", "4R7D", "3BKY", "4EDW", "3QWO"
]

# Create a directory to store the files
output_dir = "Input_Fasta"
os.makedirs(output_dir, exist_ok=True)

base_url = "https://www.rcsb.org/fasta/entry/"

print(f"Starting download of {len(pdb_ids)} FASTA files...")

for pdb_id in pdb_ids:
    url = f"{base_url}{pdb_id}"
    file_path = os.path.join(output_dir, f"{pdb_id}.fasta")

    try:
        print(f"Downloading {pdb_id}...", end=" ", flush=True)
        # Fetch and save the file
        urllib.request.urlretrieve(url, file_path)
        print("Done.")
        # Pause briefly to be polite to the PDB server
        time.sleep(0.1)
    except Exception as e:
        print(f"Failed! Error: {e}")

print(f"\nAll downloads complete! Check the '{output_dir}' folder.")
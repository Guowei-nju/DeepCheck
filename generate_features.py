import subprocess
import sys

def extract_features(input_folder, output_folder, threads=30):
    command = [
        "checkm2", "predict",
        "--threads", str(threads),
        "--input", input_folder,
        "--output-directory", output_folder
    ]
    result = subprocess.run(command, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error: {result.stderr}", file=sys.stderr)
    else:
        print(f"Output: {result.stdout}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python run_checkm2.py <folder_with_bins> <output_folder>")
        sys.exit(1)
    
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]
    extract_features(input_folder, output_folder)
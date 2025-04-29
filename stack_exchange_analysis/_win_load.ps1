# Dataset options: 'math', 'cs'
$env:DATASET = 'math'                           

# Folder options: 'math.stackexchange.com', 'stackoverflow.com'
$env:RAW_DATA_FOLDER = 'math.stackexchange.com' 

# Format options: 'csv', 'pkl' (csv is recommended)
$env:OUT_FORMAT = 'csv'                         

# Set to 0 to disable all filters
$env:FILTER_BY_VIEWCOUNT = 0
$env:MIN_ANSWERS = 0

# Activate the Python virtual environemnt
& .\.load\Scripts\Activate.ps1

# Extract, load and preprocess the dataset
python 7z_extraction.py
python load_dataset.py
python preprocessing.py

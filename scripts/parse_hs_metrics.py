import csv
import sys

# Function to parse and print each section
def parse_metrics_section(section):
    lines = section.strip().split('\n')
    if len(lines) < 2:
        return
    reader = csv.reader(lines, delimiter='\t')
    headers = next(reader)
    values = next(reader)
    for header, value in zip(headers, values):
        print(f"{header}\t{value}")
    print("\n" + "-"*40 + "\n")

# Check if a file name is provided
if len(sys.argv) < 2:
    print("Usage: python parse_hsmetrics_specific.py <filename>")
    sys.exit(1)

# Read the file and process each section
filename = sys.argv[1]
with open(filename, 'r') as file:
    data = file.read()

# Find and parse the section starting with '## METRICS CLASS picard.analysis.directed.HsMetrics'
start_marker = '## METRICS CLASS\tpicard.analysis.directed.HsMetrics'
start_index = data.find(start_marker)
if start_index != -1:
    # Extract the relevant part of the data
    metrics_section = data[start_index + len(start_marker):]
    # Find the end of the section (next '##' line or end of data)
    end_index = metrics_section.find('##')
    if end_index != -1:
        metrics_section = metrics_section[:end_index]
    # Parse and print the metrics section
    parse_metrics_section(metrics_section)
else:
    print("Metrics section not found in the file.")
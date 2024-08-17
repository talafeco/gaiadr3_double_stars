import csv
import sys

csv_file = sys.argv[1]

def remove_duplicates_from_csv(file_path):
    # Read the CSV file and store the rows in a list
    with open(file_path, mode='r', newline='') as file:
        reader = csv.reader(file)
        unique_rows = list(reader)
    
    # Remove duplicates by converting the list of rows to a set of tuples
    unique_rows_set = set(tuple(row) for row in unique_rows)
    
    # Convert the set back to a list of lists
    unique_rows = [list(row) for row in unique_rows_set]
    
    # Write the unique rows back to the CSV file
    with open(file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(unique_rows)
        
remove_duplicates_from_csv(csv_file)

# Example usage
# remove_duplicates_from_csv('path/to/your/file.csv')

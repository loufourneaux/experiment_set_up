import csv 

# Example of reading and filtering data
def read_coordinates(file_path, identifier, coord_type):
    """
    Read and filter coordinates from a CSV file based on an identifier and coordinate type.

    Args:
        file_path (str): Path to the CSV file containing the coordinates.
        identifier (str): The identifier to match in the first column of the CSV file.
        coord_type (str): The type of coordinates to match in the second column of the CSV file.

    Returns:
        List[List[float]]: A list of lists containing the filtered coordinates. Each inner list represents
                            a coordinate set, with coordinates as float values.
    """
    filtered_coords = []
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        for row in reader:
            if row[0] == identifier and row[1] == coord_type:
                filtered_coords.append([float(coord) for coord in row[2:]])
    return filtered_coords


    


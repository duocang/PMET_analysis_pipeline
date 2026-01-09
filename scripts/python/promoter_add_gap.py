import sys

def read_chromosome_ranges(chromosome_file):
    chromosome_ranges = {}
    with open(chromosome_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            chromosome_id = parts[0]
            start, end = int(parts[3]), int(parts[4])
            chromosome_ranges[chromosome_id] = (start, end)
    return chromosome_ranges

def adjust_promoter_coordinates(promoters_file, chromosome_ranges, distance):
    """
    Adjusts the coordinates of promoters based on their strand and a specified distance.
    It ensures that the adjusted coordinates do not exceed the limits of the corresponding chromosome.

    Args:
        promoters_file (str): The path to the file containing promoter information.
        chromosome_ranges (dict): A dictionary where each key is a chromosome ID and the value is a tuple
                                  (start, end) representing the range of that chromosome.
        distance (int): The distance to adjust the promoter coordinates by. For promoters on the '+' strand,
                        the distance is subtracted, and for those on the '-' strand, it is added.

    Returns:
        list of tuple: A list of tuples, where each tuple represents an adjusted promoter entry.
                       Each tuple contains the chromosome ID, adjusted start position, adjusted end position,
                       and other original fields from the promoter file.
    """
    adjusted_promoters = []
    with open(promoters_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            chromosome_id, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]

            # Check if the chromosome_id exists in chromosome_ranges
            if chromosome_id not in chromosome_ranges:
                # If not, skip this line and continue with the next
                continue

            # Get the chromosome range
            chrom_start, chrom_end = chromosome_ranges[chromosome_id]

            # Adjust coordinates based on strand
            if strand == '+':
                start, end = max(start - distance, chrom_start), max(end - distance, chrom_start)
            else:
                start, end = min(start + distance, chrom_end), min(end + distance, chrom_end)

            adjusted_promoters.append((chromosome_id, start, end, parts[3], parts[4], strand))

    return adjusted_promoters

def main():
    # Parse arguments
    distance = int(sys.argv[1]) if len(sys.argv) > 1 else 500
    chromosome_file = sys.argv[2]
    promoters_file = sys.argv[3]

    # Read chromosome ranges and adjust promoter coordinates
    chromosome_ranges = read_chromosome_ranges(chromosome_file)
    adjusted_promoters = adjust_promoter_coordinates(promoters_file, chromosome_ranges, distance)

    # Print adjusted promoters
    for promoter in adjusted_promoters:
        print('\t'.join(map(str, promoter)))

if __name__ == "__main__":
    main()

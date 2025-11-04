"""
Convert SeqGeq_demo.csv to h5mu format
"""
from converters.csv_converter import csv_to_h5mu
import os

csv_path = "/home/bdkoki/デスクトップ/csv_h5mu/SeqGeq_demo.csv"
output_path = "/home/bdkoki/デスクトップ/csv_h5mu/SeqGeq_demo.h5mu"

print("Starting conversion...")
print(f"Input: {csv_path}")
print(f"Output: {output_path}")
print(f"File size: {os.path.getsize(csv_path) / (1024*1024):.1f} MB")
print()

try:
    # For SeqGeq format:
    # - transpose=False (data is already cells x genes)
    # - has_header=True (column names are gene names)
    # - has_index=True (first column is cell IDs)
    csv_to_h5mu(
        csv_path=csv_path,
        output_path=output_path,
        transpose=False,
        has_header=True,
        has_index=True
    )

    print()
    print("Conversion completed successfully!")
    print(f"Output file: {output_path}")
    print(f"Output size: {os.path.getsize(output_path) / (1024*1024):.1f} MB")

except Exception as e:
    print(f"ERROR: {e}")
    import traceback
    traceback.print_exc()

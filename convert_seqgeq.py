"""
SeqGeq_demo.csv を h5mu 形式に変換するワンショットスクリプト。
（コマンドラインから直接実行する想定の検証用）
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
    # SeqGeq 形式向けのオプション:
    # - transpose=False ... 既に「細胞 × 遺伝子」の並びになっている
    # - has_header=True ... 列名は遺伝子名
    # - has_index=True ... 1 列目は細胞 ID
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

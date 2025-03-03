import pandas as pd
import re
import os
import csv

def check_delimiter(csv_file_path: str):
    with open(csv_file_path, 'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.readline())
        is_delimiter=dialect.delimiter
    return is_delimiter
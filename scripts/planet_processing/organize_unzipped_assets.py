# The code in this module was generated with the assistance of
# Anthropic's Claude Sonnet 5 model. Code was then reviewed by
# Cameron Scholl to ensure it met the requirements of the project.

import argparse
import os
import re
import shutil


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Sort unzipped Planet order imagery, where each date has '
                 'been extracted into its own subfolder, into two flat '
                 'folders: one for analytic surface reflectance images and '
                 'one for their matching UDM2 masks. Filenames are left '
                 'unchanged.')

parser.add_argument('--input', '-i', action='store', nargs='+', required=True,
                     help='One or more folders to search recursively for '
                          'analytic/UDM2 files.')

parser.add_argument('--analytic-output', '-a', action='store', required=True,
                     help='Destination folder for analytic SR images.')

parser.add_argument('--udm-output', '-u', action='store', required=True,
                     help='Destination folder for UDM2 masks.')

parser.add_argument('--analytic-pattern', action='store',
                     default=r'AnalyticMS_SR(_clip)?\.tif$',
                     help='Regex matched against each filename to identify '
                          'analytic SR images.')

parser.add_argument('--udm-pattern', action='store',
                     default=r'udm2(_clip)?\.tif$',
                     help='Regex matched against each filename to identify '
                          'UDM2 masks.')

parser.add_argument('--copy', action='store_true',
                     help='Copy files instead of moving them.')

parser.add_argument('--skip-existing', action='store_true',
                     help='Skip a file if one of the same name already '
                          'exists in its destination folder.')

args = parser.parse_args()

os.makedirs(args.analytic_output, exist_ok=True)
os.makedirs(args.udm_output, exist_ok=True)

analytic_re = re.compile(args.analytic_pattern)
udm_re = re.compile(args.udm_pattern)
transfer = shutil.copy2 if args.copy else shutil.move
verb = 'Copied' if args.copy else 'Moved'


##############################################################
# Find and sort matching files
##############################################################

moved = 0
skipped = 0
unmatched = 0

for folder in args.input:
    assert os.path.isdir(folder), f'{folder} is not a directory.'
    for root, _, names in os.walk(folder):
        for name in sorted(names):
            if analytic_re.search(name):
                dest_dir = args.analytic_output
            elif udm_re.search(name):
                dest_dir = args.udm_output
            else:
                unmatched += 1
                continue

            src_path = os.path.join(root, name)
            dest_path = os.path.join(dest_dir, name)

            if args.skip_existing and os.path.exists(dest_path):
                print(f'Skipped (already exists): {dest_path}')
                skipped += 1
                continue

            transfer(src_path, dest_path)
            print(f'{verb} {src_path} -> {dest_path}')
            moved += 1

print(f'{verb} {moved} file(s), skipped {skipped} existing, '
      f'{unmatched} file(s) did not match either pattern.')

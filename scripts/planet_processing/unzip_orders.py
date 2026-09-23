import argparse
import os
import zipfile


##############################################################
# Parse arguments
##############################################################

parser = argparse.ArgumentParser(
    description='Unzip Planet order archives from one or more source '
                 'folders into a single destination folder. Useful when '
                 'imagery for a single location is split across multiple '
                 'order downloads.')

parser.add_argument('--input', '-i', action='store', nargs='+', required=True,
                     help='One or more folders to search for zip files.')

parser.add_argument('--output', '-o', action='store', required=True,
                     help='Destination folder for the extracted contents.')

parser.add_argument('--skip-existing', action='store_true',
                     help='Skip a zip file if a file it would extract '
                          'already exists in the destination folder.')

args = parser.parse_args()

os.makedirs(args.output, exist_ok=True)


##############################################################
# Find zip files
##############################################################

zip_paths = []
for folder in args.input:
    assert os.path.isdir(folder), f'{folder} is not a directory.'
    for name in sorted(os.listdir(folder)):
        if name.lower().endswith('.zip'):
            zip_paths.append(os.path.join(folder, name))

assert zip_paths, f'No zip files found in {args.input}.'
print(f'Found {len(zip_paths)} zip file(s) across {len(args.input)} '
      'folder(s).')


##############################################################
# Extract zip files
##############################################################

extracted = 0
skipped = 0

for zip_path in zip_paths:
    with zipfile.ZipFile(zip_path) as archive:
        for member in archive.namelist():
            dest_path = os.path.join(args.output, member)
            if args.skip_existing and os.path.exists(dest_path):
                skipped += 1
                continue
            archive.extract(member, args.output)
            extracted += 1
    print(f'Extracted {zip_path}')

print(f'Extracted {extracted} file(s), skipped {skipped} existing file(s).')

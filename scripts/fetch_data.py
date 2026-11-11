"""Download the public CDC CSV and verify the exact analyzed input.

Only CSV serialization is normalized: nonempty cells are quoted, empty cells
are unquoted, and lines use LF. Cell values and row order remain unchanged.
"""
from pathlib import Path
import csv
import hashlib
import io
import json
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
URL = 'https://data.cdc.gov/api/views/37nu-tuw8/rows.csv?accessType=DOWNLOAD'
EXPECTED = '2DFCE56408F4411C01426556C250DBD932790CFD24293B278E63F32A8323EDED'

def canonical_bytes(payload):
    rows = list(csv.reader(io.StringIO(payload.decode('utf-8-sig'), newline='')))
    if len(rows) != 580 or any(len(row) != 40 for row in rows):
        raise ValueError('Expected one header and 579 data rows with 40 columns.')
    quote = '"'
    serialized = ''.join(','.join(quote + cell.replace(quote, quote + quote) + quote
                                  if cell else '' for cell in row) + '\n' for row in rows).encode('utf-8')
    if list(csv.reader(io.StringIO(serialized.decode(), newline=''))) != rows:
        raise ValueError('Serialization changed parsed cell values.')
    if hashlib.sha256(serialized).hexdigest().upper() != EXPECTED:
        raise ValueError('CDC content differs from the analyzed version. No analyzed input was written.')
    return serialized

def main():
    destination = ROOT / 'input' / 'data.csv'
    if destination.exists():
        if hashlib.sha256(destination.read_bytes()).hexdigest().upper() == EXPECTED:
            print('The verified analyzed input already exists.')
            return
        raise FileExistsError('Existing input differs. Preserve and inspect it before downloading.')
    with urllib.request.urlopen(URL, timeout=90) as response:
        payload = response.read()
    analyzed = canonical_bytes(payload)
    destination.parent.mkdir(parents=True, exist_ok=True)
    (destination.parent / 'cdc_download.csv').write_bytes(payload)
    destination.write_bytes(analyzed)
    provenance = {'source': URL, 'source_sha256': hashlib.sha256(payload).hexdigest().upper(),
                  'analyzed_sha256': EXPECTED, 'rows': 579, 'columns': 40,
                  'parsed_cell_values_unchanged': True}
    (destination.parent / 'download_provenance.json').write_text(json.dumps(provenance, indent=2)+'\n')
    print('Downloaded 579 records; analyzed input SHA-256 verified.')

if __name__ == '__main__':
    main()

"""Real filesystem checks for case collisions, payload identity and legacy caches."""
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
import numpy as np
from molecular_asset_io import chain_file, chain_filenames


class ChainFiles(unittest.TestCase):
    def test_distinct_payloads_and_legacy_guard(self):
        names = ['LA', 'La', 'A', '../unsafe', '链']
        mapping = chain_filenames(names)
        self.assertEqual(len({s.casefold() for s in mapping.values()}), len(names))
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder); rows = []
            for i, name in enumerate(names):
                path = root/mapping[name]
                np.savez_compressed(path, source_chain=np.array(name), atoms=np.full((i+1, 3), i),
                                    vertices=np.ones((3, 3)), faces=np.array([[0, 1, 2]]))
                rows.append({'chain': name, 'file': path.name, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                             'atoms': i+1, 'vertices': 3, 'faces': 1})
            manifest = {'chain_file_schema': 1, 'chains': rows}
            target = root/'manifest.json'; target.write_text(json.dumps(manifest))
            for name in names:
                with np.load(chain_file(root, name)) as data:
                    self.assertEqual(data['source_chain'].item(), name)
            self.assertEqual(len(list(root.glob('*.npz'))), len(names))
            # Same-size wrong identity still fails even if a new hash is recorded.
            path = root/mapping['LA']; data = dict(np.load(path)); data['source_chain'] = np.array('La')
            np.savez_compressed(path, **data)
            with self.assertRaisesRegex(ValueError, 'hash differs'): chain_file(root, 'LA')
            rows[0]['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest(); target.write_text(json.dumps(manifest))
            with self.assertRaisesRegex(ValueError, 'payload identity'): chain_file(root, 'LA')
            target.write_text(json.dumps({'chains': [{'chain': 'LA', 'atoms': 1}, {'chain': 'La', 'atoms': 2}]}))
            with self.assertRaisesRegex(ValueError, 'Ambiguous chain filenames'): chain_file(root, 'LA')

    def test_unambiguous_legacy(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            np.savez_compressed(root/'chain_A.npz', atoms=np.ones((2, 3)), vertices=np.ones((3, 3)), faces=np.array([[0, 1, 2]]))
            row = {'chain': 'A', 'atoms': 2, 'vertices': 3, 'faces': 1}
            (root/'manifest.json').write_text(json.dumps({'chains': [row]}))
            self.assertEqual(chain_file(root, 'A').name, 'chain_A.npz')
            row['atoms'] = 3; (root/'manifest.json').write_text(json.dumps({'chains': [row]}))
            with self.assertRaisesRegex(ValueError, 'inventory differs'): chain_file(root, 'A')


if __name__ == '__main__': unittest.main()

import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from public_provenance import source_identity_sha256

class PublicProvenanceTests(unittest.TestCase):
    def make(self,root):
        rel='data/rpe1/preparation/manifest.json';p=root/rel;p.parent.mkdir(parents=True);p.write_text('{"path":"reconstruction/raw/source.h5ad"}')
        r=dict(public_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),original_sha256='a'*64,source_identity_redaction=True)
        q=root/'provenance/public_artifact_mapping.json';q.parent.mkdir();q.write_text(json.dumps({'changed_files':{rel:r}}))
        return p
    def test_documented_original_identity(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=self.make(root)
            self.assertEqual(source_identity_sha256(root,p),'a'*64)
    def test_public_byte_change_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=self.make(root);p.write_text('{}')
            with self.assertRaises(AssertionError):source_identity_sha256(root,p)
    def test_unlisted_file_uses_actual_digest(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/'raw.txt';p.write_text('real bytes')
            self.assertEqual(source_identity_sha256(root,p),hashlib.sha256(p.read_bytes()).hexdigest())

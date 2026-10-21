"""Resolve explicitly redacted metadata identities without changing byte hashes.

The public package manifest checks actual bytes. Frozen study manifests contain
original metadata identities; only the two declared path-redacted manifest views
may resolve to a documented original identity after their public bytes validate.
"""
import hashlib
import json
from pathlib import Path

REDACTED_MANIFESTS = frozenset({
    'data/rpe1/preparation/manifest.json',
    'provenance/rpe1_contents_v1/manifest.json',
})

def byte_sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def source_identity_sha256(root, path):
    root, path = Path(root).resolve(), Path(path).resolve()
    rel = path.relative_to(root).as_posix()
    actual = byte_sha256(path)
    if rel not in REDACTED_MANIFESTS:
        return actual
    mapping = json.loads((root/'provenance/public_artifact_mapping.json').read_text())
    record = mapping['changed_files'][rel]
    if not record['source_identity_redaction'] or actual != record['public_sha256']:
        raise AssertionError('Redacted public metadata digest mismatch: '+rel)
    original = record['original_sha256']
    if len(original) != 64 or any(c not in '0123456789abcdef' for c in original):
        raise AssertionError('Invalid declared original metadata identity')
    return original

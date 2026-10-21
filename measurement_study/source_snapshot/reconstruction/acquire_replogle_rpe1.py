"""Acquire the author-deposited Replogle RPE1 raw single-cell H5AD.

The download is resumed when the Figshare endpoint honors HTTP Range requests.
The completed file is accepted only when both declared size and Figshare MD5
match the live article metadata.
"""
from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw"
PROVENANCE = HERE / "provenance"
SNAPSHOTS = PROVENANCE / "source_snapshots"
RECEIPTS = PROVENANCE / "download_receipts"
ARTICLE_ID = 20029387
FILE_NAME = "rpe1_raw_singlecell_01.h5ad"
FILE_ID = 35775606
API_URL = f"https://api.figshare.com/v2/articles/{ARTICLE_ID}"
EXPECTED_STUDY_MD5 = "6a2a9d0d2bf4ec147f4d1104043b268c"


def fetch_json(url: str) -> dict:
    with urlopen(Request(url, headers={"User-Agent": "CARF-provenance-audit/1.0"}), timeout=60) as response:
        return json.load(response)


def file_hashes(path: Path) -> dict[str, str]:
    md5 = hashlib.md5()
    sha256 = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            md5.update(block)
            sha256.update(block)
    return {"md5": md5.hexdigest(), "sha256": sha256.hexdigest()}


def download_resumable(url: str, path: Path, expected_size: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    offset = path.stat().st_size if path.exists() else 0
    if offset > expected_size:
        path.unlink()
        offset = 0
    headers = {"User-Agent": "CARF-provenance-audit/1.0"}
    if offset:
        headers["Range"] = f"bytes={offset}-"
    request = Request(url, headers=headers)
    try:
        response = urlopen(request, timeout=180)
    except (HTTPError, URLError, TimeoutError) as exc:
        raise RuntimeError(f"Could not open Figshare download: {exc}") from exc
    status = getattr(response, "status", response.getcode())
    mode = "ab" if offset and status == 206 else "wb"
    if mode == "ab":
        content_range = response.headers.get("Content-Range", "")
        if not content_range.startswith(f"bytes {offset}-"):
            response.close()
            raise RuntimeError(f"Unexpected range response: {content_range!r}")
    elif offset and status == 200:
        # The server ignored Range; restart cleanly instead of appending a full copy.
        offset = 0
    elif status not in (200, 206):
        response.close()
        raise RuntimeError(f"Unexpected HTTP status: {status}")
    written = offset
    with response, path.open(mode) as output:
        while True:
            block = response.read(8 * 1024 * 1024)
            if not block:
                break
            output.write(block)
            written += len(block)
            if written % (256 * 1024 * 1024) < len(block):
                output.flush()
                print(f"Downloaded {written:,} / {expected_size:,} bytes", flush=True)
    if written != expected_size:
        raise RuntimeError(f"Incomplete Figshare download: {written} of {expected_size} bytes")


def main() -> None:
    SNAPSHOTS.mkdir(parents=True, exist_ok=True)
    RECEIPTS.mkdir(parents=True, exist_ok=True)
    article = fetch_json(API_URL)
    snapshot_path = SNAPSHOTS / f"figshare_article_{ARTICLE_ID}.json"
    snapshot_path.write_text(json.dumps(article, ensure_ascii=False, indent=2) + "\n")
    record_file = next((f for f in article["files"] if f["id"] == FILE_ID), None)
    if record_file is None or record_file["name"] != FILE_NAME:
        raise RuntimeError("Expected RPE1 file was not present in the live Figshare record")
    expected_size = int(record_file["size"])
    expected_md5 = str(record_file.get("supplied_md5", "")).lower()
    if expected_size != 8_700_873_216 or expected_md5 != EXPECTED_STUDY_MD5:
        raise RuntimeError("Live Figshare size/checksum metadata differ from the pinned source")

    path = RAW / FILE_NAME
    if not path.exists() or path.stat().st_size != expected_size:
        download_resumable(record_file["download_url"], path, expected_size)
    hashes = file_hashes(path)
    if hashes["md5"] != expected_md5:
        raise RuntimeError(f"Figshare MD5 mismatch: {hashes['md5']} != {expected_md5}")

    receipt = {
        "status": "verified_author_deposit",
        "checked_utc": datetime.now(timezone.utc).isoformat(),
        "source_record": {
            "title": article.get("title"),
            "doi": article.get("doi"),
            "record_url": f"https://doi.org/{article.get('doi')}",
            "record_api": API_URL,
            "record_snapshot": str(snapshot_path.relative_to(HERE)),
            "record_snapshot_sha256": file_hashes(snapshot_path)["sha256"],
            "file_id": FILE_ID,
            "file_name": FILE_NAME,
            "download_url": record_file["download_url"],
            "bytes": expected_size,
            "expected_md5": expected_md5,
        },
        "study": {
            "citation": "Replogle et al. 2022, Cell, doi:10.1016/j.cell.2022.05.013",
            "cell_context": "RPE1 essential-scale CRISPRi Perturb-seq, day 7",
            "primary_raw_archive": "NCBI SRA BioProject PRJNA831566",
            "file_scope": "author-deposited raw single-cell H5AD; not a reprocessed reconstruction from FASTQ",
            "license": article.get("license"),
        },
        "local_file": str(path.relative_to(HERE)),
        "actual_bytes": path.stat().st_size,
        "actual_md5": hashes["md5"],
        "actual_sha256": hashes["sha256"],
        "validation": "live Figshare file ID, declared size, supplied MD5, downloaded file size, and downloaded MD5 all agree",
    }
    (RECEIPTS / "Replogle_RPE1_raw_singlecell.json").write_text(
        json.dumps(receipt, ensure_ascii=False, indent=2) + "\n"
    )
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()

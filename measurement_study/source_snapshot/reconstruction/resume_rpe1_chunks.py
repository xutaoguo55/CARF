"""Resume a pinned public deposit without curl retries truncating the prefix.

Each bounded range uses a fresh public redirect and is written to a sidecar.
Only complete, range-checked chunks are appended in byte order. Final trust
still requires the author's full-file MD5 through acquire_replogle_rpe1.
"""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from urllib.request import Request, urlopen
import time
import json

from acquire_replogle_rpe1 import (
    RAW, FILE_NAME, FILE_ID, API_URL, EXPECTED_STUDY_MD5, fetch_json, main as verify,
)

SIZE = 8_700_873_216
CHUNK = 64 * 1024 * 1024
WORKERS = 4


def fetch_range(start, end, folder):
    target = folder / f"{start}-{end}.part"
    for attempt in range(4):
        try:
            url = f"https://ndownloader.figshare.com/files/{FILE_ID}?refresh={time.time_ns()}"
            req = Request(url, headers={"Range": f"bytes={start}-{end}",
                                       "User-Agent": "CARF-source-verification/1.0"})
            with urlopen(req, timeout=60) as response:
                expected_range = f"bytes {start}-{end}/{SIZE}"
                if response.status != 206 or response.headers.get("Content-Range") != expected_range:
                    raise ValueError("Server did not return the exact pinned byte range")
                written = 0
                with target.open("wb") as out:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        written += len(block)
                        if written > end-start+1:
                            raise ValueError("Oversized range response")
                        out.write(block)
                if written != end-start+1:
                    raise ValueError("Incomplete range response")
            return start, end, target
        except Exception as exc:
            print(f"Range {start}-{end} attempt {attempt+1}: {type(exc).__name__}", flush=True)
            if attempt == 3:
                raise
            time.sleep(2)


def main():
    article = fetch_json(API_URL)
    record = next(f for f in article["files"] if f["id"] == FILE_ID)
    if record["name"] != FILE_NAME or record["size"] != SIZE or record["supplied_md5"] != EXPECTED_STUDY_MD5:
        raise ValueError("Live source metadata differ from pinned deposit")
    path = RAW / FILE_NAME
    if not path.exists() or not 0 < path.stat().st_size <= SIZE:
        raise ValueError("Expected an existing partial file of a valid length")
    folder = RAW / "rpe1_resume_chunks"
    folder.mkdir(exist_ok=True)
    start = path.stat().st_size
    print(json.dumps({"resume_offset": start, "expected_bytes": SIZE,
                      "chunk_bytes": CHUNK, "parallel_ranges": WORKERS}), flush=True)
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        while start < SIZE:
            ranges = [(a, min(a+CHUNK,SIZE)-1)
                      for a in range(start, min(start+WORKERS*CHUNK,SIZE), CHUNK)]
            jobs = [pool.submit(fetch_range,a,b,folder) for a,b in ranges]
            for job in jobs:
                a,b,chunk = job.result()
                if path.stat().st_size != a:
                    raise ValueError("Partial file changed concurrently; refusing append")
                with path.open("ab") as out, chunk.open("rb") as inp:
                    for block in iter(lambda: inp.read(4*1024*1024),b""):
                        out.write(block)
                if path.stat().st_size != b+1:
                    raise ValueError("Appended size mismatch")
                chunk.unlink()
                start=b+1
            print(f"Range-checked bytes {start}/{SIZE}; full-file checksum still pending",flush=True)
    verify()


if __name__ == "__main__":
    main()

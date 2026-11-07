#!/usr/bin/env python3
"""Archive a release of this repository on Zenodo and publish it (mints a DOI).

Called by the `zenodo` job in .gitlab-ci.yml on every git tag. Uses only the
Python standard library.

Environment variables:
    ZENODO_TOKEN          Personal access token (scopes: deposit:write, deposit:actions).
    ZENODO_DEPOSITION_ID  Optional. ID of any already published version of this
                          software on Zenodo. If set, a new version of that record
                          is created (same concept DOI). If unset, a brand new
                          record is created; set the variable afterwards.
    ZENODO_URL            Optional. Defaults to https://zenodo.org. Use
                          https://sandbox.zenodo.org for testing.

Usage:
    python ci/zenodo_release.py <archive.zip> <version>
"""

import datetime
import json
import os
import sys
import urllib.error
import urllib.parse
import urllib.request

METADATA_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".zenodo.json")


def request(method, url, token, data=None, content_type="application/json"):
    headers = {"Authorization": f"Bearer {token}"}
    if data is not None:
        if content_type == "application/json":
            data = json.dumps(data).encode()
        headers["Content-Type"] = content_type
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req) as resp:
            body = resp.read()
    except urllib.error.HTTPError as e:
        sys.exit(f"Zenodo API error: {method} {url} -> {e.code}\n{e.read().decode(errors='replace')}")
    return json.loads(body) if body else None


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    archive, version = sys.argv[1], sys.argv[2]

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        sys.exit("ZENODO_TOKEN is not set.")
    api = os.environ.get("ZENODO_URL", "https://zenodo.org").rstrip("/") + "/api"
    previous_id = os.environ.get("ZENODO_DEPOSITION_ID", "").strip()

    with open(METADATA_FILE) as f:
        metadata = json.load(f)
    metadata["version"] = version
    metadata["publication_date"] = datetime.date.today().isoformat()

    if previous_id:
        print(f"Creating new version of Zenodo record {previous_id}")
        resp = request("POST", f"{api}/deposit/depositions/{previous_id}/actions/newversion", token)
        draft = request("GET", resp["links"]["latest_draft"], token)
        # A new version starts with the files of the previous one; replace them.
        for f in draft.get("files", []):
            request("DELETE", f"{draft['links']['self']}/files/{f['id']}", token)
    else:
        print("ZENODO_DEPOSITION_ID not set: creating a new Zenodo record")
        draft = request("POST", f"{api}/deposit/depositions", token, data={})

    filename = os.path.basename(archive)
    print(f"Uploading {filename}")
    with open(archive, "rb") as f:
        request(
            "PUT",
            f"{draft['links']['bucket']}/{urllib.parse.quote(filename)}",
            token,
            data=f.read(),
            content_type="application/octet-stream",
        )

    request("PUT", draft["links"]["self"], token, data={"metadata": metadata})
    record = request("POST", draft["links"]["publish"], token)

    print(f"Published record:  {record['links'].get('record_html', record['id'])}")
    print(f"Version DOI:       {record.get('doi')}")
    print(f"Concept DOI:       {record.get('conceptdoi')}  (cite this one for all versions)")
    if not previous_id:
        print(f"\nSet the CI/CD variable ZENODO_DEPOSITION_ID={record['id']} so future tags become new versions.")


if __name__ == "__main__":
    main()

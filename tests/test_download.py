"""Tests for the retry around RCSB downloads. No network access is used."""

from __future__ import annotations

import json
import urllib.error

import pytest

from plinter import download
from plinter.chemistry import fetch_chemcomp
from plinter.structures import download_structure


class _Response:
    def __init__(self, body: bytes) -> None:
        self._body = body

    def __enter__(self) -> _Response:
        return self

    def __exit__(self, *exc_info: object) -> None:
        return None

    def read(self) -> bytes:
        return self._body


def _http_error(code: int) -> urllib.error.HTTPError:
    return urllib.error.HTTPError("https://files.rcsb.org/x", code, "error", {}, None)


RESET = urllib.error.URLError(ConnectionResetError(104, "Connection reset by peer"))
TIMEOUT = urllib.error.URLError(TimeoutError("timed out"))


@pytest.fixture
def network(monkeypatch):
    """Replace urlopen with a scripted sequence of failures and responses."""
    calls: list[str] = []
    waits: list[float] = []
    script: list[object] = []

    def fake_urlopen(url, timeout):
        calls.append(url)
        # Once the script runs out, every further attempt is reset.
        outcome = script.pop(0) if script else RESET
        if isinstance(outcome, BaseException):
            raise outcome
        return _Response(outcome)

    monkeypatch.setattr(download.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(download, "_sleep", waits.append)
    return script, calls, waits


@pytest.mark.parametrize(
    "failure",
    [RESET, TIMEOUT, ConnectionResetError(104, "reset"), _http_error(503)],
    ids=["connection-reset", "timeout", "bare-reset", "http-503"],
)
def test_a_transient_failure_is_retried_and_then_succeeds(network, tmp_path, failure):
    script, calls, waits = network
    script.extend([failure, b"HEADER    1KMV\n"])

    path = download_structure("1kmv", tmp_path / "pdb" / "1kmv.pdb")

    assert path.read_bytes() == b"HEADER    1KMV\n"
    assert calls == ["https://files.rcsb.org/download/1KMV.pdb"] * 2
    assert waits == [download.BACKOFF_SECONDS]


def test_a_persistent_failure_raises_an_actionable_error(network, tmp_path):
    _, calls, waits = network
    destination = tmp_path / "pdb" / "1kmv.pdb"

    with pytest.raises(download.DownloadError) as caught:
        download_structure("1KMV", destination)

    message = str(caught.value)
    assert f"after {download.ATTEMPTS} attempts" in message
    assert "Connection reset by peer" in message
    assert str(destination) in message, "the error should say where to put the file"
    assert len(calls) == download.ATTEMPTS
    # Exponential backoff between attempts, none after the last one.
    assert waits == [download.BACKOFF_SECONDS * 2**i for i in range(download.ATTEMPTS - 1)]
    assert not destination.exists(), "a failed download must not leave a partial file"


def test_a_client_error_is_not_retried(network, tmp_path):
    script, calls, waits = network
    script.append(_http_error(404))

    with pytest.raises(urllib.error.HTTPError):
        download_structure("0XXX", tmp_path / "0xxx.pdb")

    assert len(calls) == 1
    assert waits == []


def test_chemcomp_records_use_the_same_retry(network, tmp_path):
    script, calls, _ = network
    record = {"rcsb_id": "LII"}
    script.extend([TIMEOUT, json.dumps(record).encode()])

    assert fetch_chemcomp("lii", tmp_path / "chemcomp" / "LII.json") == record
    assert len(calls) == 2

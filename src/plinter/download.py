"""HTTP retrieval with a bounded retry for transient network failures.

RCSB occasionally resets a connection or times out. Without a retry, one such
blip fails a whole run even though the next request would succeed. Only
failures that a later attempt can plausibly fix are retried (connection resets,
timeouts, DNS hiccups, HTTP 5xx and 429). A 404 or another client error is
raised at once, because asking again cannot change the answer.
"""

from __future__ import annotations

import http.client
import time
import urllib.error
import urllib.request

ATTEMPTS = 3
BACKOFF_SECONDS = 2.0
TIMEOUT_SECONDS = 60.0

# Looked up at call time so tests can replace it and run without waiting.
_sleep = time.sleep


class DownloadError(RuntimeError):
    """A download still failed after every retry."""


def _is_transient(error: BaseException) -> bool:
    if isinstance(error, urllib.error.HTTPError):
        return error.code >= 500 or error.code == 429
    if isinstance(error, urllib.error.URLError):
        # URLError wraps the socket-level cause (reset, timeout, DNS failure)
        # in .reason. A plain string reason, such as an unknown URL scheme, is
        # a usage error rather than a network one.
        return isinstance(error.reason, OSError)
    return isinstance(error, (OSError, http.client.HTTPException))


def fetch(
    url: str,
    *,
    attempts: int = ATTEMPTS,
    backoff: float = BACKOFF_SECONDS,
    timeout: float = TIMEOUT_SECONDS,
    hint: str = "",
) -> bytes:
    """Return the body of ``url``, retrying transient failures with backoff.

    The wait doubles after each failed attempt (2 s, then 4 s by default).

    Raises:
        DownloadError: if every attempt failed with a transient error. The
            message names the URL, the attempt count, the last error and
            ``hint``, which callers use to say how to work around it.
        urllib.error.HTTPError: at once, for a client error such as a 404.
    """
    last_error: BaseException | None = None
    for attempt in range(1, attempts + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as response:
                body: bytes = response.read()
                return body
        except Exception as error:
            if not _is_transient(error):
                raise
            last_error = error
            if attempt < attempts:
                _sleep(backoff * 2 ** (attempt - 1))

    message = f"Could not download {url} after {attempts} attempts (last error: {last_error})."
    if hint:
        message = f"{message} {hint}"
    raise DownloadError(message) from last_error

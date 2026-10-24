# HTTP helpers and demonstration server

This package supplies utilities used by BBTools services. The actual taxonomy,
clade, SSU and demultiplexing services live in their respective packages.
`SimpleHttpServer` is a demonstration with fixed mock data, not one of those
production services.

## Classes

| Class | Purpose |
|---|---|
| `ServerTools` | HTTP request/response helpers, stream reading, FTP helpers, gzip and JSON utilities |
| `PercentEncoding` | Targeted URL escaping and percent decoding used by taxonomy queries and resource paths |
| `SimpleHttpServer` | Small demonstration listener on port 8321; returns a requested key/value and fixed `tax_id` and `organism` values |

`SimpleHttpServer` preserves raw URI text, including percent escapes and query
text. It retains the historical status-200 error object for paths that do not
contain two components. Its output is UTF-8 JSON, and it closes each exchange
when handling completes or fails. It does not perform biological lookups.

## Compatibility and failure behavior

`ServerTools` currently imports Java 11's `java.net.http` API. Although it has a
legacy `HttpURLConnection` implementation and a Java-version dispatch, the whole
package does not compile against the Java 8 API. The legacy path alone is not a
Java 8 compatibility guarantee.

The request helpers have different retry and timeout behavior. Their return
values and resource ownership must be checked at the individual method; an empty
or null body is not proof that a request succeeded. HTTP status codes, I/O
failures and interruption are distinct conditions. Assertions are not a substitute
for checking these outcomes when a service runs with `-da`.

## Tests and load clients

Regression fixtures and the three load clients live in `test/server`, outside
the production package. `StressTest` uses bounded asynchronous requests,
`StressTest2` uses a bounded worker pool, and `StressTest3` sends requests serially.
Each requires an explicit endpoint, waits for its work, reports failed responses,
and fails independently of JVM assertions. The clients require Java 11 or later.
Their connect timeout is 10 seconds and request timeout is 30 seconds.

The mock fixtures open no listening sockets. `ServerLoopbackTest` requires an
explicit `loopback` argument, binds only to `127.0.0.1` on an operating-system
assigned port, uses two handler threads, and enforces a 60-second lifetime.
It exercises small HTTP success/failure cases and shuts down afterward. These
finite tests do not certify every distributed deployment, proxy configuration,
network interruption, or production workload.

Never point a load test at a service without the operator's authorization.

# Server Package - HTTP Server and Network Utilities

BBTools server package providing HTTP server implementation, network utilities, and stress testing tools for web service development and performance evaluation.

## HTTP Server Implementation

### SimpleHttpServer (SimpleHttpServer.java)
**Purpose**: A lightweight HTTP server for handling RESTful-style parameter requests and generating JSON responses.
- **Core Function**: Creates a minimal HTTP server on port 8321 that parses URL path components into key-value pairs and returns them as JSON
- **Key Features**:
  - Creates HttpServer on predefined port (8321)
  - Supports RESTful-style parameter parsing (e.g., `gi/123456`)
  - Generates dynamic JSON responses based on URL components
  - Hardcodes some default parameters like tax_id and organism
  - Provides basic error handling for incorrect URL formats
- **Usage**: Used for lightweight parameter retrieval and JSON response generation, particularly in bioinformatics contexts like retrieving information about genetic identifiers
- **Implementation Notes**:
  - Manual JSON construction instead of using JSON libraries
  - Potential bug in parameter parsing (line 42 might not correctly handle trailing slashes)
  - Hardcoded values suggest a prototype or specific use case
  - Debug output of processed parameters to console

## Network Utilities

### ServerTools (ServerTools.java)
**Purpose**: Utility class for network communication, HTTP/FTP interactions, and stream processing
- **Core Function**: Provides methods for sending/receiving network requests, handling HTTP exchanges, and managing network streams
- **Key Features**:
  - Supports both HttpURLConnection and modern HttpClient APIs
  - Robust error handling and retry mechanisms for network operations
  - Stream reading with dynamic buffer resizing
  - Flexible response handling (byte arrays, gzip compression)
  - Client address and forwarding IP detection
- **Usage**: Network communication, server-side interactions, HTTP/FTP request processing in BBTools ecosystem

### PercentEncoding (PercentEncoding.java)
**Purpose**: Provides utilities for encoding and decoding URL-safe percent-encoded strings.
- **Core Function**: Converts special and common symbols to and from their percent-encoded representations
- **Key Features**:
  - Convert special symbols to percent-encoded equivalents
  - Decode percent-encoded strings back to original symbols
  - Supports both reserved and common symbol encoding
  - Efficient BitSet-based symbol detection
- **Usage**: URL parameter encoding, safe string transmission across network protocols, web request/response handling

## Performance Testing

### StressTest (StressTest.java)
**Purpose**: Performs high-volume asynchronous HTTP request stress testing against a taxonomy service endpoint.
- **Core Function**: Generates and dispatches a configurable number of HTTP requests to a specified URL without waiting for responses
- **Key Features**:
  - Configurable iteration count (default 10,000 requests)
  - Default payload targeting taxonomy accession
  - Uses HTTP/2 for request efficiency
  - Asynchronous request dispatching without blocking
  - Precise timing measurement of request generation
- **Usage**: Performance testing and benchmarking network request capabilities, simulating high-load scenarios for taxonomic service endpoints

### StressTest2 (StressTest2.java)
**Purpose**: Concurrent HTTP request stress testing tool for URL endpoint performance evaluation.
- **Core Function**: Generates and sends a configurable number of concurrent HTTP requests to a specified URL
- **Key Features**:
  - Dynamically adjusts concurrency based on available processors
  - Supports configurable iteration count via command-line arguments
  - Uses Java's ExecutorService for managing concurrent request threads
  - Utilizes HTTP/2 client for modern network communication
  - Captures and times total request execution
- **Usage**: Performance testing network endpoints, particularly taxonomy service queries
  - Allows custom payload and iteration count specification
  - Prints total query time and execution statistics
  - Useful for assessing server responsiveness and load handling capabilities

### StressTest3 (StressTest3.java)
**Purpose**: Performs synchronous HTTP request stress testing against a taxonomy service endpoint.
- **Core Function**: Generates and sends a configurable number of synchronous HTTP requests to a specific URL, waiting for each response
- **Key Features**:
  - Configurable iteration count (default 10,000 requests)
  - Default payload targeting JGI taxonomy accession service
  - Uses HTTP/2 for request communication
  - Synchronous request handling with response validation
  - Precise timing measurement of total query execution
  - Enables SSL debug handshake logging
  - Basic error handling for failed requests
- **Usage**: Performance testing and benchmarking network request capabilities with synchronous request processing, specifically for testing taxonomy service endpoints

## Package Overview

The server package provides essential infrastructure for:

1. **HTTP Service Development**: SimpleHttpServer offers a foundation for creating lightweight web services with RESTful parameter handling
2. **Network Communication**: ServerTools provides robust utilities for HTTP/FTP interactions with comprehensive error handling
3. **URL Processing**: PercentEncoding ensures safe transmission of special characters in network protocols
4. **Performance Validation**: Three different stress testing approaches (asynchronous, concurrent, synchronous) enable comprehensive load testing

The package is particularly suited for bioinformatics applications requiring web service integration, with specialized support for taxonomy services and genetic identifier processing.
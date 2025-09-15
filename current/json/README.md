# json Package - JSON Processing Utilities
*Lightweight JSON parsing and generation tools for BBTools*

## Package Overview
The `json` package provides custom JSON processing capabilities optimized for BBTools output formatting and data exchange.

---

## Core Components

### JsonParser (JsonParser.java)
**Purpose**: Thread-safe JSON parsing library for BBTools supporting nested structures
**Core Function**: Converts JSON strings/byte arrays into Java objects, handling complex nested objects and arrays with state machine parsing
**Key Features**: 
  - Supports parsing objects, arrays with nested types
  - Handles escape sequences and quote modes
  - Flexible boolean and numeric type parsing
**Usage**: `JsonObject obj = new JsonParser(jsonString).parseJsonObject()`
**Dependencies**: java.util.ArrayList, structures.ByteBuilder

### JsonObject (JsonObject.java)
**Purpose**: Custom JSON object implementation for BBTools output formatting
**Core Function**: Builds flexible JSON structures with preserved key insertion order, supporting nested objects, arrays, and literal values with precise numeric formatting
**Key Features**: 
- Supports nested objects and dynamic key-value pair addition
- Handles various data types with custom formatting
- Provides detailed JSON string conversion methods
**Usage**: `JsonObject bob = new JsonObject("name", "bob"); bob.add("age", 30, true);`
**Dependencies**: java.util.LinkedHashMap, structures.ByteBuilder

### JsonLiteral (JsonLiteral.java)
**Purpose**: Represents unquoted literal values in JSON output for precise numeric formatting
**Core Function**: Enables creation of JSON-compatible numeric literals that can be precisely formatted without string quotation, supporting both pre-formatted strings and numeric values with specific decimal precision
**Key Features**: 
- Supports direct string literals
- Formats double values with specified decimal places
- Overrides toString() for direct JSON output
**Usage**: `JsonLiteral literal = new JsonLiteral(3.14159, 2); // Creates "3.14"`
**Dependencies**: shared.Tools

---

## Package Usage
The json package enables BBTools to generate structured JSON output for results, statistics, and configuration data with precise formatting control.

---
*Documentation generated using evidence-based analysis of source code*
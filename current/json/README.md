# JSON values and text

`JsonObject` builds JSON objects, `JsonParser` reads UTF-8 objects and arrays,
and `JsonLiteral` stores trusted, unquoted fragments. The package has no external
JSON-library dependency. Use a separate mutable parser/object per thread, or
coordinate access externally.

## Objects and serialization

`add` replaces a member even when its value changes between a scalar and a nested
object. `addAndRename` finds an unused name across both maps, beginning with the
suffix ` 2`. Values and children are borrowed rather than copied. Null values
serialize as JSON `null`; nested structures must be acyclic. Numeric Float/Double
NaN and either infinity also serialize as `null`, including array elements and
formatted numeric literals. JSON has no numeric NaN or infinity tokens.

Storage remains two public insertion-ordered maps: `omap` for scalar/array values
and `jmap` for non-null nested objects. Serialization emits the first map and then
the second, preserving insertion order within each. Direct map edits must preserve
those invariants and avoid duplicate names across maps. A Java null key retains
the historical spelling `"null"`; use ordinary non-null String names for reliable
interchange.

Quoted strings escape control characters, quotes, backslashes and non-ASCII UTF-16
units. Output is ASCII, which is also valid UTF-8; escaped surrogate pairs retain
supplementary Unicode characters. Primitive arrays and collections are serialized
element by element. `toText` appends to a supplied ByteBuilder without clearing it.

`toString(Object[])` includes array brackets. The legacy
`toString(ArrayList<JsonObject>)` helper emits only comma/newline-separated objects,
without surrounding brackets. It is not a complete JSON array by itself.

## Parsing and validation

String inputs are encoded explicitly as UTF-8. Byte-array inputs are borrowed and
must remain unchanged while parsing. Parsed integers use `Long`; fractional or
exponent values use finite `Double`, with the corresponding range and precision
limits. Typed getters do not generally coerce arbitrary Number implementations;
check their method contracts.

Parsing methods retain historical token extensions such as case-insensitive
booleans/null, `t`/`f`, and Java-compatible numeric forms. They return partial values
on malformed input and set `errorState()`. Check that flag before using data from
an external source. `set(byte[])` resets position and errors for reuse.

`validate()` checks a complete value from the beginning, including separators and
trailing text, without consuming the saved position. Empty input is not a complete
value. It requires standard lowercase literals, JSON number syntax and recognized
string escapes, while parsing retains the extensions above. The same Long/finite
Double representation limits apply. Validation updates the error flag; use `set`
to clear it and restart parsing after a validation error.

`JsonLiteral(String)` is deliberately unvalidated. Supply already escaped ASCII
JSON text; the class adds no quotes and does not repair malformed fragments.
Formatting switches such as `setDecimals` and `noNewlinesInArrays` are global:
configure them before starting concurrent serialization.

Regression and timing fixtures live under `test/json`. They cover literal expected
values and output bytes rather than relying solely on parser/serializer round trips.

# Tadpole config-dispatch fixtures

These are the original inputs for `test.assemble.TadpoleConfigUnitTest`, restored from
Fischl's retained Dori test inputs on November 8, 2026, with her permission. The
copies were compared by freshly computed sha80 values at both locations.

| File | Purpose | sha80 |
|---|---|---|
| k31.config | Select the short-K implementation | 6415a901fb632a9caa84 |
| k63.config | Select the long-K implementation | 431fdfec72e55538ae63 |
| packed.config | Select packed K=64 | 1bd5bd37f59439198fdd |
| forced.config | Force Tadpole2 at K=31 | 7a5560e3608f93d6187d |
| multi.config | Request the K=63,31 phase schedule | 8a8bafb353631bada7ee |

The harness adds its read-input and small-constructor options separately. These
files are literal regression inputs, not a recommended assembly configuration.

Retained origin: `Fischl/BBTools/testdata/tadpole_config/` on Dori. No old source
tree, class files or local-build helper was copied with them.

package prot;

/**
 * Injection seam for every external-tool invocation the Candidate-C build pipeline makes (MAFFT,
 * {@code hmmbuild}, {@code hmmpress}). Injected via CONSTRUCTOR, never a CLI class-name string or
 * shell command string (Elly, 2026-09-01: "do not accept arbitrary class names or shell command
 * strings from CLI") -- the production entry point wires the one explicit
 * {@link RealCommandRunner}; a test constructs a fixture implementation directly and passes it in.
 * <p>
 * {@code argv} is always a literal array passed straight to {@link ProcessBuilder} (never built as
 * or passed through a shell string), so there is no shell-injection surface regardless of what a
 * family/member ID contains -- the exact class of risk this project has already been burned by
 * once (a pipe-bearing rep_id breaking naive shell quoting,
 * {@code project_magqc_identity_grouping_batch_argv_fix}).
 *
 * @author Eru
 */
interface CommandRunner {

	/**
	 * Runs one external command to completion, capturing stdout/stderr to the given files.
	 * <p>
	 * Frozen signature (`CANDIDATE_C_PHASE_AB_GATE_PLAN_v1.md` sec 2, UMP45/Elly co-sealed
	 * 2026-09-01): {@code workingDirectory} is passed to the process launcher explicitly -- never
	 * relied upon implicitly via the launching JVM's own ambient current directory -- so Phase B's
	 * exact {@code hmmpress db.hmmdb} command (no path argument) runs against the correct staging
	 * directory regardless of what directory the calling process happens to be in. A caller never
	 * smuggles a {@code cd} through a shell string; there is no shell string anywhere in this seam.
	 * @param argv The literal command + arguments (never shell-interpreted).
	 * @param workingDirectory The process's working directory -- the caller's own task-private staging
	 *        directory, never a shared or parent path (each builder validates this before calling).
	 * @param stdoutFile Path to write captured stdout.
	 * @param stderrFile Path to write captured stderr.
	 * @return The command's exit code. A launch or I/O failure propagates as a thrown exception --
	 *         it must never be silently mapped to a successful exit code.
	 */
	int run(String[] argv, String workingDirectory, String stdoutFile, String stderrFile) throws Exception;
}

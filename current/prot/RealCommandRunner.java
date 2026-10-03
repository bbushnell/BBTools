package prot;

import java.io.File;

/**
 * Production {@link CommandRunner}: a literal {@code argv} array passed directly to
 * {@link ProcessBuilder} (never a shell string), stdout/stderr redirected to the given files.
 * Every real {@code mafft}/{@code hmmbuild}/{@code hmmpress} invocation in production goes through
 * exactly this class; the only substitution point for tests is the constructor injection at the
 * call site (Elly/UMP45, 2026-09-01), never a runtime class-name switch.
 *
 * @author Eru
 */
final class RealCommandRunner implements CommandRunner {

	@Override
	public int run(final String[] argv, final String workingDirectory, final String stdoutFile,
			final String stderrFile) throws Exception{
		final ProcessBuilder pb=new ProcessBuilder(argv);
		pb.directory(new File(workingDirectory));
		pb.redirectOutput(new File(stdoutFile));
		pb.redirectError(new File(stderrFile));
		final Process p=pb.start();
		return p.waitFor();
	}
}

package prok;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.util.concurrent.TimeUnit;
import fileIO.ByteFile;
import map.ObjectMap;
import parse.LineParser1;
import structures.ByteBuilder;

/** Verifies declared input contents with the canonical sha80 executable.
 * The manifest is an input supplied before the run, never generated from the
 * current inputs and then presented as independent evidence of their identity.
 * Created for the earlier pipeline's declared-input verification. Its old
 * pipeline consumers were retired, but the retained text-I/O owner does not
 * replace manifest-relative path resolution or verifyAllUnder. Preserved as a
 * reusable verifier when a future workflow needs predeclared sha80 manifests.
 * @author Raiden */
public final class NcrnaBoundaryInputs {

	public NcrnaBoundaryInputs(String manifest, String command_) throws IOException{
		this(manifest, command_, new File(manifest).getCanonicalFile().getParentFile());
	}
	/** Build receipts resolve their named source/class entries beneath the installation. */
	public NcrnaBoundaryInputs(String manifest, String command_, File base) throws IOException{
		command=new File(command_).getCanonicalPath();
		require(new File(command).isFile() && new File(command).canExecute(), "Require the canonical native sha80 executable");
		final File parent=base.getCanonicalFile();require(parent.isDirectory(), "Manifest-relative root must exist");
		final ByteFile in=RrnaResourceIO.open(manifest);final LineParser1 p=new LineParser1('\t');
		try{
			for(byte[] row; (row=in.nextLine())!=null;){
				p.set(row);require(p.terms()==2, "Input manifest rows must be sha80<TAB>path");
				final String expected=p.parseString(0), name=p.parseString(1);
				require(expected.matches("[0-9a-f]{20}") && !name.isEmpty(), "Input provenance requires a complete sha80 and path");
				File file=new File(name);if(!file.isAbsolute()){file=new File(parent, name);}
				final String path=file.getCanonicalPath();
				require(file.isFile() && hashes.put(path, expected)==null, "Missing or duplicate canonical input path: "+name);
			}
		}finally{require(!in.close(), "Input manifest read failed");}
		require(hashes.size()>0, "Input provenance manifest must not be empty");
	}

	/** Bind each required original library, graph, corpus and configuration before use. */
	public String verify(String path) throws IOException, InterruptedException{
		final String canonical=new File(path).getCanonicalPath(), expected=hashes.get(canonical);
		require(expected!=null, "Required input is absent from the predeclared sha80 manifest: "+path);
		final String actual=hash(command, canonical);
		require(expected.equals(actual), "Input content differs from its declared provenance: "+path);
		return actual;
	}
	/** Verify only the manifest's explicit entries, all inside the declared tree. */
	public void verifyAllUnder(File root) throws IOException, InterruptedException{
		final java.nio.file.Path parent=root.getCanonicalFile().toPath();
		for(String path:hashes.keys()){if(path!=null){
			require(new File(path).toPath().startsWith(parent), "Build provenance cannot escape its declared installation: "+path);verify(path);
		}}
	}

	/** Control-plane subprocess: the installed native sha80 writes one short line.
	 * No shell, second hash implementation or full digest is involved. */
	public static String hash(String command, String path) throws IOException, InterruptedException{
		require(new File(command).isFile() && new File(path).isFile(), "Hash an existing file with the canonical sha80 executable");
		final Process child=new ProcessBuilder(command, new File(path).getAbsolutePath()).redirectError(ProcessBuilder.Redirect.INHERIT).start();
		try{
			final int code=child.waitFor();if(code!=0){throw new IOException("sha80 exited "+code+" for "+path);}
			final ByteBuilder result=new ByteBuilder(32);
			try(InputStream in=child.getInputStream()){
				for(int value; (value=in.read())>=0;){if(result.length>=22){throw new IOException("sha80 produced more than one short digest line");}result.append((byte)value);}
			}
			final String value=new String(result.toBytes(), StandardCharsets.US_ASCII).trim();
			require(value.matches("[0-9a-f]{20}"), "sha80 did not produce exactly twenty lowercase hexadecimal characters");
			return value;
		}finally{
			if(child.isAlive()){
				child.destroyForcibly();boolean interrupted=Thread.interrupted();
				final long deadline=System.nanoTime()+TimeUnit.SECONDS.toNanos(5);
				try{
					while(child.isAlive() && System.nanoTime()<deadline){
						try{child.waitFor(100, TimeUnit.MILLISECONDS);}catch(InterruptedException e){interrupted=true;}
					}
					if(child.isAlive()){throw new IOException("Native sha80 process did not stop after cancellation");}
				}finally{if(interrupted){Thread.currentThread().interrupt();}}
			}
		}
	}

	private final String command;
	private final ObjectMap<String, String> hashes=new ObjectMap<String, String>(String.class, String.class);
	private static void require(boolean ok, String why){if(!ok){throw new IllegalArgumentException(why);}}
}

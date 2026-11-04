package test.dna;

import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;

import dna.Matrix;

/** Finite matrix-registry controls and failed-initialization retry oracles.
 * Run each case in a fresh JVM whose working directory is the supplied empty work directory.
 * @author Shinobu
 */
public final class MatrixRegistryTest{

	/** Writes only the four named fixture files, refusing pre-existing inputs. */
	public static void main(String[] args) throws IOException{
		if(args.length!=2 || !args[0].startsWith("case=") || !args[1].startsWith("work=")){
			throw new IllegalArgumentException("Require case=controls|retry-malformed work=<current directory>");
		}
		final String mode=args[0].substring(5);
		require(Files.isSameFile(Paths.get("."), Paths.get(args[1].substring(5))),
				"Matrix loads relative filenames; the JVM must run inside the isolated fixture directory");
		if(!mode.equals("controls") && !mode.equals("retry-malformed")){
			throw new IllegalArgumentException(mode);
		}
		for(String file : FILES){require(!Files.exists(Paths.get(file)), "Refusing pre-existing fixture: "+file);}
		write(0, "First", "1", false);
		write(2, "Third", "3", false);
		write(3, "Percentile", "4", false);
		if(mode.equals("controls")){
			write(1, "Second", "2", false);
			controls(1);
		}else{
			write(1, "Second", "not-a-number", false);
			RuntimeException failure=null;
			try{Matrix.get("First");}
			catch(RuntimeException e){failure=e;}
			require(failure!=null, "Initialization must report the incomplete second matrix file");
			require(failure instanceof NumberFormatException,
					"Unexpected failure instead of the deliberately malformed numeric fixture: "+failure);
			System.out.println("Initial failure: "+failure.getClass().getSimpleName());
			write(0, "First", "9", true);
			write(1, "Second", "2", true);
			controls(9);
		}
		System.out.println("PASS "+mode);
	}

	/** Verifies complete loading, aliases, copied subgrids and retry freshness. */
	private static void controls(float first){
		final Matrix second=Matrix.get("Second");
		require(second.grid[0][0]==2 && second.prefix==0, "Second matrix must be loaded after initialization");
		final Matrix initial=Matrix.get("First");
		require(initial.grid[0][0]==first, "A retry must reread the first file rather than reuse a partial cache");
		require(Matrix.get("first")==initial, "Original and lowercase aliases must retain object identity");
		require(Matrix.get("Third").grid[0][0]==3 && Matrix.get("Percentile").grid[0][0]==4,
				"Initialization must load all four input files");
		require(Matrix.keys().size()==8, "Four distinct names must retain both original and lowercase aliases");
		final float[][] copy=initial.subGrid(0, 1);
		copy[0][0]=-1;
		require(initial.grid[0][0]==first, "subGrid must copy rows rather than share mutable cells");
		boolean rejected=false;
		try{Matrix.get("fIrSt");}catch(RuntimeException e){rejected=true;}
		require(rejected, "Existing lookup is exact, not arbitrary-case normalization");
	}

	/** Emits one tiny matrix; replacement is used only for a fixture written by this JVM. */
	private static void write(int index, String name, String value, boolean replace) throws IOException{
		final Path file=Paths.get(FILES[index]);
		final String text="//name: "+name+"\n//size: 1x4\n//prefix: 0\n//count: 1\n{"+
				value+", 2, 3, 4}\n//end\n";
		Files.write(file, text.getBytes(StandardCharsets.US_ASCII), StandardOpenOption.WRITE,
				replace ? StandardOpenOption.TRUNCATE_EXISTING : StandardOpenOption.CREATE_NEW);
	}

	/** Keeps oracle checks active in both assertion modes. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}

	/** The exact default input names used by Matrix.get. */
	private static final String[] FILES={"matrix_build37_N1.txt", "matrix_build37_N2.txt",
			"matrix_build37_N3.txt", "SplicePercentiles_b37_Nov24.txt"};
}

package test.dna;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import dna.ChromosomeArray;
import dna.ChromToFasta;
import fileIO.ReadWrite;
import shared.Shared;

/** Finite FASTA writer checks without modifying reference indexes.
 * The filename-overload producer and verifier run in separate JVM invocations:
 * the producer's non-daemon writer must finish before its JVM exits.
 * @author Sucrose
 */
public final class ChromToFastaContractTest{

	/** Requires case=whole|split|negative-trigger|full|write-overload|verify-overload and work=path. */
	public static void main(String[] args) throws Exception{
		String mode=null;
		Path work=null;
		for(String arg : args){
			if(arg.startsWith("case=")){mode=arg.substring(5);}
			else if(arg.startsWith("work=")){work=Paths.get(arg.substring(5));}
			else{throw new IllegalArgumentException(arg);}
		}
		if(mode==null || work==null){throw new IllegalArgumentException("Require case= and work=");}
		if(mode.equals("verify-overload")){
			requireText(work.resolve("output.fa"), ">7\nNNAC\nNNTG\nNN\n");
		}else{
			if(!(mode.equals("whole") || mode.equals("split") || mode.equals("negative-trigger") ||
					mode.equals("full") || mode.equals("write-overload"))){
				throw new IllegalArgumentException(mode);
			}
			Files.createDirectory(work);
			Shared.setThreads(1);
			ReadWrite.USE_PIGZ=ReadWrite.USE_UNPIGZ=false;
			ChromosomeArray.CHANGE_U_TO_T=ChromosomeArray.CHANGE_DEGENERATE_TO_N=true;
			final ChromosomeArray ca=new ChromosomeArray(7, Shared.PLUS);
			ca.set(0, "NNACNNTGNN");
			final Path output=work.resolve("output.fa");
			if(mode.equals("write-overload")){
				ChromToFasta.writeChrom(ca, output.toString(), 4);
				System.out.println("PRODUCER_RETURNED; verify only after this JVM exits");
				return;
			}
			final Path input=work.resolve("input.chrom");
			ReadWrite.write(ca, input.toString(), false);
			require(Files.isRegularFile(input), "Serialized input was not written");
			if(mode.equals("full")){
				require(Files.exists(Paths.get("/dev/full")), "Linux /dev/full is required");
				RuntimeException failure=null;
				try{ChromToFasta.main(new String[]{input.toString(), "/dev/full", "4"});}
				catch(RuntimeException expected){failure=expected;}
				System.out.println("outputFailure="+(failure==null ? "none" : failure));
				require(failure!=null, "ChromToFasta returned normally after output failure");
			}else{
				//Nonpositive split thresholds promise whole output, including its leading/trailing Ns.
				final String trigger=mode.equals("split") ? "2" : mode.equals("negative-trigger") ? "-1" : "0";
				ChromToFasta.main(new String[]{input.toString(), output.toString(), "4",
						trigger});
				requireText(output, mode.equals("split") ? ">1\nAC\n>2\nTG\n" : ">7\nNNAC\nNNTG\nNN\n");
			}
		}
		System.out.println("PASS "+mode);
	}

	/** Prints literal evidence before requiring exact FASTA bytes. */
	private static void requireText(Path path, String expected) throws Exception{
		final String actual=new String(Files.readAllBytes(path), StandardCharsets.US_ASCII);
		System.out.println("ACTUAL_BEGIN\n"+actual+"ACTUAL_END");
		require(actual.equals(expected), "FASTA differs: expected "+expected.replace('\n', '|')+
				" observed "+actual.replace('\n', '|'));
	}

	private static void require(boolean ok, String message){
		if(!ok){throw new AssertionError(message);}
	}
}

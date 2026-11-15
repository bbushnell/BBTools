package test.pacbio;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.Locale;

import dna.ChromosomeArray;
import dna.Data;
import fileIO.ReadWrite;
import fileIO.TextStreamWriter;
import pacbio.CalcCoverageFromSites;
import structures.CoverageArray;
import structures.CoverageArray3;

/**
 * Finite coverage-mode regression checks using two literal in-memory chromosomes.
 * Reflection installs an already-loaded test genome; genome-file loading is not
 * under test. Both public processing methods still parse actual SiteScoreR files,
 * aggregate and print their own summaries; serialized coverage is read back.
 * @author Brian Bushnell, Shinobu
 */
public class CalcCoverageFromSitesTest{

	/** Runs historical observations or the approved identical-mode contract. */
	public static void main(String[] args) throws Exception{
		File out=null;
		boolean fixed=false;
		for(String arg : args){
			if(arg.startsWith("outdir=")){out=new File(arg.substring(7));}
			else if(arg.equals("phase=fixed")){fixed=true;}
			else if(arg.equals("phase=baseline")){fixed=false;}
			else{throw new IllegalArgumentException(arg);}
		}
		require(out!=null && out.mkdir(), "Use a fresh output directory");
		Locale.setDefault(Locale.ROOT);
		CalcCoverageFromSites.MIN_END_DIST=4;
		System.out.println("case\tmode\ttotal_depth\tcorrect_depth\tresult");
		check(out, "ends", 0, 11, 1, 1, 1, new int[] {0,0,0,0,1,1,1,1,0,0,0,0}, 4, fixed);
		check(out, "deep", 0, 11, 70000, 66000, 1,
			new int[] {0,0,0,0,70000,70000,70000,70000,0,0,0,0}, 264000, fixed);
		check(out, "clipped", -4, 15, 1, 1, 1, new int[] {1,1,1,1,1,1,1,1,1,1,1,1}, 12, fixed);
		check(out, "short", 3, 7, 1, 1, 1, new int[12], 0, fixed);
		check(out, "threshold", 0, 11, 1, 1, 2, new int[] {0,0,0,0,1,1,1,1,0,0,0,0}, 4, fixed);
		check(out, "empty", 0, 11, 0, 0, 1, new int[12], 0, fixed);
		if(fixed){
			check(out, "outside-high", Integer.MAX_VALUE-1, Integer.MAX_VALUE, 1, 1, 1, new int[12], 0, true);
			check(out, "outside-low", Integer.MIN_VALUE, Integer.MIN_VALUE+1, 1, 1, 1, new int[12], 0, true);
			checkOverflow();
		}
	}

	/** Exercises the integer boundary without constructing billions of input sites. */
	private static void checkOverflow() throws Exception{
		final Method increment=CalcCoverageFromSites.class.getDeclaredMethod("incrementCoverage",
			CoverageArray3.class, CoverageArray3.class, boolean.class, int.class, int.class, int.class);
		increment.setAccessible(true);
		final CoverageArray3 total=new CoverageArray3(1, 1), correct=new CoverageArray3(1, 1);
		total.set(0, Integer.MAX_VALUE-1);
		correct.set(0, Integer.MAX_VALUE-2);
		increment.invoke(null, total, correct, true, 0, 0, 1);
		require(total.get(0)==Integer.MAX_VALUE && correct.get(0)==Integer.MAX_VALUE-1, "Last representable depth was not retained");
		boolean threw=false;
		try{increment.invoke(null, total, correct, true, 0, 0, 1);}
		catch(InvocationTargetException e){
			if(!(e.getCause() instanceof ArithmeticException)){throw e;}
			threw=true;
		}
		require(threw && total.get(0)==Integer.MAX_VALUE && correct.get(0)==Integer.MAX_VALUE-1,
			"Depth overflow must throw before wrapping or saturating a new observation");
		System.out.println("overflow\thelper\t2147483647\t2147483646\tPASS");
	}

	/** Writes a literal repeated site, then compares independently specified counts. */
	private static void check(File root, String name, int start, int stop, int count, int goodCount,
		int mincoverage, int[] expectedDepth, long expectedGood, boolean fixed) throws Exception{
		final File input=new File(root, name+".sites");
		final TextStreamWriter writer=new TextStreamWriter(input.getPath(), false, false, false);
		writer.start();
		for(int i=0; i<count; i++){
			writer.println((i<goodCount ? "*" : "")+"1,0,"+start+","+stop+","+(stop-start+1)+","+i+",0,00,0,0");
		}
		require(!writer.poisonAndWait(), "Could not write "+input);
		String previous=null;
		long expectedTotal=0;
		for(int depth : expectedDepth){expectedTotal+=depth;}
		for(boolean output : new boolean[] {false, true}){
			installGenome();
			final String mode=output ? "write" : "summary";
			final String pattern=new File(root, name+"-chr#.ca").getPath();
			final ByteArrayOutputStream bytes=new ByteArrayOutputStream();
			final PrintStream original=System.out;
			try(PrintStream capture=new PrintStream(bytes, true, "UTF-8")){
				System.setOut(capture);
				if(output){CalcCoverageFromSites.processAndWrite(input.getPath(), 1, mincoverage, pattern);}
				else{CalcCoverageFromSites.process(input.getPath(), 1, mincoverage);}
			}finally{System.setOut(original);}
			ReadWrite.waitForWritingToFinish();
			final String report=new String(bytes.toByteArray(), StandardCharsets.UTF_8);
			Files.write(new File(root, name+"-"+mode+".txt").toPath(), bytes.toByteArray());
			final long total=total(report, "Avg Coverage:");
			final long good=total(report, "Correct Coverage:");
			if(fixed){
				require(!report.contains("NaN") && !report.contains("Infinity"), name+": nonfinite report");
				require(total==expectedTotal && good==expectedGood, name+": depth totals "+total+","+good);
				require(value(report, "Avg Coverage:").equals(String.format(Locale.ROOT, "%.3f", expectedTotal/20d)), name+": denominator is not 20 bases");
				if(previous!=null){require(previous.equals(report), name+": modes disagree");}
				if(name.equals("deep")){
					require(value(report, "mostlyCorrect:").equals("33.333%"), "Deep mostly-correct bases misclassified");
					require(value(report, "onlyCorrect:").equals("0.000%"), "Deep mixture called only correct");
				}
				if(name.equals("threshold")){
					require(total(report, "Avg Coverage Base:")==0 && total(report, "Avg Coverage N:")==0,
						"Subthreshold positions entered covered-only totals");
				}
			}else if(name.equals("ends")){
				require(total==(output ? 3 : 12) && good==total, "Historical endpoint baseline changed");
			}else if(name.equals("deep")){
				require(total==(output ? 196605 : 1524) && good==(output ? 381 : 1524), "Historical saturation baseline changed");
			}
			if(output){
				for(int chrom=1; chrom<=2; chrom++){
					final CoverageArray coverage=ReadWrite.read(CoverageArray.class, pattern.replace("#", ""+chrom), false);
					final int length=chrom==1 ? 12 : 8;
					for(int i=0; i<length; i++){
						final int expected=chrom==2 ? 0 : expectedDepth[i];
						require(coverage.get(i)==(fixed ? expected : Math.min(expected, 65535)), name+": serialized depth at "+chrom+":"+i);
					}
				}
			}
			previous=report;
			System.out.println(name+"\t"+mode+"\t"+total+"\t"+good+"\t"+(fixed ? "PASS" : "OBSERVED"));
		}
	}

	/** Installs exactly 20 positions: nine defined bases and eleven Ns. */
	private static void installGenome() throws Exception{
		final Field loaded=Data.class.getDeclaredField("genome_set_to");
		loaded.setAccessible(true);
		loaded.setInt(null, 1);
		Data.GENOME_BUILD=1;
		Data.numChroms=2;
		Data.chromLengths=new int[] {0,12,8};
		Data.chromosomePlusMatrix=new ChromosomeArray[3];
		final String[] bases={"", "NAAAANAAAAAN", "NNNNNNNN"};
		for(int chrom=1; chrom<=2; chrom++){
			final ChromosomeArray array=new ChromosomeArray(chrom, (byte)0, 0, bases[chrom].length()-1);
			Arrays.fill(array.array, (byte)'N');
			final byte[] literal=bases[chrom].getBytes(StandardCharsets.US_ASCII);
			System.arraycopy(literal, 0, array.array, 0, literal.length);
			Data.chromosomePlusMatrix[chrom]=array;
		}
		assert(Data.chromLengths[1]+Data.chromLengths[2]==20) : "Fixture oracle uses a 20-position denominator";
	}

	/** Gets the first numeric report field for an exact label. */
	private static String value(String report, String label){return fields(report, label)[0];}

	/** Gets the integer count after the formatted average. */
	private static long total(String report, String label){return Long.parseLong(fields(report, label)[1]);}

	/** Splits only the selected human-readable report line. */
	private static String[] fields(String report, String label){
		for(String line : report.split("\\R")){
			if(line.startsWith(label)){return line.substring(label.length()).trim().split("\\s+");}
		}
		throw new AssertionError("Missing report label "+label);
	}

	/** Throws independently of assertion-enable flags. */
	private static void require(boolean condition, String message){
		if(!condition){throw new AssertionError(message);}
	}
}

package test.prok;

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import fileIO.ByteFile;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import prok.NcrnaMemberLengthBounds;
import prok.RrnaResourceIO;

/** Export extrema rather than model length/median; keep singleton models and
 * reject a changed final library instead of attaching bounds to row numbers.
 * @author Brian Bushnell, Ganyu
 */
public final class NcrnaMemberLengthBoundsTest {
	public static void main(String[] args)throws Exception{
		check(args.length==1, "Supply a fresh fixture directory");final Path root=Paths.get(args[0]);Files.createDirectory(root);
		write(root.resolve("ref.fa"), ">cobalamin_RF00174_original n=3\nACGTAC\n>cobalamin_RF00174_singleton n=1\nACGT\n");
		final String original="cobalamin_RF00174_original\tRF00174\t6\t3\t11\t11\t12\t12\t19\t1\tPURE\t3\t1\t0\t0\t0\t0\t0\t0\n";
		final String singleton="cobalamin_RF00174_singleton\tRF00174\t4\t1\t4\t4\t4\t4\t4\t1\tPURE\t1\t1\t0\t0\t0\t0\t0\t0\n";
		write(root.resolve("members.tsv"), NcrnaMemberLengthBounds.HEADER+'\n'+singleton+original);
		run(root, "valid.tsv");final ByteFile in=RrnaResourceIO.open(root.resolve("valid.tsv").toString());
		try{check(new String(in.nextLine()).equals("model\tmodel_length\tmin_length\tmax_length"), "Resource uses the existing length_bounds schema");
			check(new String(in.nextLine()).equals("cobalamin_RF00174_original n=3\t6\t8\t23"), "Bounds use existing outward floor(0.8*11)=8 and ceil(1.2*19)=23, not median12 or consensus length6");
			check(new String(in.nextLine()).equals("cobalamin_RF00174_singleton n=1\t4\t3\t5") && in.nextLine()==null, "Added singleton remains covered and output follows complete caller identity");
		}finally{check(!in.close(), "Fixture output read failed");}
		final prok.NcrnaFamily loaded=NcrnaMemberIdentityCutoffsTest.loadProfile(root.resolve("ref.fa"), root.resolve("valid.tsv"));
		check(loaded.modelThresholds.minimum(0)==8 && loaded.modelThresholds.maximum(1)==5, "Actual caller name loader must bind exported length limits with metadata intact");
		write(root.resolve("members.tsv"), NcrnaMemberLengthBounds.HEADER+'\n'+original);fails(root, "missing.tsv");
		write(root.resolve("members.tsv"), NcrnaMemberLengthBounds.HEADER+'\n'+original+original);fails(root, "duplicate.tsv");
		write(root.resolve("members.tsv"), NcrnaMemberLengthBounds.HEADER+'\n'+original+singleton.replace("\t4\t1\t", "\t5\t1\t"));fails(root, "changed.tsv");
		System.out.println("NCRNA_MEMBER_LENGTH_BOUNDS_TEST_PASS");
	}
	private static void run(Path root, String output){NcrnaMemberLengthBounds.main(new String[]{"in="+root.resolve("members.tsv"), "ref="+root.resolve("ref.fa"), "out="+root.resolve(output), "family=RF00174", "expected=2"});}
	private static void fails(Path root, String output){boolean failed=false;try{run(root, output);}catch(IllegalArgumentException expected){failed=true;}check(failed && !Files.exists(root.resolve(output)), "Invalid model memberships must fail before emitting a usable bounds resource");}
	/** Rewrites are deliberate fixture mutations inside the newly created test directory. */
	private static void write(Path file, String value){
		final ByteStreamWriter out=new ByteStreamWriter(FileFormat.testOutput(file.toString(), FileFormat.TEXT, null, true, true, false, false));
		out.start();try{out.print(value);}finally{RrnaResourceIO.finish(out);}
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}

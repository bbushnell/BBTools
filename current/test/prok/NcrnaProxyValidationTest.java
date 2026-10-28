package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;

import fileIO.ByteFile;
import parse.LineParser1;
import prok.NcrnaProxyValidation;

/** Verifies arm conservation, no-hit retention and exact source-metadata joins.
 * Uses the already-executed synthetic sampler fixture, never a real sample.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaProxyValidationTest {
	public static void main(String[] args)throws Exception{
		if(args.length!=2){throw new IllegalArgumentException("Supply fresh fixture directory and synthetic sample ledger");}
		final Path root=Paths.get(args[0]), sample=Paths.get(args[1]);Files.createDirectory(root);
		final ArrayList<String> rows=new ArrayList<String>();
		final java.util.List<String> selected=Files.readAllLines(sample, StandardCharsets.US_ASCII);
		check(selected.size()==9, "Expected the sampler fixture's2+6 selected loci");
		for(int i=1;i<selected.size();i++){
			final String[] fields=selected.get(i).split("\t", -1);check(fields.length==29, "Synthetic selection schema changed");
			rows.add(String.join("\t", Arrays.copyOf(fields, 21))+"\tNO_REPORTED_HIT\t0\t0\tNA\tNA\tNA\tNA\tNA\tNA\t0\tNA\tNA\tNONE\n");
		}
		final Path audit=root.resolve("audit.tsv");write(audit, rows);run(sample, audit, root.resolve("good"));
		final ByteFile in=ByteFile.makeByteFile(root.resolve("good.joined.tsv.gz").toString(), false);final LineParser1 p=new LineParser1('\t');int np=0, nc=0;
		try{
			check(in.nextLine()!=null, "Joined evidence header missing");
			for(byte[] row;(row=in.nextLine())!=null;){
				p.set(row);check(p.terms()==42 && p.termEquals("NO_REPORTED_HIT", 21), "Every no-hit source must retain its evidence and selection fields");
				if(p.termEquals("PROBABILITY", 38)){np++;}else{check(p.termEquals("CHALLENGE", 38), "Unknown joined arm");nc++;}
			}
		}finally{check(!in.close(), "Joined evidence read failed");}
		check(np==2 && nc==6, "Probability and enriched challenge denominators must remain separate");
		for(String name:new String[]{"good.summary.tsv", "good.phyla.tsv"}){
			final ByteFile table=ByteFile.makeByteFile(root.resolve(name).toString(), false);int total=0;
			try{table.nextLine();for(byte[] row;(row=table.nextLine())!=null;){
				p.set(row);final int offset=name.contains("phyla")?2:0;
				final int n=p.parseInt(offset+5), groups=p.parseInt(offset+6);check(groups>0 && groups<=n, "Source-group counts cannot exceed contributing loci");total+=n;
			}}finally{check(!table.close(), "Summary read failed");}
			check(total==8, "Each summary view must conserve all selected loci");
		}
		final String original=rows.get(0);final String[] altered=original.split("\t", -1);altered[1]="wrong_source_genome";
		rows.set(0, String.join("\t", altered));write(audit, rows);
		boolean failed=false;try{run(sample, audit, root.resolve("wrong_source"));}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "A matching locus ID cannot conceal different original source metadata");
		rows.set(0, original);rows.remove(rows.size()-1);write(audit, rows);
		failed=false;try{run(sample, audit, root.resolve("missing"));}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "Missing no-hit rows must fail rather than improve apparent support");
		System.out.println("NCRNA_PROXY_VALIDATION_TEST_PASS");
	}
	private static void run(Path sample, Path audit, Path prefix){
		NcrnaProxyValidation.main(new String[]{"sample="+sample, "audit="+audit, "outprefix="+prefix, "expectedprobability=2", "expectedchallenge=6"});
	}
	private static void write(Path path, ArrayList<String> rows)throws Exception{
		final StringBuilder b=new StringBuilder(NcrnaProxyValidation.AUDIT_HEADER).append('\n');for(String row:rows){b.append(row);}
		Files.write(path, b.toString().getBytes(StandardCharsets.US_ASCII));
	}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}

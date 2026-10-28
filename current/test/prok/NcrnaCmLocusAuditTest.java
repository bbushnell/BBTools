package test.prok;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;

import fileIO.ByteFile;
import parse.LineParser1;
import prok.NcrnaCmLocusAudit;

/** Attribution fixtures distinguish a neighboring hit from support of its target.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaCmLocusAuditTest {

	public static void main(String[] args)throws Exception{
		if(args.length!=1){throw new IllegalArgumentException("Supply an unused fixture directory");}
		final Path root=java.nio.file.Paths.get(args[0]);Files.createDirectory(root);
		final String[] ids={"none", "neighbor", "opposite", "partial", "truncated", "full", "multiple", "low", "touch", "long"};
		final StringBuilder rows=new StringBuilder(NcrnaCmLocusAudit.HEADER).append('\n');
		for(String id:ids){
			final int length=id.equals("long")?477:190;
			rows.append(id).append("\tgenome\tGCF_1.1\t1\tbacteria\tp\tcontig\t-\t1001-").append(1000+length).append('\t').append(length)
				.append("\t1\tfeature\tfalse\tfalse\tANNOTATION_FULL_COPY\t2000\tfalse\t0\t200\t200\t1\n");
		}
		write(root.resolve("loci.tsv"), rows.toString());
		final String hits=hit("neighbor", 1, 189, 1, 189, "+", "no", "50")
			+hit("opposite", 1, 189, 389, 201, "-", "no", "50")
			+hit("partial", 2, 189, 201, 389, "+", "no", "50")
			+hit("truncated", 1, 189, 201, 389, "+", "5'", "50")
			+hit("full", 1, 189, 201, 389, "+", "no", "31.5")
			+hit("multiple", 1, 189, 201, 389, "+", "no", "50")
			+hit("multiple", 1, 189, 210, 398, "+", "no", "100")
			+hit("low", 1, 189, 201, 389, "+", "no", "31.4")
			+hit("touch", 1, 189, 390, 578, "+", "no", "80")
			+hit("long", 1, 189, 201, 406, "+", "no", "55.8")
			+"# Program: cmsearch\n# Version: 1.1.5 (Sep 2023)\n# [ok]\n";
		write(root.resolve("hits.tbl"), hits);
		run(root, "audit", 10);
		final String[] states={"NO_REPORTED_HIT", "NO_OVERLAPPING_SENSE_HIT", "NO_OVERLAPPING_SENSE_HIT",
			"CM_SUPPORTED_INCOMPLETE_OR_UNRESOLVED", "CM_SUPPORTED_INCOMPLETE_OR_UNRESOLVED", "CM_FULL_SPAN_SUPPORTED",
			"MULTIPLE_OVERLAPPING_HITS", "BELOW_GA_OR_NOT_INCLUDED", "PARTIAL_ANNOTATION_OVERLAP_UNRESOLVED", "CM_FULL_SPAN_SUPPORTED"};
		final ByteFile in=ByteFile.makeByteFile(root.resolve("audit.loci.tsv.gz").toString(), false);
		final LineParser1 p=new LineParser1('\t');int count=0;
		try{
			check(in.nextLine()!=null, "Output header missing");
			for(byte[] line;(line=in.nextLine())!=null;){
				p.set(line);check(count<ids.length && p.termEquals(ids[count], 0), "Every original locus must retain order and identity");
				check(p.termEquals(states[count], 21), "Wrong evidence class for "+ids[count]);
				if(count==5){check(p.parseInt(27)==389 && p.parseInt(30)==189 && p.parseInt(31)==0 && p.parseInt(32)==1,
					"Inclusive sense-oriented coordinates and excess lengths must not be flipped by original minus strand");}
				if(count==8){check(p.parseInt(30)==1 && p.termEquals("PARTIAL_OVERLAP", 33),
					"A neighboring full model touching one annotated base must remain unresolved");}
				count++;
			}
		}finally{check(!in.close(), "Fixture output read failed");}
		check(count==10, "No-hit and unsupported annotations must not disappear");
		final ByteFile profile=ByteFile.makeByteFile(root.resolve("audit.lengths.tsv").toString(), false);int profiles=0;
		try{
			check(profile.nextLine()!=null, "Length profile header missing");
			for(byte[] line;(line=profile.nextLine())!=null;){
				p.set(line);check(p.terms()==16, "Length profile preserves its16-field schema");
				if(p.termEquals("SOURCE_FULL_CM_FULL", 0) && p.termEquals("HIT_WITHIN_ANNOTATION", 2)){
					check(p.parseInt(4)==2 && p.termEquals("long", 13) && p.parseInt(15)==271, "Reference extrema retain the long annotation's identity and271-base excess");
					check(p.parseInt(9)==(p.termEquals("SOURCE", 3)?477:206), "Source length and CM-hit length must remain separate metrics");profiles++;
				}
			}
		}finally{check(!profile.close(), "Length profile read failed");}
		check(profiles==2, "Both source and CM-hit reference distributions are required");
		boolean failed=false;
		try{run(root, "wrong_denominator", 11);}catch(IllegalArgumentException expected){failed=true;}
		check(failed, "A changed declared input population must fail loudly");
		System.out.println("NCRNA_CM_LOCUS_AUDIT_TEST_PASS");
	}

	private static void run(final Path root, final String name, final int count){
		NcrnaCmLocusAudit.main(new String[]{"loci="+root.resolve("loci.tsv"), "hits="+root.resolve("hits.tbl"),
			"outprefix="+root.resolve(name), "expected="+count, "modelid=RF00174", "family=cobalamin", "clen=189", "gabits=31.5"});
	}
	private static String hit(String id, int from, int to, int start, int stop, String strand, String trunc, String score){
		return id+" - Cobalamin RF00174 cm "+from+' '+to+' '+start+' '+stop+' '+strand+' '+trunc+" 1 0.5 0.0 "+score+" 1e-12 ! fixture\n";
	}
	private static void write(Path path, String text)throws Exception{Files.write(path, text.getBytes(StandardCharsets.US_ASCII));}
	private static void check(boolean ok, String why){if(!ok){throw new AssertionError(why);}}
}

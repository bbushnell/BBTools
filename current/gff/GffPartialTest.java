package gff;

import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashSet;
import java.util.Arrays;

import fileIO.ByteFile;
import prok.AnalyzeGenes;
import prok.GeneModel;
import shared.Shared;

/** Partial-attribute parsing, real RNA-training input and CutGff output verification.
 * Does not claim the later caller/writer or in-memory Orf migration is implemented.
 * @author Keqing
 */
public final class GffPartialTest {

	public static void main(String[] args) throws Exception{
		Shared.setThreads(1);
		ByteFile.FORCE_MODE_BF1=true;
		GffLine.parseAttributes=true;
		if(args.length==2 && args[0].equals("prepare")){
			testReader();
			prepare(Paths.get(args[1]));
			testTraining(Paths.get(args[1]));
			System.out.println("PASS GffPartialTest reader and RNA training input");
		}else if(args.length==3 && args[0].equals("verify")){
			verify(Paths.get(args[1]), args[2]);
		}else{throw new IllegalArgumentException("Usage: GffPartialTest prepare DIR | verify FASTA filtered|all|banned|required");}
	}

	private static void testReader(){
		for(String attr : new String[]{".", "ID=x", "partial=00", "partial=false", "Note=partial=true",
			"notpartial=true", "partiality=true", "Note=partial%3Dtrue", "partial=00;Note=partial=true"}){
			check(!row("test", '+', attr).partial(), "Complete/unrelated attribute read as partial: "+attr);
		}
		for(String flag : new String[]{"true", "10", "01", "11"}){
			for(char strand : new char[]{'+', '-'}){
				check(row("test", strand, "ID=x;partial="+flag+";product=5S ribosomal RNA").partial(), "Lost partial flag "+flag);
			}
		}
		final GffLine missing=row("test", '+', "ID=x");
		missing.attributes=null;
		check(!missing.partial(), "Absent attributes must not imply partial");
		for(String attr : new String[]{"partial=", "partial=1", "partial=02", "partial=100", "partial=trueExtra",
			"partial=TRUE", "partial=10,01", "partial=00;partial=10", "partial=true;partial=true"}){
			boolean rejected=false;
			try{row("test", '+', attr).partial();}
			catch(IllegalArgumentException e){
				check(e.getMessage().contains("partial"), "Failure must identify invalid partiality: "+e);
				rejected=true;
			}
			check(rejected, "Malformed partiality must not enter training silently: "+attr);
		}
	}

	private static GffLine row(final String seq, final char strand, final String attributes){
		return new GffLine((seq+"\ttest\trRNA\t51\t170\t.\t"+strand+"\t.\t"+attributes).getBytes(StandardCharsets.US_ASCII));
	}

	private static void prepare(final Path dir) throws Exception{
		Files.createDirectories(dir);
		final StringBuilder fasta=new StringBuilder(), gff=new StringBuilder("##gff-version 3\n");
		final String[] flags={"", ";partial=00", ";partial=false", ";partial=10", ";partial=01", ";partial=11", ";partial=true", ";Note=partial=true"};
		for(int i=0; i<flags.length; i++){
			for(char strand : new char[]{'+', '-'}){
				final String id="f"+i+(strand=='+' ? "p" : "m");
				fasta.append('>').append(id).append('\n');
				for(int j=0; j<100; j++){fasta.append("ACGT");}
				fasta.append('\n');
				final String marker=(i==0 ? ";keepout=yes" : i==1 ? ";wanted=yes" : "");
				gff.append(row(id, strand, "ID="+id+";product=5S ribosomal RNA"+flags[i]+marker).toString()).append('\n');
			}
		}
		Files.write(dir.resolve("input.fna"), fasta.toString().getBytes(StandardCharsets.US_ASCII));
		Files.write(dir.resolve("input.gff"), gff.toString().getBytes(StandardCharsets.US_ASCII));
	}

	/** The real GeneModel.processRNA path consumes literal GFF. Caller output remains a later gate. */
	private static void testTraining(final Path dir) throws Exception{
		final Field align=AnalyzeGenes.class.getDeclaredField("alignRibo");
		align.setAccessible(true);
		final boolean old=align.getBoolean(null);
		final byte[] before=Files.readAllBytes(dir.resolve("input.gff"));
		try{
			// This fixture tests annotation eligibility, not homology of the invented RNA sequence.
			align.setBoolean(null, false);
			final GeneModel model=new GeneModel(true);
			check(!model.process(dir.resolve("input.fna").toString(), dir.resolve("input.gff").toString()), "RNA training reported an input error");
			final Object stats=GeneModel.class.getField("stats5S").get(model);
			final Field count=stats.getClass().getDeclaredField("lengthCount");
			count.setAccessible(true);
			check(count.getLong(stats)==8, "RNA trainer must retain 8 complete rows and exclude the 8 partial rows; got "+count.getLong(stats));
			check(Arrays.equals(before, Files.readAllBytes(dir.resolve("input.gff"))), "Training exclusion must not remove partial calls from the source GFF");
		}finally{align.setBoolean(null, old);}
	}

	/** Assert identities, not merely record counts; legacy ST's known extraction-length bug is separate. */
	private static void verify(final Path fasta, final String mode) throws Exception{
		check(mode.equals("filtered") || mode.equals("all") || mode.equals("banned") || mode.equals("required"), "Unknown verification mode "+mode);
		final HashSet<String> expected=new HashSet<String>(), observed=new HashSet<String>();
		for(int i=0; i<8; i++){
			boolean retain=(mode.equals("all") || i==0 || i==1 || i==2 || i==7);
			if(mode.equals("banned") && i==0){retain=false;}
			if(mode.equals("required")){retain=(i==1);}
			if(retain){expected.add("f"+i+"p"); expected.add("f"+i+"m");}
		}
		final ByteFile bf=ByteFile.makeByteFile1(fasta.toString(), false);
		try{
			for(byte[] line=bf.nextLine(); line!=null; line=bf.nextLine()){
				if(line.length==0 || line[0]!='>'){continue;}
				final String header=new String(line, StandardCharsets.US_ASCII);
				final int start=header.indexOf("ID=");
				check(start>=0, "CutGff header lost fixture identity: "+header);
				final int end=header.indexOf(';', start);
				check(end>start, "Malformed fixture identity: "+header);
				final String id=header.substring(start+3, end);
				check(observed.add(id), "Duplicate extracted feature: "+id);
			}
		}finally{check(!bf.close(), "Error closing extraction output "+fasta);}
		check(observed.equals(expected), "Wrong "+mode+" extraction: observed="+observed+", expected="+expected);
		System.out.println("PASS CutGff partial filtering "+fasta.getFileName()+" mode="+mode+" records="+observed.size());
	}

	private static void check(final boolean condition, final String reason){
		if(!condition){throw new AssertionError(reason);}
	}
}

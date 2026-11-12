package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.SummarizeNcrnaSweep;
import shared.Shared;
import structures.ByteBuilder;

/** Literal metadata oracles for multicopy coverage and overlapping support.
 * FASTA and GFF paths deliberately do not exist. No caller or grader is run.
 * @author Brian Bushnell, Raiden
 */
public final class SummarizeNcrnaSweepTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, SummarizeNcrnaSweepTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null,
				"Supply out=fresh_fixture_directory with no input or overwrite");
			new SummarizeNcrnaSweepTest(new File(parser.out1).getAbsoluteFile()).run();
			pp.outstream.println("NCRNA_SWEEP_COVERAGE_TEST_PASS");
		}finally{
			Shared.closeStream(pp.outstream);
		}
	}

	private SummarizeNcrnaSweepTest(final File root_){
		root=root_;
		check(root.mkdir(), "Fixture directory must be fresh under an existing parent");
	}

	/*--------------------------------------------------------------*/
	/*----------------           Fixtures           ----------------*/
	/*--------------------------------------------------------------*/

	private void run(){
		final ByteBuilder reference=new ByteBuilder("input_id\tfasta\tgff\n");
		for(final String id:IDS){reference.append(id).tab().append(fasta(id)).tab().append(path("missing"+id+".gff")).nl();}
		write("reference.tsv", utf8(reference));
		write("members.tsv", "family\tmember\treference_model\nfam\tmajor\tRF00001\nfam\tminor%α\tRF00002\n");
		final ByteBuilder changes=new ByteBuilder("input_id\tpolicy\tchange\treference_id\tseqid\tstart\tstop\tstrand\n");
		final ByteBuilder loci=new ByteBuilder("input_id\tparent_id\tfamily\tseqid\tstart\tstop\tstrand\tstatus\tscored_members\tmarginal_members\n");
		locus(changes, loci, "1", "p%α", "LOST", "major,minor%α", ".");
		locus(changes, loci, "1", "p2", "LOST", "major", "minor%α");
		locus(changes, loci, "2", "p%α", "GAINED", "minor%α", ".");
		locus(changes, loci, "2", "p2", "GAINED", "major", ".");
		locus(changes, loci, "3", "before", "LOST", "major", ".");
		locus(changes, loci, "3", "after", "GAINED", "major", ".");
		locus(changes, loci, IDS[3], "only", "LOST", "minor%α", ".");
		final String changeText=utf8(changes), locusText=utf8(loci);
		final String counts="input_id\tpolicy\tTP\tlost_baseline_TP\tgained_TP\n"
			+"1\tbaseline\t2\t0\t0\n1\tchanged\t0\t2\t0\n1\tunchanged\t2\t0\t0\n"
			+"2\tbaseline\t0\t0\t0\n2\tchanged\t2\t0\t2\n2\tunchanged\t0\t0\t0\n"
			+"3\tbaseline\t2\t0\t0\n3\tchanged\t2\t1\t1\n3\tunchanged\t2\t0\t0\n"
			+IDS[3]+"\tbaseline\t1\t0\t0\n"+IDS[3]+"\tchanged\t0\t1\t0\n"+IDS[3]+"\tunchanged\t1\t0\t0\n";
		final String summary="# native metadata fixture\npolicy\tTP\tlost_baseline_TP\tgained_TP\n"
			+"baseline\t5\t0\t0\nchanged\t4\t4\t3\nunchanged\t5\t0\t0\n";
		final String prefix=invoke(counts, changeText, locusText, summary, null);
		final String report=read(prefix+".summary.tsv");
		check(report.endsWith("changed\t4\t5\t4\t4\t3\t4\t1\t1\t3\t2\t2\t1\t1\t1\n"
			+"unchanged\t4\t5\t5\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"),
			"Literal totals must separate seven incidences from four inputs and two lost-coverage inputs");
		final String inputs=read(prefix+".inputs.tsv");
		check(inputs.contains("changed\t1\t2\t0\t2\t0\ttrue\tfalse\t"+fasta("1")+"\n"), "Two losses erase one input's coverage");
		check(inputs.contains("changed\t2\t0\t2\t0\t2\tfalse\ttrue\t"+fasta("2")+"\n"), "Two gains recover one input's coverage");
		check(inputs.contains("changed\t3\t2\t2\t1\t1\tfalse\tfalse\t"+fasta("3")+"\n"), "Equal TP can still contain paired changes");
		final String members=read(prefix+".members.tsv");
		check(members.contains("changed\tLOST\tmajor\tRF00001\tfalse\t3\t2\t0\t0\tNONADDITIVE_REFERENCE_SUPPORT\n"), "Major lost incidences overlap minority support");
		check(members.contains("changed\tLOST\tminor%25α\tRF00002\ttrue\t2\t2\t1\t1\tNONADDITIVE_REFERENCE_SUPPORT\n"), "Marginal evidence must not inflate scored minority support");
		final String parents=read(prefix+".loci.tsv");
		check(parents.contains("\tp%25α\tfam\tseq%25α\t10\t30\t+\t"), "Percent and UTF-8 identities must survive the join exactly once");
		check(parents.contains("FIRST_DETECTION_RECOVERED") && parents.contains("OTHER_TARGET_DETECTIONS_REMAIN")
			&& parents.contains("LAST_DETECTED_COPY_LOST") && parents.contains("ADDITIONAL_DETECTION_IN_COVERED_INPUT"), "All four coverage classes must be exercised");
		check(!new File(prefix+".strata.tsv").exists(), "Default reporting must create no optional category output");
		final ByteBuilder taxonomy=new ByteBuilder("input_id\tfasta\tcategory\n");
		final String[] categories={"CORPUS_LOWER_QUARTILE", "OTHER_CLASSIFIED", "UNKNOWN", "UNCLASSIFIED"};
		for(int i=0; i<IDS.length; i++){taxonomy.append(IDS[i]).tab().append(fasta(IDS[i])).tab().append(categories[i]).nl();}
		taxonomy.append("off_panel\tmissing.fa\tUNKNOWN\n");
		final String taxonomyText=utf8(taxonomy);
		final String categorized=invoke(counts, changeText, locusText, summary, null, taxonomyText);
		for(final String suffix:new String[]{".loci.tsv", ".inputs.tsv", ".summary.tsv", ".members.tsv"}){
			check(read(prefix+suffix).equals(read(categorized+suffix)), "Taxonomy observation must preserve every default report: "+suffix);
		}
		final String strata=read(categorized+".strata.tsv");
		check(strata.contains("changed\tCORPUS_LOWER_QUARTILE\t0\t1\n")
			&& strata.contains("changed\tOTHER_CLASSIFIED\t1\t0\n")
			&& strata.contains("changed\tUNKNOWN\t0\t0\n")
			&& strata.contains("changed\tUNCLASSIFIED\t0\t1\n"), "Literal category counts must count inputs, not lost/gained copies");
		for(final String category:categories){
			check(strata.contains("baseline\t"+category+"\t0\t0\n")
				&& strata.contains("unchanged\t"+category+"\t0\t0\n"), "Every policy/category needs an explicit zero row");
		}
		check(!new File(categorized+".comparison.summary.tsv").exists(), "Taxonomy alone must not create comparison reports");
		final String compared=invoke(counts, changeText, locusText, summary, null, taxonomyText, "changed,unchanged");
		for(final String suffix:new String[]{".loci.tsv", ".inputs.tsv", ".summary.tsv", ".members.tsv", ".strata.tsv"}){
			check(read(categorized+suffix).equals(read(compared+suffix)), "Comparison must preserve existing reports: "+suffix);
		}
		final String comparison=read(compared+".comparison.summary.tsv");
		check(comparison.contains("changed\tunchanged\tCORPUS_LOWER_QUARTILE\t0\t2\t2\t0\t0\t0\t0\t2\ttrue\tfalse\n")
			&& comparison.contains("changed\tunchanged\tOTHER_CLASSIFIED\t2\t0\t0\t2\t0\t0\t2\t0\tfalse\ttrue\n")
			&& comparison.contains("changed\tunchanged\tUNKNOWN\t2\t2\t1\t1\t0\t0\t1\t1\tfalse\tfalse\n")
			&& comparison.contains("changed\tunchanged\tUNCLASSIFIED\t0\t1\t1\t0\t0\t0\t0\t1\ttrue\tfalse\n")
			&& comparison.endsWith("changed\tunchanged\tALL\t4\t5\t4\t3\t0\t0\t3\t4\tfalse\tfalse\n"),
			"Literal set differences must distinguish net one TP from three exclusive gains and four exclusive baseline detections");
		check(read(compared+".comparison.loci.tsv").contains("1\tp%25α\tCORPUS_LOWER_QUARTILE\tfam\tseq%25α\t10\t30\t+\ttrue\tchanged\tLOST\tfalse\tunchanged\tUNCHANGED\ttrue\n"),
			"Comparison must preserve encoded exact keys and inferred unchanged baseline membership");
		final String reversed=invoke(counts, changeText, locusText, summary, null, taxonomyText, "unchanged,changed");
		check(read(reversed+".comparison.summary.tsv").endsWith("unchanged\tchanged\tALL\t5\t4\t0\t0\t4\t3\t4\t3\tfalse\tfalse\n"), "Reversing policies must reverse exclusive sets");
		final String equal=invoke(counts, changeText, locusText, summary, null, taxonomyText, "baseline,unchanged");
		check(read(equal+".comparison.summary.tsv").endsWith("baseline\tunchanged\tALL\t5\t5\t0\t0\t0\t0\t0\t0\ttrue\ttrue\n"),
			"Changes under an unselected third policy cannot make equal compared sets differ");
		final String losses=changeText.replace("2\tchanged\tGAINED\tp%25α\tseq%25α\t10\t30\t+\n", "")
			.replace("2\tchanged\tGAINED\tp2\tseq%25α\t10\t30\t+\n", "")
			.replace("3\tchanged\tGAINED\tafter\tseq%25α\t10\t30\t+\n", "");
		final String fewer=invoke(counts.replace("2\tchanged\t2\t0\t2", "2\tchanged\t0\t0\t0")
			.replace("3\tchanged\t2\t1\t1", "3\tchanged\t1\t1\t0"), losses, locusText,
			summary.replace("changed\t4\t4\t3", "changed\t1\t4\t0"), null, taxonomyText, "changed,unchanged");
		check(read(fewer+".comparison.summary.tsv").endsWith("changed\tunchanged\tALL\t1\t5\t4\t0\t0\t0\t0\t4\ttrue\tfalse\n"),
			"A true strict subset must be recognized without expanding unchanged scored references");
		invoke(counts, changeText, locusText, summary, "compare requires taxonomy", null, "changed,unchanged");
		invoke(counts, changeText, locusText, summary, "two named policies", taxonomyText, "changed,missing");
		invoke(counts, changeText, locusText, summary, "must be distinct", taxonomyText, "changed,changed");
		invoke(counts, changeText, locusText, summary, "exact FASTA once", taxonomyText.replace(fasta("1"), "wrong.fa"));
		invoke(counts, changeText, locusText, summary, "Missing saved taxonomy input", taxonomyText.replace("1\t"+fasta("1")+"\tCORPUS_LOWER_QUARTILE\n", ""));
		invoke(counts, changeText, locusText, summary, "Unknown saved taxonomy category", taxonomyText.replace("UNCLASSIFIED", "INVENTED"));
		invoke(counts, changeText, locusText, summary, "exact FASTA once", taxonomyText+"1\t"+fasta("1")+"\tCORPUS_LOWER_QUARTILE\n");

		invoke(counts+"1\tchanged\t0\t2\t0\n", changeText, locusText, summary, "Repeated input/policy");
		invoke(counts.replace("3\tunchanged\t2\t0\t0\n", ""), changeText, locusText, summary, "Missing input/policy");
		invoke(counts.replace("3\tchanged\t2\t1\t1", "3\tchanged\t3\t1\t1"), changeText, locusText, summary, "minus lost plus gained");
		invoke(counts, changeText+changeText.substring(changeText.indexOf('\n')+1), locusText, summary, "Repeated changed parent/policy");
		invoke(counts, changeText, locusText.replace("\t10\t30\t", "\t11\t30\t"), summary, "geometry disagrees");
		invoke(counts, changeText, locusText.replace("\tscored\t", "\tdont_care\t"), summary, "scored inventory parent");
		invoke(counts, changeText, locusText.substring(0, locusText.lastIndexOf('\n', locusText.length()-2)+1), summary, "No scored inventory binding");
		invoke(counts, changeText, locusText+locusText.substring(locusText.indexOf('\n')+1), summary, "Duplicate inventory binding");
		invoke(counts, changeText, locusText.replace("major", "unknown"), summary, "undeclared member");
		invoke(counts, changeText, locusText, summary.replace("changed\t4\t4\t3", "changed\t4\t4\t2"), "disagree with policy summary");
		invoke(counts, changeText, locusText, summary+"changed\t4\t4\t3\n", "Repeated policy summary");
		invoke(counts, changeText, locusText.replace("\tscored\tmajor,minor%25α\t", "\tscored\t.\t"), summary, "requires qualified member support");
		invoke(counts, changeText, locusText.replace("\tscored\tmajor,minor%25α\t", "\tscored\t\t"), summary, "requires qualified member support");
		invoke(counts.replace("1\tunchanged\t2\t0\t0", "1\tunchanged\t3\t0\t1"),
			changeText+"1\tunchanged\tGAINED\tp%25α\tseq%25α\t10\t30\t+\n", locusText,
			summary.replace("unchanged\t5\t0\t0", "unchanged\t6\t0\t1"), "both LOST and GAINED");
	}

	/** Appends matching literal geometry to the two independent metadata schemas. */
	private void locus(final ByteBuilder changes, final ByteBuilder loci, final String id, final String parent,
		final String direction, final String scored, final String marginal){
		changes.append(id).append("\tchanged\t").append(direction).tab();
		field(changes, parent).tab();
		field(changes, "seq%α").append("\t10\t30\t+\n");
		loci.append(id).tab();
		field(loci, parent).append("\tfam\t");
		field(loci, "seq%α").append("\t10\t30\t+\tscored\t");
		field(loci, scored).tab();
		field(loci, marginal).nl();
	}

	/** A negative fixture passes only when the expected diagnostic is observed. */
	private String invoke(final String counts, final String changes, final String loci, final String summary, final String expected){
		return invoke(counts, changes, loci, summary, expected, null);
	}

	/** Optional category metadata must validate before any output is opened. */
	private String invoke(final String counts, final String changes, final String loci, final String summary,
			final String expected, final String taxonomy){
		return invoke(counts, changes, loci, summary, expected, taxonomy, null);
	}

	/** Optional set comparison shares the same native counts and literal key oracle. */
	private String invoke(final String counts, final String changes, final String loci, final String summary,
			final String expected, final String taxonomy, final String comparison){
		final String stem="case"+(cases++);
		write(stem+".counts", counts);
		write(stem+".changes", changes);
		write(stem+".loci", loci);
		write(stem+".summary", summary);
		if(taxonomy!=null){write(stem+".taxonomy", taxonomy);}
		final String prefix=path(stem+".report");
		try{
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList(
				"in="+path(stem+".counts"), "changes="+path(stem+".changes"),
				"loci="+path(stem+".loci"), "summary="+path(stem+".summary"), "reference="+path("reference.tsv"),
				"members="+path("members.tsv"), "family=fam", "policies=changed,unchanged", "highlight=minor%α", "out="+prefix));
			if(taxonomy!=null){args.add("taxonomy="+path(stem+".taxonomy"));}
			if(comparison!=null){args.add("compare="+comparison);}
			SummarizeNcrnaSweep.main(args.toArray(new String[args.size()]));
		}catch(IllegalArgumentException failure){
			check(expected!=null && failure.getMessage().contains(expected), "Unexpected fixture failure: "+failure);
			check(!new File(prefix+".summary.tsv").exists(), "Invalid input must fail before reports are opened");
			return prefix;
		}
		check(expected==null, "Expected a failure containing: "+expected);
		return prefix;
	}

	/*--------------------------------------------------------------*/
	/*----------------           Helpers            ----------------*/
	/*--------------------------------------------------------------*/

	private void write(final String name, final String text){
		try(Output out=new Output(path(name))){out.writer.print(text.getBytes(StandardCharsets.UTF_8));}
	}
	private static String read(final String path){
		final ByteBuilder buffer=new ByteBuilder();
		try(Input in=new Input(path)){
			for(byte[] line; (line=in.reader.nextLine())!=null;){buffer.append(line).nl();}
		}
		return utf8(buffer);
	}
	private static String utf8(final ByteBuilder buffer){return new String(buffer.array, 0, buffer.length(), StandardCharsets.UTF_8);}
	private String path(final String name){return new File(root, name).toString();}
	private String fasta(final String id){return path("missing"+id+".fa");}
	private static ByteBuilder field(final ByteBuilder buffer, final String text){return BuildNcrnaReference.field(buffer, text);}
	private static void check(final boolean valid, final String message){
		if(!valid){throw new IllegalStateException(message);}
	}
	private final File root;
	private int cases;
	private static final String[] IDS={"1", "2", "3", "999999999999999999999"};
}

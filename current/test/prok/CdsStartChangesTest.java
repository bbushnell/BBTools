package test.prok;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.GradeCdsStartChanges;
import prok.GradeCdsStartChanges.Reference;
import prok.GradeCdsStartChanges.ReferenceContig;
import prok.NcrnaFamilyRegistry;
import prok.NcrnaCdsStartChanges;
import prok.NcrnaGffFamilyGrader.Locus;
import prok.NcrnaGffFamilyGrader;
import prok.NcrnaGffFamilyGrader.MatchResult;
import prok.PanelAnnotationGrade;
import prok.PanelAnnotationGrade.Feature;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.IntList;

/** Literal paired annotations: start correction on both strands, wrong shortening,
 * lost/added CDS, origin-spanning truth, uncertain starts, duplicates and an empty
 * input. Reused contig names and IDs must remain scoped to their exact input.
 * @author Brian, Raiden
 */
public final class CdsStartChangesTest{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, CdsStartChangesTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			String testcase="normal";
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("testcase")){testcase=b;}
				else{check(parser.parse(arg, a, b), "Use standard fixture out=/thread/JVM parameters");}
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"The fixture needs a fresh output directory and never reads external genomes");
			final String root=new File(parser.out1).getAbsolutePath();
			check(new File(root).mkdir(), "Preserve earlier fixture results rather than overwriting them");
			if(testcase.equals("virtual_origin")){
				virtualOrigin(root);
				pp.outstream.println("CDS_VIRTUAL_ORIGIN_FIXTURES_PASS");
				return;
			}
			check(testcase.equals("normal") || testcase.equals("unknown_contig") || testcase.equals("missing_length")
				|| testcase.equals("declaration_conflict") || testcase.equals("repeated_region_conflict")
				|| testcase.equals("taxid_alias") || testcase.equals("taxid_alias_conflict"), "Use a declared finite fixture case");
			fixture(root, testcase);
			if(testcase.startsWith("taxid_alias")){taxidAliases(root, testcase.equals("taxid_alias_conflict"));}
			GradeCdsStartChanges.main(new String[]{"in="+root+"/before.tsv", "in2="+root+"/after.tsv",
				"corpus="+root+"/corpus.tsv", "expectedinputs=5", "out="+root+"/result", "t=2"});
			check(testcase.equals("normal") || testcase.equals("taxid_alias"), "Expected validation failure for "+testcase);
			verify(root+"/result");
			referenceViews(root, testcase.equals("taxid_alias"));
			repeatedStartChanges(root);
			if(testcase.equals("normal")){
				reportedRna(root);
				virtualOrigin(root);
			}
			pp.outstream.println("CDS_START_CHANGES_FIXTURES_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/*--------------------------------------------------------------*/
	/*----------------         Literal Inputs       ----------------*/
	/*--------------------------------------------------------------*/

	private static void fixture(final String root, final String testcase){
		try(Output before=new Output(root+"/before.tsv"); Output after=new Output(root+"/after.tsv");
			Output corpus=new Output(root+"/corpus.tsv")){
			before.writer.println("input_id\tfasta\tgff");
			after.writer.println("input_id\tfasta\tgff");
			corpus.writer.println(CORPUS);
			for(final String id:new String[]{"A", "B", "U", "E", "D"}){
				before.writer.println(id+"\t"+root+"/"+id+".fa\t"+root+"/"+id+".before.gff");
				after.writer.println(id+"\t"+root+"/"+id+".fa\t"+root+"/"+id+".after.gff");
				corpus.writer.println("assembly_"+id+"\t"+id+"\tfixture\t"+root+"/"+id+".fa\t/new_"+id+
					".gff\t\t\t"+root+"/"+id+".fna.gz\t/new_"+id+".gff\tfasta_gzip\tncrna\tcore");
			}
		}
		write(root+"/A.gff.gz", HEADER, testcase.equals("missing_length") ? "# no declared length" : "##sequence-region ctg 1 1000",
			region(false, testcase.equals("declaration_conflict") ? 999 : 1000),
			testcase.equals("repeated_region_conflict") ? region(false, 999) : "# no conflicting repeated region",
			ref(201, 500, '+', "ID=plus"), ref(600, 899, '-', "ID=minus"), ref(50, 90, '+', "ID=lost"));
		final String clan=NcrnaFamilyRegistry.displayName("cobalamin_rf00174");
		write(root+"/A.before.gff", HEADER, call(101, 500, '+'), call(600, 999, '-'), call(50, 90, '+'),
			rna(10, 40, "RNA", "RNA,model:cobalamin_fixture,family:"+clan),
			rna(920, 960, "RNA", "RNA,model:rnasep_fixture,family:RNaseP"));
		write(root+"/A.after.gff", HEADER, call(201, 500, '+'), call(600, 899, '-'), call(700, 760, '+'),
			rna(10, 40, "RNA", "RNA,model:cobalamin_fixture,family:"+clan),
			rna(60, 80, "RNA", "RNA,model:cobalamin_fixture,family:"+clan),
			rna(920, 960, "RNA", "RNA,model:rnasep_fixture,family:RNaseP"),
			testcase.equals("unknown_contig") ? "ghost\tBBTools\tCDS\t10\t100\t.\t+\t0\tCDS,startScr:1" : "# no undeclared query");

		// Circumference1000 is not divisible by3. The wrapped reference remains ONE row.
		write(root+"/B.gff.gz", HEADER, ref(951, 1100, '+', "ID=wrapped"), ref(400, 699, '+', "ID=plus"),
			region(true, 1000), "##sequence-region ctg 1 1000");
		write(root+"/B.before.gff", HEADER, call(1, 100, '+'), call(400, 699, '+'), rna(150, 200, "rRNA", "16S,startScr:0"));
		write(root+"/B.after.gff", HEADER, call(1, 100, '+'), call(499, 699, '+'), rna(150, 200, "rRNA", "16S,startScr:0"),
			rna(250, 300, "tRNA", "tRNA,startScr:0"));

		write(root+"/U.gff.gz", HEADER, "##sequence-region ctg 1 2000", region(false, 2000),
			ref(201, 500, '+', "ID=partial;partial=true"), ref(600, 899, '+', "ID=joined"),
			ref(950, 999, '+', "ID=joined"), ref(1100, 1399, '+', "ID=ambiguous_a"), ref(1199, 1399, '+', "ID=ambiguous_b"));
		write(root+"/U.before.gff", HEADER, call(101, 500, '+'), call(600, 899, '+'), call(950, 999, '+'),
			call(1100, 1399, '+'), call(1400, 1500, '+'));
		write(root+"/U.after.gff", HEADER, call(201, 500, '+'), call(699, 899, '+'), call(950, 999, '+'),
			call(1199, 1399, '+'), call(1410, 1500, '+'));

		write(root+"/E.gff.gz", HEADER, "##sequence-region ctg 1 1000", region(false, 1000));
		write(root+"/E.before.gff", HEADER);
		write(root+"/E.after.gff", HEADER);
		write(root+"/D.gff.gz", HEADER, "##sequence-region ctg 1 2000", region(false, 2000),
			ref(201, 500, '+', "ID=duplicate_truth"), ref(600, 899, '+', "ID=unique_to_multiple"),
			ref(1100, 1399, '+', "ID=multiple_both"));
		write(root+"/D.before.gff", HEADER, call(201, 500, '+'), call(201, 500, '+'), call(600, 899, '+'),
			call(1100, 1399, '+'), call(1100, 1399, '+'));
		write(root+"/D.after.gff", HEADER, call(201, 500, '+'), call(600, 899, '+'), call(699, 899, '+'),
			call(1100, 1399, '+'), call(1199, 1399, '+'));
	}

	/*--------------------------------------------------------------*/
	/*----------------          Literal Oracles      ----------------*/
	/*--------------------------------------------------------------*/

	private static void verify(final String prefix){
		final LineParser1 p=new LineParser1('\t');
		int inputRows=0;
		try(Input input=new Input(prefix+".inputs.tsv")){
			input.reader.nextLine();
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==16, "Each input row retains all13native metrics and its qualification");
				final String id=p.parseString(0);
				final boolean before=p.termEquals("before", 1);
				if(id.equals("A")){
					check(p.parseLong(2)==(before ? 3 : 2) && p.parseLong(3)==(before ? 0 : 1)
						&& p.parseLong(4)==(before ? 0 : 1) && p.parseLong(5)==(before ? 1 : 2),
						"Both corrected strands improve starts while the lost gene and added unmatched CDS remain visible");
				}else if(id.equals("B")){
					check(p.parseLong(2)==2 && p.parseLong(3)==0 && p.parseLong(4)==0
						&& p.parseLong(5)==(before ? 1 : 0) && p.parseLong(6)==2,
						"Wrapped truth is matched once with modular3prime exactness, while a wrong shortened start is still a TP");
				}else if(id.equals("U")){
					check(p.parseLong(2)==4 && p.parseLong(3)==1 && p.parseLong(4)==1,
						"Uncertain starts remain in whole-source feature-row detection grading");
				}else if(id.equals("E")){
					for(int i=2; i<15; i++){check(p.parseLong(i)==0, "The empty input retains an explicit zero row for each arm");}
				}else{
					check(id.equals("D") && p.parseLong(2)==3 && p.parseLong(3)==2
						&& p.parseLong(8)==(before ? 2 : 0), "Duplicate native calls and alternative starts remain unmatched extras, never deduplicated before matching");
				}
				inputRows++;
			}
		}
		check(inputRows==10, "Five complete input pairs must survive ordered multiworker output");
		int changes=0, uncertain=0, ambiguous=0, unknown=0, duplicate=0, ambiguousBefore=0, ambiguousAfter=0;
		boolean sawNew899=false, sawNew1399=false;
		try(Input input=new Input(prefix+".changes.tsv")){
			input.reader.nextLine();
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				check(p.terms()==14, "Changed and lost calls retain coordinates, uncertainty and scope");
				if(p.termEquals("UNCERTAIN_PARTIAL_OR_REPEATED_ID", 9)){
					check(p.termEquals("NA", 10) && p.termEquals("NA", 11), "Uncertain annotation cannot become numeric start truth");
					uncertain++;
				}
				if(p.termEquals("AMBIGUOUS_ANNOTATED_START", 9)){ambiguous++;}
				if(p.termEquals("NO_SAME_STOP_ANNOTATION", 9)){unknown++;}
				if(p.termEquals("AMBIGUOUS_BEFORE_STOP_MEMBER", 3)){
					check(!p.termEquals("NA", 4) && p.termEquals("NA", 6), "Each ambiguous baseline member retains only its own coordinates");
					ambiguousBefore++; duplicate++;
				}
				if(p.termEquals("AMBIGUOUS_AFTER_STOP_MEMBER", 3)){
					check(p.termEquals("NA", 4) && !p.termEquals("NA", 6), "Each ambiguous candidate member retains only its own coordinates");
					sawNew899|=p.parseInt(6)==699 && p.parseInt(7)==899;
					sawNew1399|=p.parseInt(6)==1199 && p.parseInt(7)==1399;
					ambiguousAfter++; duplicate++;
				}
				if(p.termEquals("B", 0)){
					check(p.termEquals("UNIQUE_COMPLETE_ANNOTATED_START", 9) && p.parseInt(10)==0 && p.parseInt(11)==99,
						"A reused ID on another input is not multipart; wrong circular start shortening has positive downstream error");
				}
				changes++;
			}
		}
		check(changes==19 && uncertain==2 && ambiguous==1 && unknown==2 && duplicate==10,
			"Every replacement/lost/added/ambiguous call remains in a distinct, honest detailed row");
		check(ambiguousBefore==5 && ambiguousAfter==5 && sawNew899 && sawNew1399,
			"Unique-to-multiple and multiple-both-side groups retain ALL candidate coordinates without arbitrary pairing");
		metric(prefix, "paired", "reliable_changed_start_pairs", 3);
		metric(prefix, "paired", "same_stop_changed_start_exact_before", 1);
		metric(prefix, "paired", "same_stop_changed_start_exact_after", 2);
		metric(prefix, "paired", "changed_start_error_worsened", 1);
		metric(prefix, "paired", "shortened_past_annotated_start", 1);
		metric(prefix, "paired", "annotated_CDS_detection_lost", 1);
		metric(prefix, "paired", "other_RNA_family_count_changes", 1);
		metric(prefix, "before", "TP", 12);
		metric(prefix, "after", "TP", 11);
	}

	private static void metric(final String prefix, final String arm, final String name, final long expected){
		final LineParser1 p=new LineParser1('\t');
		boolean found=false;
		try(Input input=new Input(prefix+".summary.tsv")){
			input.reader.nextLine();
			for(byte[] row; (row=input.reader.nextLine())!=null;){
				p.set(row);
				if(p.termEquals(arm, 0) && p.termEquals(name, 1)){
					check(!found && p.parseLong(2)==expected, "Literal expected metric "+arm+"/"+name+"="+expected);
					found=true;
				}
			}
		}
		check(found, "Required literal metric was not emitted: "+name);
	}
	/** Opt-in accounting preserves all five default tables, then exercises exact
	 * losses concealed by equal totals, two physical contigs, and an empty input.
	 */
	private static void reportedRna(final String root){
		final String target=NcrnaFamilyRegistry.displayName("cobalamin_rf00174");
		write(root+"/rna_families.tsv", "family\tmodel_prefix\tdisplay_label\tinclude",
			"cobalamin\tcobalamin\t"+target+"\ttrue", "RNaseP\trnasep\tRNaseP\tfalse", "FMN\tFMN\tFMN\tfalse",
			"SRP\tsrp_small\tSRP RNA (small)\tfalse", "SRP\tsrp_large\tSRP RNA (large)\tfalse",
			"6S\t6S_RF00013\t6S\tfalse", "6S\t6S_RF01685\t6S-Flavo\tfalse");
		GradeCdsStartChanges.main(new String[]{"in="+root+"/before.tsv", "in2="+root+"/after.tsv",
			"corpus="+root+"/corpus.tsv", "expectedinputs=5", "out="+root+"/exact", "rnacallchanges=t",
			"familymap="+root+"/rna_families.tsv", "t=2"});
		for(final String suffix:new String[]{".summary.tsv", ".inputs.tsv", ".changes.tsv", ".rna_counts.tsv", ".sources.tsv"}){
			equalFiles(root+"/result"+suffix, root+"/exact"+suffix);
		}
		check(!new File(root+"/result.rna_call_changes.tsv").exists(), "The historical default must not emit the opt-in sidecar");
		final String fixture=root+"/rna";
		check(new File(fixture).mkdir(), "Use a fresh, finite two-input saved-RNA fixture");
		try(Output before=new Output(fixture+"/before.tsv"); Output after=new Output(fixture+"/after.tsv"); Output corpus=new Output(fixture+"/corpus.tsv")){
			before.writer.println("input_id\tfasta\tgff");
			after.writer.println("input_id\tfasta\tgff");
			corpus.writer.println(CORPUS);
			for(final String id:new String[]{"X", "E"}){
				before.writer.println(id+"\t"+fixture+"/"+id+".fa\t"+fixture+"/"+id+".before.gff");
				after.writer.println(id+"\t"+fixture+"/"+id+".fa\t"+fixture+"/"+id+".after.gff");
				corpus.writer.println("assembly_"+id+"\t"+id+"\tfixture\t"+fixture+"/"+id+".fa\t/unused.gff\t\t\t"
					+fixture+"/"+id+".fna.gz\t/unused.gff\tfasta_gzip\tncrna\tcore");
			}
		}
		write(fixture+"/X.gff.gz", HEADER, "##sequence-region ctg 1 1000", "##sequence-region other 1 1000");
		write(fixture+"/E.gff.gz", HEADER, "##sequence-region ctg 1 1000");
		write(fixture+"/E.before.gff", HEADER);
		write(fixture+"/E.after.gff", HEADER);
		write(fixture+"/X.before.gff", HEADER,
			rna(10, 20, "tRNA", "tRNA,startScr:0"), rna(10, 20, "tRNA", "tRNA,startScr:0"),
			rna(30, 40, "rRNA", "16S,startScr:0"), rna(50, 60, "RNA", "RNA,family:RNaseP,model:rnasep_old"),
			rna(70, 80, "RNA", "RNA,family:FMN,model:FMN_m"), rna(110, 120, "RNA", "RNA,family:RNaseP,model:rnasep_old"),
			rna(90, 100, "RNA", "RNA,family:"+target+",model:cobalamin_old"),
			rna(300, 350, "RNA", "RNA,family:SRP RNA (small),model:srp_small_old"),
			rna(400, 450, "RNA", "RNA,family:6S,model:6S_RF00013_old"),
			"other\tBBTools\ttRNA\t10\t20\t.\t+\t0\ttRNA,startScr:0");
		write(fixture+"/X.after.gff", HEADER,
			rna(10, 20, "tRNA", "tRNA,startScr:0"), rna(30, 40, "rRNA", "23S,startScr:0"),
			"ctg\tBBTools\tRNA\t50\t60\t.\t-\t0\tRNA,family:RNaseP,model:rnasep_new",
			rna(71, 80, "RNA", "RNA,family:FMN,model:FMN_m"), rna(110, 120, "RNA", "RNA,family:RNaseP,model:rnasep_new"),
			rna(80, 90, "RNA", "RNA,family:"+target+",model:cobalamin_replacement"),
			rna(300, 350, "RNA", "RNA,family:SRP RNA (large),model:srp_large_new"),
			rna(400, 450, "RNA", "RNA,family:6S-Flavo,model:6S_RF01685_new"),
			"other\tBBTools\ttRNA\t10\t20\t.\t+\t0\ttRNA,startScr:0",
			"other\tBBTools\ttRNA\t200\t210\t.\t+\t0\ttRNA,startScr:0");
		GradeCdsStartChanges.main(new String[]{"in="+fixture+"/before.tsv", "in2="+fixture+"/after.tsv",
			"corpus="+fixture+"/corpus.tsv", "expectedinputs=2", "out="+fixture+"/result", "rnacallchanges=t",
			"familymap="+root+"/rna_families.tsv", "t=2"});
		try(Input input=new Input(fixture+"/result.rna_call_changes.tsv")){
			check(Arrays.equals(input.reader.nextLine(), "input_id\tbefore_calls\tafter_calls\tlost_calls\tgained_calls\tqualification".getBytes(java.nio.charset.StandardCharsets.UTF_8)),
				"The exact saved-RNA sidecar has a declared complete schema");
			final LineParser1 p=new LineParser1('\t');
			p.set(input.reader.nextLine());
			check(p.terms()==6 && p.termEquals("X", 0) && p.parseLong(1)==9 && p.parseLong(2)==9
				&& p.parseLong(3)==4 && p.parseLong(4)==4,
				"Nine calls conceal four losses; same-scoring-family SRP/6S subtype switches retain both reports like training/coarse");
			p.set(input.reader.nextLine());
			check(p.terms()==6 && p.termEquals("E", 0) && p.parseLong(1)==0 && p.parseLong(2)==0
				&& p.parseLong(3)==0 && p.parseLong(4)==0 && input.reader.nextLine()==null,
				"Zero-call inputs retain an explicit ordered row; target-family movement is excluded");
		}
	}
	private static void equalFiles(final String before, final String after){
		try(Input old=new Input(before); Input next=new Input(after)){
			for(byte[] row; (row=old.reader.nextLine())!=null;){
				check(Arrays.equals(row, next.reader.nextLine()), "Opting into exact RNA accounting must preserve default table bytes: "+before);
			}
			check(next.reader.nextLine()==null, "Opt-in default parity cannot append rows: "+after);
		}
	}
	/** Literal source-only API checks reuse the same original GFF reader. */
	private static void referenceViews(final String root, final boolean qualified){
		final Reference a=new Reference("A", root+"/A.gff.gz"), b=new Reference("B", root+"/B.gff.gz");
		final ReferenceContig ac=a.contig("ctg description", 1000), bc=b.contig("ctg", 1000);
		check(a.input.equals("A") && a.originalGff.equals(root+"/A.gff.gz") && a.contigs()==1 && a.cdsRows()==3
			&& ac.input().equals("A") && ac.name().equals("ctg") && ac.length()==1000
			&& ac.circularLength()==0 && ac.rows()==3 && a.contig(0).rows()==3,
			"A source reference retains input identity, declared physical geometry and all original-order CDS rows without queries");
		check(ac.feature(0).locus.start==201 && ac.feature(0).locus.stop==500 && ac.completeUniqueRow(0)
			&& bc.circularLength()==1000 && bc.feature(0).locus.start==951 && bc.feature(0).locus.stop==1100
			&& bc.completeUniqueRow(1), "Wrapped CDS stays one unrolled row; IDs reused in another input do not become multipart");
		final ReferenceContig u=new Reference("U", root+"/U.gff.gz").contig(0);
		check(!u.completeUniqueRow(0) && !u.completeUniqueRow(1) && !u.completeUniqueRow(2)
			&& u.completeUniqueRow(3) && u.completeUniqueRow(4)
			&& u.feature(3).locus.stop==u.feature(4).locus.stop
			&& u.feature(3).locus.start!=u.feature(4).locus.start,
			"Row-level partial/repeated-ID eligibility must remain separate from grouped same-stop start ambiguity");
		final Reference empty=new Reference("E", root+"/E.gff.gz");
		check(empty.contigs()==1 && empty.cdsRows()==0 && empty.contig(0).rows()==0,
			"A declared zero-CDS contig remains present in the source reference");
		referenceFailure(()->a.contig("ghost", 1000), "declared reference identity");
		referenceFailure(()->a.contig("ctg", 999), "length");
		if(qualified){
			check(a.contig("tid|42|ctg description", 1000).name().equals("ctg"), "Only the declared taxID alias resolves to the physical reference contig");
		}else{referenceFailure(()->a.contig("tid|42|ctg", 1000), "qualified identity");}
		referenceFailure(()->a.contig("tid|43|ctg", 1000), "qualified identity");
		referenceFailure(()->b.contig("tid|42|ctg", 1000), "qualified identity");
		final String source=root+"/reference_only.gff";
		write(source, HEADER, "##sequence-region seen 1 500", "##sequence-region reference_only 1 500",
			"seen\tNCBI\tCDS\t10\t99\t.\t+\t0\tID=seen",
			"reference_only\tNCBI\tCDS\t100\t399\t.\t-\t0\tID=only");
		final Reference only=new Reference("only", source);
		check(only.contigs()==2 && only.cdsRows()==2 && only.contig(1).name().equals("reference_only")
			&& only.contig(1).rows()==1, "Enumeration conserves reference-only contigs so a path consumer can retain their false negatives");
	}
	private static void referenceFailure(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException failure){
			check(failure.getMessage().contains(fragment), "Unexpected source-reference failure: "+failure);
			return;
		}
		throw new AssertionError("Expected source-reference failure containing "+fragment);
	}
	/** Literal geometry from the circular Pandoraea transposase parts that blocked
	 * saved-CDS reporting. Expectations come from inclusive physical segments:
	 * part1=416..665; part2=126554..126976 plus1..416. Never merge their rows.
	 */
	private static void virtualOrigin(final String root){
		final int circle=126976;
		final String path=root+"/virtual_origin.gff";
		final String part1="ctg\tRefSeq\tCDS\t127392\t127641\t.\t-\t0\tID=joined;exception=ribosomal slippage;part=1";
		final String part2="ctg\tRefSeq\tCDS\t126554\t127392\t.\t-\t2\tID=joined;exception=ribosomal slippage;part=2";
		write(path, HEADER, "##sequence-region ctg 1 126976", region(true, circle), part1, part2);
		final Reference reference=new Reference("virtual", path);
		final ReferenceContig contig=reference.contig("ctg", circle);
		check(reference.cdsRows()==2 && contig.rows()==2 && contig.circularLength()==circle
			&& contig.feature(0).locus.start==127392 && contig.feature(0).locus.stop==127641
			&& contig.feature(1).locus.start==126554 && contig.feature(1).locus.stop==127392
			&& !contig.completeUniqueRow(0) && !contig.completeUniqueRow(1),
			"Keep both original multipart rows in order and exclude their repeated IDs from reliable start truth");
		final int[][] cases={{127392,127641,416,665,250}, {126554,127392,126554,127392,839},
			{126554,127392,126554,126976,423}, {126554,127392,1,416,416},
			{127392,127641,126554,127392,1}, {127392,127641,666,700,0},
			{254368,254617,416,665,250}, {127392,254367,1,126976,126976}};
		for(final char strand:new char[]{'+', '-'}){
			for(final int[] row:cases){
				final Locus a=new Locus("ctg", "CDS", row[0], row[1], strand);
				final Locus b=new Locus("ctg", "CDS", row[2], row[3], strand);
				check(NcrnaGffFamilyGrader.circularOverlap(a, b, circle)==row[4]
					&& NcrnaGffFamilyGrader.circularOverlap(b, a, circle)==row[4],
					"Circular overlap must match literal inclusive segments, symmetrically on both strands: "+Arrays.toString(row));
				check(a.start==row[0] && a.stop==row[1] && b.start==row[2] && b.stop==row[3],
					"Overlap normalization borrows immutable loci; it never rewrites source coordinates");
			}
			final ArrayList<Feature> truth=new ArrayList<Feature>(), calls=new ArrayList<Feature>();
			truth.add(new Feature(new Locus("ctg", "CDS", 127392, 127641, strand), 2, 4));
			calls.add(new Feature(new Locus("ctg", "CDS", 416, 665, strand), 2, 9));
			final MatchResult match=PanelAnnotationGrade.match(calls, truth, circle);
			final long[] stats=PanelAnnotationGrade.statistics(calls, truth, match, true, circle);
			check(match.matched==1 && match.callMatchedTruthIdx[0]==0 && match.truthMatchedCallIdx[0]==0
				&& stats[3]==1 && stats[4]==1 && stats[5]==1 && stats[8]==1 && stats[9]==0,
				"Whole-turn copies match once with both endpoints exact and phase preserved even when circumference is not divisible by3");
		}
		final ArrayList<Feature> calls=new ArrayList<Feature>(), truth=new ArrayList<Feature>();
		truth.add(contig.feature(0));
		truth.add(contig.feature(1));
		calls.add(cds(126554, 126976, '-'));
		calls.add(cds(416, 665, '-'));
		final MatchResult match=PanelAnnotationGrade.match(calls, truth, circle);
		check(match.matched==2 && match.callMatchedTruthIdx[0]==1 && match.callMatchedTruthIdx[1]==0,
			"Virtual parts retain original truth indices through one-to-one matching, not their translated coordinate order");
		write(root+"/virtual_linear.gff", HEADER, "##sequence-region ctg 1 126976", region(false, circle), part1);
		referenceFailure(()->new Reference("linear", root+"/virtual_linear.gff"), "declared circular geometry");
		write(root+"/virtual_unknown.gff", HEADER, "##sequence-region ctg 1 126976", part1);
		referenceFailure(()->new Reference("unknown", root+"/virtual_unknown.gff"), "declared circular geometry");
		write(root+"/virtual_long.gff", HEADER, "##sequence-region ctg 1 126976", region(true, circle),
			ref(1, 126977, '+', "ID=too_long"));
		referenceFailure(()->new Reference("long", root+"/virtual_long.gff"), "at most one circumference");
		referenceFailure(()->NcrnaGffFamilyGrader.circularOverlap(new Locus("ctg", "CDS", 1, 126977, '+'),
			new Locus("ctg", "CDS", 1, 10, '+'), circle), "at most one turn");
	}
	/** Reuses the literal five-input annotations as an independent saved-tool
	 * oracle; new uncertainty counters describe only unique changed-start pairs.
	 */
	private static void repeatedStartChanges(final String root){
		final long[] totals=new long[NcrnaCdsStartChanges.statisticsCount()];
		add(totals, startTrial(root, "A", new Feature[]{cds(101, 500, '+'), cds(600, 999, '-'), cds(50, 90, '+')},
			new Feature[]{cds(201, 500, '+'), cds(600, 899, '-'), cds(700, 760, '+')},
			new long[]{2, 1, 1, 2, 0, 2, 0, 0, 0, 0, 0, 0}));
		add(totals, startTrial(root, "B", new Feature[]{cds(1, 100, '+'), cds(400, 699, '+')},
			new Feature[]{cds(1, 100, '+'), cds(499, 699, '+')},
			new long[]{1, 0, 0, 1, 1, 0, 1, 1, 0, 0, 0, 0}));
		add(totals, startTrial(root, "U", new Feature[]{cds(101, 500, '+'), cds(600, 899, '+'), cds(950, 999, '+'),
			cds(1100, 1399, '+'), cds(1400, 1500, '+')},
			new Feature[]{cds(201, 500, '+'), cds(699, 899, '+'), cds(950, 999, '+'), cds(1199, 1399, '+'), cds(1410, 1500, '+')},
			new long[]{4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2}));
		add(totals, startTrial(root, "D", new Feature[]{cds(201, 500, '+'), cds(201, 500, '+'), cds(600, 899, '+'),
			cds(1100, 1399, '+'), cds(1100, 1399, '+')},
			new Feature[]{cds(201, 500, '+'), cds(600, 899, '+'), cds(699, 899, '+'), cds(1100, 1399, '+'), cds(1199, 1399, '+')},
			new long[]{0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 0}));
		add(totals, startTrial(root, "E", new Feature[0], new Feature[0], new long[12]));
		check(Arrays.equals(totals, new long[]{7, 1, 1, 3, 1, 2, 1, 1, 10, 1, 1, 2}),
			"Repeated start statistics must retain all reliable, uncertain, ambiguous and absent-reference outcomes");
		for(int i=0; i<9; i++){metric(root+"/result", "paired", NcrnaCdsStartChanges.name(i), totals[i]);}
	}
	private static long[] startTrial(final String root, final String input, final Feature[] before, final Feature[] after,
			final long[] expected){
		final Feature[] pool=new Feature[before.length+after.length];
		System.arraycopy(before, 0, pool, 0, before.length);
		System.arraycopy(after, 0, pool, before.length, after.length);
		final Reference reference=new Reference(input, root+"/"+input+".gff.gz");
		final NcrnaCdsStartChanges grade=new NcrnaCdsStartChanges(pool, reference.contig(0));
		final IntList old=new IntList(), next=new IntList();
		for(int i=0; i<before.length; i++){old.add(i);}
		for(int i=0; i<after.length; i++){next.add(before.length+i);}
		grade.grade(old);
		grade.freezeBaseline("current multiplier");
		final long[] borrowed=grade.grade(next);
		check(Arrays.equals(borrowed, expected), "Literal same-stop start statistics differ for "+input+": "+Arrays.toString(borrowed));
		grade.grade(old);
		for(int i=0; i<borrowed.length; i++){
			if(i!=8){check(borrowed[i]==0, "Returning to baseline must undo changed-start metrics, not retain the preceding trial");}
		}
		check(grade.grade(next)==borrowed && Arrays.equals(borrowed, expected), "Start trials reuse their result array and fixed baseline");
		final IntList outside=new IntList();
		outside.add(pool.length);
		boolean rejected=false;
		try{grade.grade(outside);}
		catch(IllegalArgumentException failure){rejected=failure.getMessage().contains("inside its fixed pool");}
		check(rejected, "An out-of-pool CDS start candidate must fail");
		boolean stale=false;
		try{grade.statistics();}
		catch(IllegalArgumentException failure){stale=failure.getMessage().contains("completed selected trial");}
		check(stale, "A failed start trial invalidates borrowed results");
		check(Arrays.equals(grade.grade(next), expected), "A failed trial cannot replace or corrupt the frozen start baseline");
		return grade.statistics().clone();
	}
	private static Feature cds(final int start, final int stop, final char strand){
		return new Feature(new Locus("ctg", "CDS", start, stop, strand), 0, start);
	}
	private static void add(final long[] total, final long[] values){
		check(total.length==values.length, "Start statistics must retain one shared schema");
		for(int i=0; i<total.length; i++){total[i]+=values[i];}
	}
	private static void write(final String path, final String... rows){try(Output out=new Output(path)){for(final String row:rows){out.writer.println(row);}}}
	private static String region(final boolean circular, final int length){return "ctg\tNCBI\tregion\t1\t"+length+"\t.\t+\t.\tID=region;Is_circular="+circular;}
	private static String ref(final int start, final int stop, final char strand, final String attrs){return "ctg\tNCBI\tCDS\t"+start+"\t"+stop+"\t.\t"+strand+"\t0\t"+attrs;}
	private static String call(final int start, final int stop, final char strand){return "ctg\tBBTools\tCDS\t"+start+"\t"+stop+"\t.\t"+strand+"\t0\tCDS,startScr:1";}
	private static String rna(final int start, final int stop, final String type, final String attrs){return "ctg\tBBTools\t"+type+"\t"+start+"\t"+stop+"\t.\t+\t0\t"+attrs;}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String HEADER="##gff-version 3";
	private static final String CORPUS="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
	/** Preserve bare sequence-region names but qualify all A feature/query rows,
	 * matching the actual older-source failure. A second taxID alias must fail. */
	private static void taxidAliases(final String root, final boolean conflict){
		for(final String suffix:new String[]{".gff.gz", ".before.gff", ".after.gff"}){
			final String path=root+"/A"+suffix;
			final java.util.ArrayList<byte[]> rows=new java.util.ArrayList<byte[]>();
			try(Input input=new Input(path)){
				for(byte[] row; (row=input.reader.nextLine())!=null;){rows.add(row);}
			}
			final String qualified=path+".qualified";
			try(Output output=new Output(qualified)){
				for(final byte[] row:rows){
					if(row.length>0 && row[0]!='#'){
						output.writer.print(new structures.ByteBuilder("tid|42|").append(row).nl());
					}else{output.writer.println(row);}
				}
				if(conflict && suffix.equals(".gff.gz")){
					output.writer.println("tid|43|ctg\tNCBI\tgene\t1\t10\t.\t+\t.\tID=conflicting_alias");
				}
			}
			check(new File(path).renameTo(new File(path+".bare")), "Preserve the original finite fixture bytes");
			// The input reader selects decompression by suffix, so rewrite through
			// its native writer rather than renaming plain bytes to a gzip suffix.
			try(Input input=new Input(qualified); Output output=new Output(path)){
				for(byte[] row; (row=input.reader.nextLine())!=null;){output.writer.println(row);}
			}
		}
	}
}

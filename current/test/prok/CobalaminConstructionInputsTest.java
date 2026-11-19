package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.FileFormat;
import gff.CutGff;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.CobalaminConstructionInputs;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import stream.Read;
import stream.ReadInputStream;

/** Finite independent extraction oracle: shared contig names, equal DNA with
 * distinct source IDs, both strands, unknown taxonomy and native truncation.
 * Negative controls keep record counts unchanged while substituting IDs, and
 * distinguish missing, duplicate, length and absent-contig observations.
 * @author Brian Bushnell, Raiden
 */
public final class CobalaminConstructionInputsTest{

	public static void main(final String[] args) throws Exception{
		final PreParser pp=new PreParser(args, CobalaminConstructionInputsTest.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				check(parser.parse(arg, a, b), "Unknown fixture argument: "+arg);
			}
			check(parser.out1!=null && !parser.overwrite && !parser.append && parser.in1==null,
				"Supply out=fresh_fixture_directory without overwrite or input");
			check(new File(parser.out1).mkdir(), "Fixture needs one fresh output directory");
			final CobalaminConstructionInputsTest test=new CobalaminConstructionInputsTest(new File(parser.out1).getAbsolutePath());
			test.run();
			pp.outstream.println("COBALAMIN_CONSTRUCTION_INPUTS_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private CobalaminConstructionInputsTest(final String path){root=path;}

	private void run() throws Exception{
		write("a.fa", ">contig\nAACCGGTTAACCGGTT\n");
		write("b.fa", ">contig\nAACCGGTTAACCGGTT\n");
		write("c.fa", ">contig\nTTTT\n");
		write("manifest.tsv", MANIFEST+"\n"+assembly(0, "a", "123")+assembly(1, "b", ".")+assembly(2, "c", "."));
		final String raw=LOCUS+"\n"+locus(0, "a", "RF00174", "contig", 2, 5, "+", true)
			+locus(0, "a", "RF01482", "contig", 2, 5, "-", true)
			+locus(0, "a", "RF01689", "contig", 6, 9, "+", false)
			+locus(1, "b", "RF00174", "contig", 2, 5, "+", true)
			+locus(1, "b", "RF01689", "contig", 10, 13, "-", true);
		write("loci.tsv", raw);
		invoke("prepare", "prepared", null, "loci.tsv");
		CutGff.main(new String[]{"config="+path("prepared/chunks/0.config"), "types=ncRNA", "flank=0", "align=f",
			"banpartial=f", "oneperfile=f", "pickbest=f", "minlen=1", "maxlen=2147483647", "maxns=-1", "maxnfraction=-1",
			"filename=f", "gccontig=f", "renamebytaxid=f", "adjustendpoints=f", "overwrite=f", "t=4"});
		invoke("validate", "valid", "prepared/chunks/0.fa", "loci.tsv");
		final ArrayList<Read> main=reads("valid/RF00174.fa"), first=reads("valid/RF01482.fa"), second=reads("valid/RF01689.fa");
		check(main.size()==2 && first.size()==1 && second.size()==1, "Every original vote survives; native truncation alone is excluded");
		boolean known=false, unknown=false;
		for(final Read read:main){
			check(dna(read).equals("ACCG"), "Plus-strand inclusive span2..5 must be ACCG");
			known|=read.id.equals("RF00174_input_0_locus_0 tid_123");
			unknown|=read.id.equals("RF00174_input_1_locus_3");
		}
		check(known && unknown, "Equal DNA retains both distinct IDs and unknown taxonomy");
		check(first.get(0).id.equals("RF01482_input_0_locus_1 tid_123") && dna(first.get(0)).equals("CGGT"),
			"Minus-strand2..5 must reverse-complement ACCG to CGGT exactly once");
		check(second.get(0).id.equals("RF01689_input_1_locus_4") && dna(second.get(0)).equals("CGGT"),
			"The second minor-family original must retain its own ID and sense orientation");
		final String a=">RF00174_input_0_locus_0 tid_123\nACCG\n";
		final String b=">RF01482_input_0_locus_1 tid_123\nCGGT\n";
		final String c=">RF00174_input_1_locus_3\nACCG\n";
		final String d=">RF01689_input_1_locus_4\nCGGT\n";
		write("substitution.fa", a+b+d+">foreign_id\nACCG\n");
		reject("substitution", "substitution.fa", "loci.tsv", "MISSING\tRF00174_input_1_locus_3", "UNEXPECTED\tforeign_id");
		write("duplicate.fa", a+b+a+d);
		reject("duplicate", "duplicate.fa", "loci.tsv", "DUPLICATE\tRF00174_input_0_locus_0", "MISSING\tRF00174_input_1_locus_3");
		write("length.fa", a.replace("ACCG", "ACC")+b+c+d);
		reject("length", "length.fa", "loci.tsv", "LENGTH_MISMATCH\tRF00174_input_0_locus_0", "\t4\t3\t1\t16");
		write("absent_loci.tsv", raw.replace("RF00174\tcontig\t2\t5\t+\t1\t1", "RF00174\tabsent\t2\t5\t+\t1\t1"));
		write("absent.fa", b+d);
		reject("absent", "absent.fa", "absent_loci.tsv", "MISSING\tRF00174_input_0_locus_0", "\tABSENT");
		boundary();
		selectedBoundary();
	}

	/** Exact member selection preserves default IDs, sense bytes and clipped flanks.
	 * Complete-ledger validation must still reject bad unselected rows/counts.
	 */
	private void selectedBoundary() throws Exception{
		final String[] members={"RF00174", "RF01482", "RF01689"};
		for(final String member:members){
			final String prepared="selected_"+member, valid=prepared+"_valid";
			invoke("prepare", prepared, null, "loci.tsv", "flank=3", "member="+member);
			CutGff.main(new String[]{"config="+path(prepared+"/chunks/0.config"), "types=ncRNA", "flank=3", "align=f",
				"banpartial=f", "oneperfile=f", "pickbest=f", "minlen=1", "maxlen=2147483647", "maxns=-1", "maxnfraction=-1",
				"filename=f", "gccontig=caller", "renamebytaxid=f", "adjustendpoints=f", "overwrite=f", "t=4"});
			invoke("validate", valid, prepared+"/chunks/0.fa", "loci.tsv", "flank=3", "member="+member,
				"cores="+path("valid/"+member+".fa"));
			final ArrayList<Read> expected=reads("boundary_valid/"+member+".fa"), selected=reads(valid+"/"+member+".fa");
			check(selected.size()==expected.size(), "Selecting one member must retain every one of its original IDs");
			for(final Read read:selected){
				boolean found=false;
				for(final Read original:expected){
					if(!original.id.equals(read.id)){continue;}
					check(!found && Arrays.equals(original.bases, read.bases),
						"Selected contexts must preserve the default source/header/flank/GC and complete sense-oriented sequence");
					found=true;
				}
				check(found, "Member selection must not renumber IDs after intervening minor or truncated loci");
			}
			for(final String other:members){
				if(other.equals(member)){continue;}
				check(!new File(path(valid+"/"+other+".fa")).exists()
					&& !new File(path(valid+"/"+other+".genomes.tsv")).exists()
					&& !new File(path(valid+"/"+other+".loci.tsv")).exists(),
					"A selected-member validation must not write unselected FASTA or learning metadata");
			}
			final String geometry=new String(Files.readAllBytes(new File(path(valid+"/context_geometry.tsv")).toPath()), StandardCharsets.US_ASCII);
			final String row=member.equals("RF00174") ? "RF00174\t3\t2\t2\t0\t0\t0\n"
				: member.equals("RF01482") ? "RF01482\t3\t1\t0\t1\t0\t0\n" : "RF01689\t3\t1\t0\t0\t0\t1\n";
			check(geometry.endsWith(row), "Observed sense-left/right physical clipping has independent literal counts for each strand/member");
		}
		final String raw=new String(Files.readAllBytes(new File(path("loci.tsv")).toPath()), StandardCharsets.US_ASCII);
		write("bad_unselected.tsv", raw.replace("RF01482\tcontig\t2\t5", "RF01482\tcontig\t5\t2"));
		boolean rejected=false;
		try{invoke("prepare", "bad_unselected", null, "bad_unselected.tsv", "member=RF00174");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path("bad_unselected")).exists(), "Selecting main cannot hide malformed unselected original loci");
		rejected=false;
		try{invoke("prepare", "bad_minor_count", null, "loci.tsv", "member=RF00174", "expected=2,0,1");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path("bad_minor_count")).exists(), "All three authoritative counts must be checked before member selection");
		rejected=false;
		try{invoke("prepare", "unknown_member", null, "loci.tsv", "member=RF99999");}
		catch(IllegalArgumentException expected){rejected=true;}
		check(rejected && !new File(path("unknown_member")).exists(), "Unknown selected members cannot create output");
	}

	/** Literal strand/edge oracle over the original 16-base source contigs, plus
	 * rejected core substitution, source substitution, duplicate and missing records.
	 */
	private void boundary() throws Exception{
		invoke("prepare", "boundary_prepared", null, "loci.tsv", "flank=3");
		CutGff.main(new String[]{"config="+path("boundary_prepared/chunks/0.config"), "types=ncRNA", "flank=3", "align=f",
			"banpartial=f", "oneperfile=f", "pickbest=f", "minlen=1", "maxlen=2147483647", "maxns=-1", "maxnfraction=-1",
			"filename=f", "gccontig=caller", "renamebytaxid=f", "adjustendpoints=f", "overwrite=f", "t=4"});
		final String cores="cores="+path("valid/RF00174.fa")+","+path("valid/RF01482.fa")+","+path("valid/RF01689.fa");
		invoke("validate", "boundary_valid", "boundary_prepared/chunks/0.fa", "loci.tsv", "flank=3", cores);
		final Read plus=reads("boundary_valid/RF00174.fa").get(0);
		final Read minus=reads("boundary_valid/RF01482.fa").get(0);
		final Read last=reads("boundary_valid/RF01689.fa").get(0);
		check(dna(plus).equals("AACCGGTT") && plus.id.contains(" lflank=1 rflank=3 contig_gc=0.5"),
			"Plus2..5 clips its genomic-left flank to one base; whole-contig GC is one half");
		check(dna(minus).equals("AACCGGTT") && minus.id.contains("source=genome_0 lflank=3 rflank=1 contig_gc=0.5"),
			"Minus2..5 swaps the clipped flank into the sense-oriented right side");
		check(dna(last).equals("AACCGGTTAA") && last.id.contains("source=genome_1 lflank=3 rflank=3 contig_gc=0.5"),
			"Minus10..13 retains source identity, both full flanks and exactly one reverse complement");
		final String context=new String(Files.readAllBytes(new File(path("boundary_prepared/chunks/0.fa")).toPath()), StandardCharsets.US_ASCII);
		final int next=context.indexOf('>', 1);
		check(next>0, "The independent negative controls require multiple extracted contexts");
		write("boundary_bad_core.fa", context.replace("AACCGGTT", "AATCGGTT"));
		write("boundary_bad_source.fa", context.replace("source=genome_0", "source=genome_1"));
		write("boundary_duplicate.fa", context+context.substring(0, next));
		write("boundary_missing.fa", context.substring(next));
		for(final String name:new String[]{"boundary_bad_core", "boundary_bad_source", "boundary_duplicate", "boundary_missing"}){
			boolean rejected=false;
			try{invoke("validate", name, name+".fa", "loci.tsv", "flank=3", cores);}
			catch(IllegalArgumentException expected){rejected=true;}
			check(rejected && !new File(path(name+"/complete.tsv")).exists(), "Invalid native context must not receive an acceptance marker: "+name);
		}
	}

	/** Expected failure is independent of an output marker; verify the named rows. */
	private void reject(final String name, final String fasta, final String raw, final String first, final String second) throws Exception{
		boolean rejected=false;
		try{invoke("validate", name, fasta, raw);}
		catch(IllegalArgumentException expected){rejected=expected.getMessage().contains("source-ID set/length validation failed");}
		check(rejected && !new File(path(name+"/complete.tsv")).exists(), "Malformed extraction must fail without an acceptance marker: "+name);
		final String differences=new String(Files.readAllBytes(new File(path(name+"/id_differences.tsv")).toPath()), StandardCharsets.UTF_8);
		check(differences.contains(first) && differences.contains(second), "Missing independent expected diagnostic rows: "+name);
	}
	private void invoke(final String mode, final String out, final String fasta, final String raw, final String... extra) throws Exception{
		final ArrayList<String> args=new ArrayList<String>();
		args.add("mode="+mode);
		args.add("in="+path("manifest.tsv"));
		args.add("loci="+path(raw));
		args.add("out="+path(out));
		args.add("expectedgenomes=3");
		args.add("expected=2,1,1");
		args.add("chunkgenomes=3");
		args.add("seed=17");
		args.add("chunk=0");
		args.add("t=1");
		if(fasta!=null){args.add("sequences="+path(fasta));}
		for(final String option:extra){args.add(option);}
		CobalaminConstructionInputs.main(args.toArray(new String[0]));
	}
	private String assembly(final int id, final String name, final String taxid){
		return "genome_"+id+"\t"+id+"\ttest\t"+path(name+".fa")+"\t"+path(name+".gff")+"\t"+taxid+"\t.\t"+path(name+".fa")+"\t"+path(name+".gff")+"\tFASTA\tncRNA\ttest\n";
	}
	private String locus(final int id, final String name, final String member, final String contig, final int start, final int stop, final String strand, final boolean full){
		return "genome_"+id+"\t"+id+"\ttest\t"+path(name+".fa")+"\t"+member+"\t"+contig+"\t"+start+"\t"+stop+"\t"+strand+"\t1\t"+(full ? 1 : 0)+"\t"+full+"\t0\n";
	}
	private ArrayList<Read> reads(final String file){return ReadInputStream.toReads(path(file), FileFormat.FASTA, -1);}
	private static String dna(final Read read){return new String(read.bases, StandardCharsets.US_ASCII);}
	private void write(final String file, final String contents){try(Output out=new Output(path(file))){out.writer.print(contents);}}
	private String path(final String file){return new File(root, file).getPath();}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private final String root;
	private static final String MANIFEST="assembly_id\tinput_id\tcohort\tassembly\tnew_gff\ttaxids\tassembly_accession\tsource_fasta\tsource_new_gff\tassembly_format\tgff_scope\tselection_role";
	private static final String LOCUS="assembly_id\tinput_id\tcohort\tsource_fasta\tmember\tcontig\tstart\tstop\tstrand\tnative_events\tuntruncated_events\tconstruction_eligible\texact_scored_links";
}

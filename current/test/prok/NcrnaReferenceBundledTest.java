package test.prok;

import java.io.File;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.BuildNcrnaReference;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaReferenceGff;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import prok.PrepareCmSearchBatches;
import shared.Shared;

/** Native bundled-GA schema fixtures. No CM computation or real digest claim.
 * Both scan groups deliberately use batch zero and the same original contig name.
 * Original/prepared FASTAs are absent: import and tiers must use captured maps.
 * @author Ganyu
 */
public final class NcrnaReferenceBundledTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaReferenceBundledTest.class, false);
		try{
			boolean assertions=false;
			assert(assertions=true);
			check(assertions, "Run with -ea so the worker invariants are exercised");
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] fields=Parse.splitOnFirst(arg, '=');
				check(parser.parse(arg, fields[0].toLowerCase(Locale.ROOT), fields.length>1 ? fields[1] : null),
					"Unknown test parameter: "+arg);
			}
			check(parser.out1!=null && parser.in1==null && !parser.overwrite && !parser.append,
				"Supply out=fresh_fixture_directory without input or overwrite");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture root must be fresh");
			run(root.getPath());
			System.err.println("NCRNA_REFERENCE_BUNDLED_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final String root){
		final Fixture good=new Fixture(root, "good");
		good.importEvidence();
		metric(good.dir+"/evidence.summary.tsv", 0, 2);
		metric(good.dir+"/evidence.summary.tsv", 4, 2);
		metric(good.dir+"/evidence.summary.tsv", 13, 4);
		metric(good.dir+"/evidence.summary.tsv", 14, 2);
		final LineParser1 p=new LineParser1('\t');
		final boolean[] seen=new boolean[2];
		try(Input in=new Input(good.dir+"/evidence.hits.tsv.gz")){
			in.reader.nextLine();
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				final int id=p.parseInt(0)-1;
				check(id>=0 && id<2 && !seen[id], "Each original input gets its own hit despite equal contig and batch names");
				seen[id]=true;
				check(p.parseLong(2)==2 && p.termEquals("shared", 3) && p.parseInt(4)==20 && p.parseInt(5)==119
					&& p.termEquals(id==0 ? "+" : "-", 6) && p.termEquals("true", 22)
					&& p.termEquals(good.dir+'/'+(id==0 ? "canary" : "fleet")+"/native/hits.tbl", 1),
					"Native table/line, original target, strand and literal coordinates must survive import");
			}
		}
		check(seen[0] && seen[1], "Both scan groups must retain their independent original input");
		try(Input in=new Input(good.dir+"/evidence.batch_sources.tsv")){
			in.reader.nextLine();
			p.set(in.reader.nextLine());
			check(p.termEquals("canary:0", 0) && p.termEquals("0", 12), "Qualified identity retains local native batch zero");
			p.set(in.reader.nextLine());
			check(p.termEquals("fleet:0", 0) && p.termEquals("0", 12), "Second scan's zero must not collide");
			check(in.reader.nextLine()==null, "Exactly two native batch sources");
		}
		good.invoke("attributes", good.dir+"/attributes");
		good.invoke("tiers", good.dir+"/tiers100", "attributes="+good.dir+"/attributes",
			"priors="+good.dir+"/priors.tsv", "marginalpercent=100");
		final NcrnaReferenceGff codec=new NcrnaReferenceGff();
		for(int id=1; id<=2; id++){
			final ArrayList<NcrnaReferenceGff.Row> rows=codec.read(good.dir+"/tiers100/input_"+id+".gff.gz", Integer.toString(id));
			check(rows.size()==1 && !rows.get(0).mask && rows.get(0).locus.start==20 && rows.get(0).locus.stop==119,
				"An intact included hit exactly at GA remains a scored native-supported locus");
		}
		contains(good.dir+"/tiers100/parameters.tsv", "captured_reporting_floor_percent\t100");
		contains(good.dir+"/tiers100/parameters.tsv", "ga_only_marginal_meaning\tNATIVE_GA_INCLUDED_BUT_FAILED_FULL_SPAN_TRUNCATION_COMPETITOR_OR_SOURCE_CHECKS_AND_CLEARED_PRINTED_FLOOR_UNCERTAINTY");
		contains(good.dir+"/tiers100/parameters.tsv", "construction_eligibility\tUSE_NAMED_INCLUSION_SPAN_TRUNCATION_COMPETITOR_SOURCE_AND_PRIOR_CONDITIONS_NOT_TIER_WORD_ALONE");
		reject(()->good.invoke("tiers", good.dir+"/tiers90", "attributes="+good.dir+"/attributes",
			"priors="+good.dir+"/priors.tsv", "marginalpercent=90"), "sub-GA reporting coverage is absent");
		check(!new File(good.dir+"/tiers90").exists(), "Unsupported coverage must fail before creating tiers output");
		final Fixture empty=new Fixture(root, "empty");
		empty.importEvidence();
		metric(empty.dir+"/evidence.summary.tsv", 4, 0);
		metric(empty.dir+"/evidence.summary.tsv", 14, 4);
		empty.invoke("attributes", empty.dir+"/attributes");
		empty.invoke("tiers", empty.dir+"/tiers100", "attributes="+empty.dir+"/attributes",
			"priors="+empty.dir+"/priors.tsv", "marginalpercent=100");
		check(codec.read(empty.dir+"/tiers100/input_1.gff.gz", "1").isEmpty(), "An empty complete batch must not invent loci");
		final String[] faults={"missing_model", "report_count", "unknown_target", "wrong_space", "wrong_receipt", "wrong_ga", "excluded_hit", "repeated_group", "missing_query", "table_count",
			"truncated_row", "missing_ok", "out_of_range", "duplicate_target"};
		final String[] reasons={"Missing model completion", "completion hit counts differ", "Unknown synthetic native target",
			"search space differs", "foreign native input receipt", "CM metadata differs", "natively excluded hit", "Repeated scan-qualified batch",
			"Missing native query", "Bundled native hit count differs", "CM row has fewer than17 fields",
			"Native report lacks", "Native coordinates leave the captured contig", "Repeated synthetic target"};
		for(int i=0; i<faults.length; i++){
			final Fixture bad=new Fixture(root, faults[i]);
			reject(()->bad.importEvidence(), reasons[i]);
			check(!new File(bad.dir+"/evidence.summary.tsv").exists(), "Invalid native input must not publish a completed summary");
		}
	}

	private static final class Fixture{
		Fixture(final String root, final String fault){
			final boolean empty=fault.equals("empty");
			dir=root+'/'+fault;
			check(new File(dir).mkdir(), "Fresh case directory");
			write(dir+"/scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+"\nRF00050\tFMN\tPROK\tfixture\tfmn\tTARGET_CLASS\nRF00001\t5S\tPROK\tfixture\t5S\tRELATED_CLASS\n");
			write(dir+"/assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+"\nRF00050\t15.1\t100\t20\t"+SHA+'\t'+dir+"/RF00050.cm\nRF00001\t15.1\t100\t30\t"+SHA+'\t'+dir+"/RF00001.cm\n");
			write(dir+"/queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+"\nRF00050\tFMN\t100\nRF00001\t5S\t100\n");
			write(dir+"/bundle.cm", "INFERNAL1/a [fixture]\nNAME FMN\nACC RF00050\nCLEN 100\nGA "+(fault.equals("wrong_ga") ? 19 : 20)
				+"\n//\nHMMER3/f [fixture]\nNAME FMN\nACC RF00050\n//\nINFERNAL1/a [fixture]\nNAME 5S\nACC RF00001\nCLEN 100\nGA 30\n//\n");
			write(dir+"/priors.tsv", PRIOR_HEADER+'\n');
			final StringBuilder census=new StringBuilder(NcrnaEvidenceCatalog.CENSUS_HEADER+'\n');
			final StringBuilder members=new StringBuilder(NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n');
			final StringBuilder collector=new StringBuilder(BuildNcrnaReference.BUNDLED_COLLECTOR_HEADER+'\n');
			for(int id=1; id<=2; id++){
				final String group=id==1 ? "canary" : "fleet";
				final String prepared=dir+'/'+group+"/prepared";
				final String nativeDir=dir+'/'+group+"/native";
				check(new File(prepared).mkdirs() && new File(nativeDir).mkdir(), "Fresh native/preparation directories");
				final String fasta=dir+"/absent_"+id+".fa";
				final String censusRow=id+"\t"+fasta+"\t10\t300\t1\n";
				census.append(censusRow);
				members.append("fixture\t"+id+'\t'+id+'\t'+fasta+"\tg"+id+"\ta"+id+'\t'+id+"\tPROK\tNA\tNA\n");
				collector.append(fault.equals("repeated_group") ? "canary" : group).append("\t0\t"+prepared+'\t'+nativeDir+"\t1\t300\t1\tRAW_NATIVE_GA_COMPLETE\n");
				write(prepared+"/members.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+censusRow);
				write(prepared+"/summary.tsv", "batch_id\ttask_id\tgenomes\tnucleotide_bases\tfasta_records\tstatus\n0\t0\t1\t300\t1\tPREPARED_NOT_SEARCHED\n");
				final String target="0\ti"+id+"c1\t"+id+"\t1\tshared\t300\n";
				write(prepared+"/targets.tsv", PrepareCmSearchBatches.TARGET_HEADER+'\n'+target+(fault.equals("duplicate_target") ? target : ""));
				write(nativeDir+"/input.tsv", "role\tpath\tsha80\nbatch.fa\t"+prepared+"/batch.fa\t"+SHA+"\ntargets.tsv\t"+prepared
					+(fault.equals("wrong_receipt") ? "/foreign.tsv" : "/targets.tsv")+'\t'+SHA+"\nmembers.tsv\t"+prepared+"/members.tsv\t"+SHA+'\n');
				write(nativeDir+"/exit.txt", "0\n");
				write(nativeDir+"/completion.tsv", COMPLETION+"\n0\tRF00050\t1\t300\t600\t"+(empty ? 0 : 1)+"\tCOMPLETE\n"
					+(fault.equals("missing_model") ? "" : "0\tRF00001\t1\t300\t600\t0\tCOMPLETE\n"));
				final String command="cmsearch --cpu 1 --cut_ga --nohmmonly --noali --notextw --tformat fasta -Z "
					+(fault.equals("wrong_space") ? "0.002" : "0.001")+" --tblout "+nativeDir+"/hits.tbl -o "+nativeDir+"/search.out "+dir+"/bundle.cm "+prepared+"/batch.fa";
				final String hit=(fault.equals("unknown_target") ? "i99c1" : "i"+id+"c1")+" - FMN RF00050 cm 1 100 "
					+(fault.equals("out_of_range") ? "220 319 +" : id==1 ? "20 119 +" : "119 20 -")+" no 1 0.5 0.0 20.0 1e-5 "+(fault.equals("excluded_hit") ? "?" : "!")+" fixture\n";
				write(nativeDir+"/hits.tbl", "# synthetic report\n"+(empty ? "" : hit)+(fault.equals("table_count") ? hit : "")
					+"# Program: cmsearch\n# Version: 1.1.5 (fixture)\n# Option settings: "+command+'\n'
					+(fault.equals("truncated_row") ? "i"+id+"c1 - FMN RF00050 cm 1" : fault.equals("missing_ok") ? "" : "# [ok]\n"));
				write(nativeDir+"/search.out", "Query: FMN [CLEN=100]\nAccession: RF00050\nTarget sequences: 1 (600 residues searched)\nTotal CM hits reported: "
					+(fault.equals("report_count") ? 2 : empty ? 0 : 1)
					+(fault.equals("missing_query") ? "" : "\nQuery: 5S [CLEN=100]\nAccession: RF00001\nTarget sequences: 1 (600 residues searched)\nTotal CM hits reported: 0")+"\n[ok]\n");
			}
			write(dir+"/census.tsv", census.toString());
			write(dir+"/memberships.tsv", members.toString());
			write(dir+"/collector.tsv", collector.toString());
		}
		void importEvidence(){invoke("batchevidence", dir+"/evidence", "cmbundle="+dir+"/bundle.cm", "searchspacemb=0.001");}
		void invoke(final String mode, final String output, final String... extra){
			final ArrayList<String> args=new ArrayList<String>(Arrays.asList("mode="+mode,
				"in="+dir+(mode.equals("batchevidence") ? "/collector.tsv" : "/evidence"), "out="+output,
				"census="+dir+"/census.tsv", "memberships="+dir+"/memberships.tsv", "scope="+dir+"/scope.tsv",
				"assets="+dir+"/assets.tsv", "queries="+dir+"/queries.tsv", "t=1", "pigz=f", "bgzip=f"));
			args.addAll(Arrays.asList(extra));
			BuildNcrnaReference.main(args.toArray(new String[args.size()]));
		}
		final String dir;
	}
	private static void write(final String path, final String value){
		try(Output out=new Output(path)){out.writer.print(value.getBytes(StandardCharsets.UTF_8));}
	}
	private static void metric(final String path, final int column, final long expected){
		try(Input in=new Input(path)){
			in.reader.nextLine();
			final LineParser1 p=new LineParser1('\t');
			p.set(in.reader.nextLine());
			check(p.parseLong(column)==expected && in.reader.nextLine()==null, "Unexpected summary value at column "+column);
		}
	}
	private static void contains(final String path, final String expected){
		try(Input in=new Input(path)){
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				if(expected.equals(new String(row, StandardCharsets.UTF_8))){return;}
			}
		}
		throw new AssertionError("Missing literal metadata: "+expected);
	}
	private static void reject(final Runnable action, final String reason){
		try{action.run();}
		catch(RuntimeException failure){
			for(Throwable cause=failure; cause!=null; cause=cause.getCause()){
				if(cause.getMessage()!=null && cause.getMessage().contains(reason)){return;}
			}
			throw new AssertionError("Wrong rejection: "+reason, failure);
		}
		throw new AssertionError("Expected rejection: "+reason);
	}
	private static void check(final boolean valid, final String why){if(!valid){throw new AssertionError(why);}}
	private static final String SHA="00000000000000000001";
	private static final String COMPLETION="batch_id\tmodel\tFASTA_records\tinput_bases\tsearched_residues\treported_hits\tstatus";
	private static final String PRIOR_HEADER="source_id\tdecision_file\tdecision_line\toriginal_locus\toriginal_genome\toriginal_fasta\tinput_id\tresolved_fasta\tliteral_target\tcontig_length\tsource_family\tassignment\tsource_start1\tsource_stop1\tcore_start1\tcore_stop1\tstrand\tinterval_semantics\tsource_authority\trow_authority\toriginal_row";
}

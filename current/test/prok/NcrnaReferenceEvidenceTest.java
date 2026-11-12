package test.prok;

import java.io.File;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaEvidenceCatalog.Observation;
import prok.NcrnaReferenceEvidence;
import prok.NcrnaReferenceEvidence.Evidence;
import prok.NcrnaReferenceEvidence.Roster;
import prok.NcrnaTextIO.Output;
import shared.Shared;

/** Independent finite expectations for native CM import and evidence-only model roles.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaReferenceEvidenceTest{
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaReferenceEvidenceTest.class, false);
		try{
			final Parser parser=new Parser();
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
			}
			check(parser.out1!=null, "Supply out=new_fixture_directory");
			final File root=new File(parser.out1).getAbsoluteFile();
			check(root.mkdir(), "Fixture directory must be new: "+root);
			Shared.setThreads(1);
			run(root);
			pp.outstream.println("NCRNA_REFERENCE_EVIDENCE_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String fasta=write(root, "one.fa", ">seq description\nAAAAAAAAAAAAAAAAAAAA\n");
		final String other=write(root, "two.fa", ">seq\nAAAAAAAAAAAAAAAAAAAA\n");
		final String census=write(root, "census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'
			+"1\t"+fasta+'\t'+new File(fasta).length()+"\t20\t1\n"
			+"2\t"+other+'\t'+new File(other).length()+"\t20\t1\n");
		final String members=write(root, "members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'
			+"fixture\t1\t1\t"+fasta+"\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n"
			+"fixture\t2\t2\t"+other+"\tNA\tNA\tNA\tprokaryotic\tNA\tNA\n");
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, members);
		final Observation first=catalog.observe(catalog.get("1")), second=catalog.observe(catalog.get("2"));
		final String scope=write(root, "scope.tsv", NcrnaReferenceEvidence.SCOPE_HEADER+'\n'
			+"RF00504\tglycine_aptamer\tprimary_glycine_evidence\taptamer_only\tglycine\tTARGET_CLASS\n"
			+"RF03166\tglycine_GGAnGA\trelated_class_arbitration\tnever_glycine_truth\tglycine_GGAnGA\tRELATED_CLASS\n"
			+"RF00030\tRNase_MRP\trelated_class_arbitration\tnever_RNaseP_truth\tRNaseMRP\tRELATED_CLASS\n");
		final String assets=write(root, "assets.tsv", NcrnaReferenceEvidence.ASSET_HEADER+'\n'
			+"RF00504\tfixture\t10\t20\t00000000000000000001\t/fixture/glycine.cm\n"
			+"RF03166\tfixture\t6\t30\t00000000000000000002\t/fixture/related.cm\n"
			+"RF00030\tfixture\t8\t15\t00000000000000000003\t/fixture/mrp.cm\n");
		final String queries=write(root, "queries.tsv", NcrnaReferenceEvidence.QUERY_HEADER+'\n'
			+"RF00504\tglycine\t10\nRF03166\trelated\t6\n");
		final Roster roster=NcrnaReferenceEvidence.loadModels(scope, assets, queries);
		check(roster.size()==3 && roster.queried()==2, "Captured models and actually searched models are distinct denominators");
		check(roster.model("RF00030").relatedClass() && roster.model("RF03166").relatedClass()
			&& !roster.model("RF00504").relatedClass(), "Related classes must never be relabeled as target classes");
		check(roster.model("RF03166").scope.scoringScope.equals("never_glycine_truth"), "Model restrictions must survive normalization");

		final String intact=hit("seq", "glycine", "RF00504", 1, 10, 1, 10, "+", "no", 25, "1e-8", "!");
		final String related=hit("seq", "related", "RF03166", 1, 6, 12, 7, "-", "no", 70, "1e-12", "!");
		final String partial=hit("seq", "glycine", "RF00504", 2, 8, 15, 21, "+", "3'", 19, "0.1", "?");
		final String unresolved=hit("seq_tid_7", "glycine", "RF00504", 1, 10, 1, 10, "+", "no", 30, "1e-9", "!");
		final String rows=intact+related+partial+unresolved;
		final String counts=completion("1", 3, 1);
		final Evidence evidence=read(root, "all", rows, counts, first, roster);
		check(evidence.size()==4 && evidence.hitsForQuery(0)==3 && evidence.hitsForQuery(1)==1,
			"Raw import retains below-GA, related-class and unresolved observations");
		check(evidence.get(0).fullModelSpan && evidence.get(0).atLeftEdge && !evidence.get(0).atRightEdge,
			"Full model span and physical contig-edge contact are separate observations");
		check(evidence.get(1).strand=='-' && evidence.get(1).start==7 && evidence.get(1).stop==12
			&& evidence.get(1).gaPass && evidence.get(1).model.relatedClass(), "Strong reverse-strand related-class evidence stays related");
		check(evidence.get(2).range.equals("OUT_OF_RANGE") && !evidence.get(2).included && !evidence.get(2).gaPass
			&& evidence.get(2).truncation==2 && evidence.get(2).seqTo==21,
			"Invalid genomic ranges, below-GA scores and truncation are retained without clipping or admission");
		check(evidence.get(3).range.equals("UNRESOLVED_LITERAL_TARGET"), "A taxonomic-looking suffix must not silently alias a native target");
		final Evidence zero=read(root, "zero", intact, completion("1", 1, 0), first, roster);
		check(zero.size()==1 && zero.hitsForQuery(1)==0, "A completed zero-hit search must remain represented");

		reject(root, "missing_zero", intact, NcrnaReferenceEvidence.COMPLETION_HEADER+"\n1\tRF00504\t1\t20\t40\t1\tCOMPLETE\n", first, roster, "Missing populated or zero-hit");
		reject(root, "duplicate_completion", intact, completion("1", 1, 0)+"1\tRF03166\t1\t20\t40\t0\tCOMPLETE\n", first, roster, "repeated completion");
		reject(root, "wrong_input", intact, completion("2", 1, 0), first, roster, "identify this input");
		reject(root, "wrong_dimensions", intact, completion("1", 1, 0).replace("\t20\t40\t", "\t21\t42\t"), first, roster, "dimensions disagree");
		reject(root, "count_mismatch", intact, completion("1", 0, 0), first, roster, "retained hits disagree");
		reject(root, "unqueried", intact.replace("glycine RF00504", "mrp RF00030"), completion("1", 1, 0), first, roster, "unqueried");
		reject(root, "unknown", intact.replace("RF00504", "RF99999"), completion("1", 1, 0), first, roster, "unknown");
		reject(root, "query_name", intact.replace("glycine RF00504", "wrong RF00504"), completion("1", 1, 0), first, roster, "differently named");
		reject(root, "strand", hit("seq", "glycine", "RF00504", 1, 10, 1, 10, "-", "no", 25, "1e-8", "!"), completion("1", 1, 0), first, roster, "coordinate order");
		final String unfinished=write(root, "unfinished.tbl", HEADER+intact);
		final String done=write(root, "unfinished.completion.tsv", completion("1", 1, 0));
		expect(()->NcrnaReferenceEvidence.read(unfinished, done, first, roster), "completion footer");

		final String ranks=hit("seq", "glycine", "RF00504", 1, 10, 2, 11, "+", "no", 30, "1e-8", "!")
			+hit("seq", "glycine", "RF00504", 2, 9, 2, 11, "+", "no", 100, "1e-50", "!")
			+hit("seq", "glycine", "RF00504", 1, 10, 2, 11, "+", "no", 35, "1e-5", "!")
			+hit("seq", "glycine", "RF00504", 1, 10, 2, 11, "+", "no", 35, "1e-8", "!")
			+hit("seq", "glycine", "RF00504", 1, 10, 8, 17, "+", "no", 80, "1e-20", "!");
		final Evidence rank=read(root, "rank", ranks, completion("1", 5, 0), first, roster);
		check(NcrnaReferenceEvidence.compareSupport(rank.get(0), rank.get(1), 2, 11)<0, "Model coverage precedes bit score");
		check(NcrnaReferenceEvidence.compareSupport(rank.get(2), rank.get(0), 2, 11)<0, "Bit score precedes E-value");
		check(NcrnaReferenceEvidence.compareSupport(rank.get(3), rank.get(2), 2, 11)<0, "Lower E-value breaks equal coverage/score ties");
		check(NcrnaReferenceEvidence.compareSupport(rank.get(0), rank.get(4), 2, 11)<0, "Actual historical overlap precedes stronger scores elsewhere");
		final String differentLengths=hit("seq", "glycine", "RF00504", 1, 6, 2, 11, "+", "no", 50, "1e-8", "!")
			+hit("seq", "related", "RF03166", 1, 5, 2, 11, "+", "no", 40, "1e-8", "!")
			+hit("seq", "related", "RF03166", 1, 6, 2, 11, "+", "no", 40, "1e-8", "!")
			+hit("seq", "glycine", "RF00504", 1, 10, 2, 11, "+", "no", 40, "1e-8", "!");
		final Evidence normalized=read(root, "normalized", differentLengths, completion("1", 2, 2), first, roster);
		check(NcrnaReferenceEvidence.compareSupport(normalized.get(1), normalized.get(0), 2, 11)<0
			&& NcrnaReferenceEvidence.compareSupport(normalized.get(0), normalized.get(1), 2, 11)>0,
			"Five of six model positions outrank six of ten, despite the smaller raw span and lower bit score");
		check(NcrnaReferenceEvidence.compareSupport(normalized.get(3), normalized.get(2), 2, 11)<0
			&& NcrnaReferenceEvidence.compareSupport(normalized.get(2), normalized.get(3), 2, 11)>0,
			"Equal-support models use accession lexical order even when it opposes native row order");
		final Evidence foreign=read(root, "foreign", intact, completion("2", 1, 0), second, roster);
		expect(()->NcrnaReferenceEvidence.compareSupport(evidence.get(0), foreign.get(0), 1, 10), "one input");
	}

	private static Evidence read(final File root, final String name, final String rows, final String counts, final Observation observed, final Roster roster){
		final String table=write(root, name+".tbl", HEADER+rows+"# [ok]\n");
		final String completion=write(root, name+".completion.tsv", counts);
		return NcrnaReferenceEvidence.read(table, completion, observed, roster);
	}
	private static void reject(final File root, final String name, final String rows, final String counts, final Observation observed, final Roster roster, final String error){
		expect(()->read(root, name, rows, counts, observed, roster), error);
	}
	private static String completion(final String input, final int target, final int related){
		return NcrnaReferenceEvidence.COMPLETION_HEADER+'\n'+input+"\tRF00504\t1\t20\t40\t"+target+"\tCOMPLETE\n"
			+input+"\tRF03166\t1\t20\t40\t"+related+"\tCOMPLETE\n";
	}
	private static String hit(final String target, final String name, final String accession, final int mf, final int mt,
		final int from, final int to, final String strand, final String truncation, final int bits, final String evalue, final String inclusion){
		return target+" - "+name+' '+accession+" cm "+mf+' '+mt+' '+from+' '+to+' '+strand+' '+truncation+" 1 0.5 0 "+bits+' '+evalue+' '+inclusion+" synthetic_fixture\n";
	}
	private static String write(final File root, final String name, final String text){
		final String path=new File(root, name).getPath();
		try(Output output=new Output(path)){output.writer.print(text);}
		return path;
	}
	private static void expect(final Runnable action, final String fragment){
		try{action.run();}
		catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(fragment), "Unexpected failure: "+expected);
			return;
		}
		throw new AssertionError("Expected failure containing: "+fragment);
	}
	private static void check(final boolean valid, final String reason){if(!valid){throw new AssertionError(reason);}}
	private static final String HEADER="# Program: cmsearch\n# Version: 1.1.5 (synthetic fixture)\n";
}

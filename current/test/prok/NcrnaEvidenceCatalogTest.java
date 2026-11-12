package test.prok;

import java.io.File;
import java.util.Locale;

import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaEvidenceCatalog;
import prok.NcrnaEvidenceCatalog.Observation;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import stream.FastaStreamerST;
import stream.Read;
import stream.StreamerFactory;
import structures.ByteBuilder;

/** Finite contracts for path identity, source claims and native contig observations.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaEvidenceCatalogTest{

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaEvidenceCatalogTest.class, false);
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
			pp.outstream.println("NCRNA_EVIDENCE_CATALOG_TEST_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	private static void run(final File root){
		final String first=write(root, "first.fa", ">same description\nACGT\n>last\nNNN");
		final String second=write(root, "second.fa.gz", ">same\nAA\n");
		final String one=censusRow("1", first, 7, 2);
		final String two=censusRow("2", second, 2, 1);
		final String census=write(root, "census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+two+one);
		final String a=member("sourceA", 1, "1", first, "/annotations/a.gff");
		final String b=member("sourceA", 2, "2", second, "NA");
		final String c=member("sourceB", 1, "1", first, "relative/other.gff");
		final String members=write(root, "members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+a+b+c);
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, members);
		check(catalog.size()==2 && catalog.memberships()==3 && catalog.get(0).id.equals("2"),
			"Distinct paths with the same assembly claim stay separate and retain census order");
		check(catalog.get("1").memberships()==2 && catalog.get("2").memberships()==1,
			"Multiple source memberships must neither disappear nor multiply genomes");
		check(catalog.get("1").membership(1).gff.equals("relative/other.gff")
			&& catalog.get("2").membership(0).gff.equals("NA")
			&& catalog.get("1").membership(0).reportedSha80.equals("NA"),
			"Unresolved GFF text, absent annotations and unrecorded hashes must remain explicit");
		final Observation x=catalog.observe(catalog.get("1"));
		final Observation y=catalog.observe(catalog.get("2"));
		check(x.size()==2 && x.length("same")==4 && x.length("last")==3 && x.length("missing")==-1,
			"Every literal FASTA record, including the final record without a newline, needs its true length");
		check(y.size()==1 && y.length("same")==2,
			"Gzip decoding and equal contig names in separate genomes must preserve separate observations");

		reject(root, "duplicate_id", two+one+one, a+b+c);
		reject(root, "duplicate_path", one+censusRow("2", first, 7, 2), a+b+c);
		reject(root, "missing_membership", one+two, a);
		reject(root, "unknown_input", one+two, a+b+member("sourceB", 1, "3", first, "NA"));
		reject(root, "duplicate_membership", one+two, a+b+a);
		reject(root, "source_index_gap", one+two, a+member("sourceA", 3, "2", second, "NA"));

		final String repeated=write(root, "repeated.fa", ">same\nACGT\n>same\nAAA\n");
		final String duplicateCensus=write(root, "repeated_census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+censusRow("1", repeated, 7, 2));
		final String duplicateMembers=write(root, "repeated_members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+member("sourceA", 1, "1", repeated, "NA"));
		final NcrnaEvidenceCatalog duplicate=NcrnaEvidenceCatalog.load(duplicateCensus, duplicateMembers);
		expectMessage(()->duplicate.observe(duplicate.get("1")), "Repeated literal contig name",
			"Duplicate literal names within one genome must fail after draining input");
		queuedDuplicate(root);
		final boolean fixHeader=Read.FIX_HEADER;
		try{
			Read.FIX_HEADER=true;
			expectMessage(()->catalog.observe(catalog.get("1")), "Read.FIX_HEADER=false",
				"Header-rewriting mode cannot silently change native CM target identifiers");
			check(Read.FIX_HEADER, "The observer must reject the setting without changing shared configuration");
		}finally{Read.FIX_HEADER=fixHeader;}
		final String wrongCensus=write(root, "wrong_bases.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+censusRow("1", first, 8, 2));
		final String firstMember=write(root, "first_member.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+a);
		final NcrnaEvidenceCatalog wrong=NcrnaEvidenceCatalog.load(wrongCensus, firstMember);
		expectFailure(()->wrong.observe(wrong.get("1")), "Observed bases must match the independent census denominator");
	}

	/** A duplicate in the first batch leaves sixteen more batches behind a two-entry queue. */
	private static void queuedDuplicate(final File root){
		final ByteBuilder b=new ByteBuilder().append(">repeat\nA\n>repeat\nA\n");
		for(int i=0; i<32; i++){b.append(">tail").append(i).append("\nA\n");}
		final String fasta=write(root, "queued.fa", b.toString());
		final String census=write(root, "queued.census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+censusRow("1", fasta, 34, 34));
		final String members=write(root, "queued.members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+member("queue", 1, "1", fasta, "NA"));
		final NcrnaEvidenceCatalog catalog=NcrnaEvidenceCatalog.load(census, members);
		final int threads=Shared.threads(), batchSize=FastaStreamerST.TARGET_LIST_SIZE;
		final boolean lowMemory=Shared.LOW_MEMORY, secondReader=StreamerFactory.FASTA_STREAMER_2;
		try{
			Shared.setThreads(4);
			Shared.LOW_MEMORY=false;
			StreamerFactory.FASTA_STREAMER_2=false;
			FastaStreamerST.TARGET_LIST_SIZE=2;
			expectMessage(()->catalog.observe(catalog.get("1")), "Repeated literal contig name",
				"The one-producer reader must drain its queued input and return the intended duplicate-ID failure");
		}finally{
			Shared.setThreads(threads);
			Shared.LOW_MEMORY=lowMemory;
			StreamerFactory.FASTA_STREAMER_2=secondReader;
			FastaStreamerST.TARGET_LIST_SIZE=batchSize;
		}
	}

	private static void reject(final File root, final String name, final String census, final String members){
		final String counts=write(root, name+".census.tsv", NcrnaEvidenceCatalog.CENSUS_HEADER+'\n'+census);
		final String claims=write(root, name+".members.tsv", NcrnaEvidenceCatalog.MEMBERSHIP_HEADER+'\n'+members);
		expectFailure(()->NcrnaEvidenceCatalog.load(counts, claims), "Malformed catalog accepted: "+name);
	}

	private static String censusRow(final String id, final String fasta, final long bases, final long records){
		return id+'\t'+fasta+'\t'+new File(fasta).length()+'\t'+bases+'\t'+records+'\n';
	}

	private static String member(final String source, final int index, final String id, final String fasta, final String gff){
		return source+'\t'+index+'\t'+id+'\t'+fasta+"\tshared_name\tGCF_000001.1\t7\tprokaryotic\t"+gff+"\tNA\n";
	}

	private static String write(final File root, final String name, final String text){
		final String path=new File(root, name).getPath();
		try(Output output=new Output(path)){output.writer.print(text);}
		return path;
	}

	private static void expectFailure(final Runnable action, final String reason){
		try{action.run();}
		catch(IllegalArgumentException expected){return;}
		throw new AssertionError(reason);
	}

	private static void expectMessage(final Runnable action, final String fragment, final String reason){
		try{action.run();}
		catch(IllegalArgumentException expected){
			check(expected.getMessage()!=null && expected.getMessage().contains(fragment), reason+": "+expected);
			return;
		}
		throw new AssertionError(reason);
	}

	private static void check(final boolean valid, final String reason){
		if(!valid){throw new AssertionError(reason);}
	}
}

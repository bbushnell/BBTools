package test.assemble;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Locale;

import dna.AminoAcid;
import fileIO.ByteStreamWriter;
import fileIO.FileFormat;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import shared.Shared;
import shared.Tools;
import stream.FastaReadInputStream;
import stream.Read;
import stream.ReadInputStream;
import structures.ByteBuilder;
import ukmer.HollowAnchorTable;

/**
 * Read-only truth audit for hand-selected 32-base assembly anchors.
 *
 * <p>The diagnostic loads complete contigs and references, disables FASTA split
 * reads for the process without restoring that global flag, and emits every
 * literal occurrence of each requested anchor in the supplied reference without
 * consulting or mutating the Hollow joiner. Palindromic 32-base matches produce
 * one positive-strand row per physical position; reported coordinates are
 * zero-based evidence locations and do not feed back into join selection.</p>
 *
 * @author Brian Bushnell, Fischl
 */
public final class HollowTipAnchorAudit{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	/** Parses the bounded query list, then scans complete contig and reference FASTA records. */
	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HollowTipAnchorAudit.class, false);
		try{new HollowTipAnchorAudit(pp.args, pp.outstream).process();}
		finally{Shared.closeStream(pp.outstream);}
	}

	/** Retains only contig-index:start0 query pairs selected from a candidate report. */
	private HollowTipAnchorAudit(final String[] args, final PrintStream log_){
		log=log_;
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(String arg : args){
			final String[] pair=Parse.splitOnFirst(arg, '=');
			final String a=pair[0].toLowerCase(Locale.ROOT), b=pair.length>1 ? pair[1] : null;
			if(a.equals("ref")){reference=b;}
			else if(a.equals("positions")){positions=b;}
			else if(a.equals("test")){test=Parse.parseBoolean(b);}
			else if(!parser.parse(arg, a, b)){throw new IllegalArgumentException("Unknown parameter: "+arg);}
		}
		input=parser.in1;
		output=parser.out1;
		if(parser.overwrite || parser.append){throw new IllegalArgumentException("Audit reports must be fresh.");}
		if(!test && (input==null || reference==null || output==null || positions==null)){
			throw new IllegalArgumentException("Require in=contigs.fa ref=reference.fa positions=index:start0,... out=hits.tsv.");
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------           Methods            ----------------*/
	/*--------------------------------------------------------------*/

	/** Builds selected-arm queries, then scans unsplit reference records on both strands. */
	private void process(){
		if(test){selfTest(); return;}
		Tools.testInputFiles(false, true, input, reference);
		if(!Tools.testOutputFiles(false, false, false, output)){throw new IllegalArgumentException("Output must be new: "+output);}
		FastaReadInputStream.SPLIT_READS=false;
		final ArrayList<Read> contigs=ReadInputStream.toReads(input, FileFormat.FASTA, -1);
		final String[] terms=positions.split(",");
		if(terms.length>32){throw new IllegalArgumentException("This bounded audit accepts at most32 selected anchors.");}
		final Query[] queries=new Query[terms.length];
		for(int i=0; i<terms.length; i++){
			final String[] pair=terms[i].split(":");
			if(pair.length!=2){throw new IllegalArgumentException("Bad anchor position: "+terms[i]);}
			final int id=Integer.parseInt(pair[0]), start=Integer.parseInt(pair[1]);
			if(id<0 || id>=contigs.size() || start<0 || start>contigs.get(id).length()-32){
				throw new IllegalArgumentException("Anchor position outside assembly: "+terms[i]);
			}
			queries[i]=new Query(id, start, encode(contigs.get(id).bases, start));
		}
		final ByteStreamWriter writer=new ByteStreamWriter(output, false, false, true);
		writer.start();
		final ByteBuilder row=new ByteBuilder();
		try{
			writer.println("#query_index\tcontig_index\tanchor_start0\treference_index\treference_start0\tstrand");
			final ArrayList<Read> references=ReadInputStream.toReads(reference, FileFormat.FASTA, -1);
			if(references.isEmpty()){throw new IllegalArgumentException("Reference contains no sequences.");}
			for(int r=0; r<references.size(); r++){
				final byte[] bases=references.get(r).bases;
				long word=0;
				int valid=0;
				for(int p=0; p<bases.length; p++){
					final int base=AminoAcid.baseToNumber[bases[p]];
					if(base<0){word=0; valid=0; continue;}
					word=(word<<2)|base;
					if(++valid<32){continue;}
					for(int q=0; q<queries.length; q++){
						final Query query=queries[q];
						if(word!=query.forward && word!=query.reverse){continue;}
						query.hits++;
						row.clear().append(q).tab().append(query.contig).tab().append(query.start).tab();
						row.append(r).tab().append(p-31).tab().append(word==query.forward ? '+' : '-');
						writer.println(row);
					}
				}
			}
			for(int q=0; q<queries.length; q++){
				assert(queries[q].hits>=0) : "Physical reference occurrence totals cannot be negative.";
				log.println("Anchor "+q+" ("+queries[q].contig+":"+queries[q].start+"): "+queries[q].hits+" reference occurrences");
			}
		}finally{if(writer.poisonAndWait()){throw new RuntimeException("Audit output failed: "+output);}}
	}

	/** Encodes one complete defined 32-base word; every bit pattern is a legal identity. */
	private static long encode(final byte[] bases, final int start){
		assert(start>=0 && start<=bases.length-32) : "Selected positions must contain a full32bp arm.";
		long word=0;
		for(int i=0; i<32; i++){
			final int base=AminoAcid.baseToNumber[bases[start+i]];
			if(base<0){throw new IllegalArgumentException("Selected anchor contains an undefined base.");}
			word=(word<<2)|base;
		}
		return word;
	}

	/** Checks the diagnostic word and strand conventions with known 32-base identities. */
	private void selfTest(){
		final byte[] bases=new byte[32];
		java.util.Arrays.fill(bases, (byte)'A');
		if(encode(bases, 0)!=0 || HollowAnchorTable.reverseComplement(0)!=-1L){throw new AssertionError("All-A/all-T strand convention changed.");}
		java.util.Arrays.fill(bases, (byte)'T');
		if(encode(bases, 0)!=-1L){throw new AssertionError("All-one32bp key must remain valid.");}
		bases[0]='N';
		boolean rejected=false;
		try{encode(bases, 0);}catch(IllegalArgumentException e){rejected=true;}
		if(!rejected){throw new AssertionError("Undefined query must fail rather than invent an identity.");}
		log.println("HOLLOW_TIP_ANCHOR_AUDIT_TEST_PASS");
	}

	/*--------------------------------------------------------------*/
	/*----------------         Inner classes        ----------------*/
	/*--------------------------------------------------------------*/

	/** One exact anchor query and its physical occurrence total. */
	private static final class Query{
		/** Stores one contig/start query and precomputes both strand encodings. */
		Query(final int contig_, final int start_, final long word){
			contig=contig_;
			start=start_;
			forward=word;
			reverse=HollowAnchorTable.reverseComplement(word);
		}
		/** Selected assembly contig index and zero-based 32-base start. */
		final int contig, start;
		/** Forward and reverse-complement 32-base encodings of the same query arm. */
		final long forward, reverse;
		/** Number of literal reference positions observed on either strand. */
		long hits;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	/** Input contigs, reference FASTA, fresh output TSV, and selected position list. */
	private String input, reference, output, positions;
	/** True for the built-in convention self-test that skips file I/O. */
	private boolean test;
	/** Diagnostic stream supplied by PreParser; main delegates cleanup to Shared.closeStream. */
	private final PrintStream log;
}

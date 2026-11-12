package prok;

import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.Locale;

import map.ObjectIntMap;
import map.LongObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Extracts original stage/GFF evidence around every saved changed reference locus.
 * Uses the native input-report order and explicit policy roster to find sidecars;
 * bare contig names never join different inputs. Stage spans become one-based
 * genomic closed intervals. Both strands are retained, including neighboring CDS.
 * This is an evidence view, not a new grader or a causal stage classifier. ACCEPT
 * remains pre-snapshot/pre-path; unobserved reused work remains unobserved.
 * One thread owns native byte readers, parsers, retained loci and output buffer.
 * Adapted from template.A_SampleByteFile and the shared native text owners.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaChangeEvidence{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final NcrnaChangeEvidence tool=new NcrnaChangeEvidence(args);
		try{tool.process();}
		finally{Shared.closeStream(tool.outstream);}
	}

	private NcrnaChangeEvidence(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaChangeEvidence.class, false);
		outstream=pp.outstream;
		Shared.setThreads(1);
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:pp.args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
			if(a.equals("inputs")){inputs=b;}
			else if(a.equals("policies")){policies=b;}
			else if(a.equals("calls")){calls=b;}
			else if(a.equals("family")){family=b;}
			else if(a.equals("padding")){padding=Parse.parseIntKMG(b);}
			else if(a.equals("models")){includeModels=Parse.parseBoolean(b);}
			else{require(parser.parse(arg, a, b), "Unknown flag: "+arg);}
		}
		changes=parser.in1;
		output=parser.out1;
		require(changes!=null && inputs!=null && policies!=null && calls!=null
			&& family!=null && !family.isEmpty() && output!=null && padding>=0,
			"Supply in=changes.tsv inputs=sweep.inputs.tsv policies=roster calls=sweep_prefix family=prefix out=fresh.tsv padding=nonnegative");
		require(!parser.overwrite && !parser.append, "Evidence output must be fresh");
		require(Tools.testInputFiles(false, true, changes, inputs, policies), "Cannot read exact saved declarations");
		require(Tools.testOutputFiles(false, false, false, output), "Cannot create fresh evidence output");
		require(Tools.testForDuplicateFiles(true, changes, inputs, policies, output), "Evidence declarations must be distinct files");
	}

	/*--------------------------------------------------------------*/
	/*----------------        Saved Evidence        ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		final NcrnaJointPolicy[] roster=NcrnaJointPolicy.load(policies);
		final ObjectIntMap<String> policyIndex=new ObjectIntMap<String>(String.class);
		for(int i=0; i<roster.length; i++){policyIndex.put(roster[i].name, i);}
		final ObjectIntMap<String> inputIndex=bindInputs(roster);
		final ArrayList<Change> loci=readChanges(inputIndex, policyIndex);
		final ArrayList<Batch> batches=batches(loci, calls);
		preflight(batches);
		long contexts=0, emptyContexts=0;
		try(Output out=new Output(output)){
			out.writer.println("# SAVED_ROW retains an original source row; NO_SAVED_EVIDENCE records an inspected empty context with both searched paths and source_line=0. Coordinates are one-based genomic closed; padding is inspection context, not a matching rule. Missing observations are not rejection proof. MODEL rows require models=t.");
			out.writer.println("change_row\tinput_id\tpolicy\tchange\treference_id\tseqid\tref_start\tref_stop\tref_strand\tsource_policy\tsearched_stage\tsearched_gff\tsource_path\tsource_line\tgenomic_start\tgenomic_stop\torf_genomic_start\torf_genomic_stop\tevidence_kind\tencoded_original_row");
			for(final Batch batch:batches){
				final long[] rows=new long[batch.loci.size()];
				stages(out, batch, roster[batch.policy].name, batch.stagePath, rows);
				gff(out, batch, roster[batch.policy].name, batch.gffPath, rows);
				for(int i=0; i<rows.length; i++){
					final Change locus=batch.loci.get(i);
					contexts++;
					locus.evidenceRows=Math.addExact(locus.evidenceRows, rows[i]);
					if(rows[i]==0){
						emptyContexts++;
						emit(out, batch, locus, roster[batch.policy].name, null, 0, -1, -1, -1, -1, null);
					}
				}
			}
		}
		long emptyLoci=0;
		for(final Change locus:loci){if(locus.evidenceRows==0){emptyLoci++;}}
		outstream.println("NCRNA_CHANGE_EVIDENCE_COMPLETE changes="+loci.size()+" contexts="+contexts
			+" no_saved_evidence_contexts="+emptyContexts+" no_saved_evidence_loci="+emptyLoci+"; manual classification remains required");
	}

	/** Read each original input/policy sidecar once, even when many loci changed. */
	private static ArrayList<Batch> batches(final ArrayList<Change> loci, final String prefix){
		assert(loci!=null && prefix!=null) : "Only parsed changed loci and the explicit prefix may select native evidence files";
		final LongObjectMap<Batch> byKey=new LongObjectMap<Batch>(Batch.class);
		final ArrayList<Batch> result=new ArrayList<Batch>();
		for(final Change locus:loci){
			for(int phase=0; phase<2; phase++){
				final int policy=phase==0 ? 0 : locus.policy;
				final long key=((long)locus.input<<32)|(policy&0xffffffffL);
				Batch batch=byKey.get(key);
				if(batch==null){
					batch=new Batch(locus.input, policy, prefix);
					byKey.put(key, batch);
					result.add(batch);
				}
				batch.loci.add(locus);
			}
		}
		return result;
	}

	/** Check every selected source pair before creating any output, including later batches. */
	private static void preflight(final ArrayList<Batch> batches){
		assert(batches!=null) : "Evidence file checks require the complete selected input/policy population";
		for(final Batch batch:batches){
			if(!Tools.testInputFiles(false, false, batch.stagePath, batch.gffPath)){
				throw new IllegalArgumentException("Cannot read native evidence files: "+batch.stagePath+" and "+batch.gffPath);
			}
		}
	}

	/** SweepNcrnaShortlist.writeReports emits every input in ordinal order, one row per policy. */
	private ObjectIntMap<String> bindInputs(final NcrnaJointPolicy[] roster){
		assert(roster.length>0) : "Native ordinal binding needs the exact nonempty policy roster";
		final ObjectIntMap<String> index=new ObjectIntMap<String>(String.class);
		final byte[][] policyTokens=new byte[roster.length][];
		for(int i=0; i<roster.length; i++){
			policyTokens[i]=BuildNcrnaReference.field(new ByteBuilder(), roster[i].name).toBytes();
		}
		String id=null;
		byte[] idToken=null;
		int policy=0, ordinal=-1;
		try(Input in=new Input(inputs)){
			final byte[] header=in.reader.nextLine();
			require(header!=null, "Missing native input report header");
			p.set(header);
			require(p.termEquals("input_id", 0) && p.termEquals("policy", 1), "Native report must retain input/policy identity columns");
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()>=2 && matches(p, 1, policyTokens[policy]), "Missing, repeated or out-of-order native policy row");
				if(policy==0){
					id=text(0);
					idToken=p.parseByteArray(0);
					require(!index.contains(id), "Input ordinal cannot repeat: "+id);
					index.put(id, ++ordinal);
				}else{require(matches(p, 0, idToken), "One input must finish all policy rows before another input begins");}
				policy=(policy+1)%roster.length;
			}
		}
		require(policy==0 && !index.isEmpty(), "Native input report ended before a complete input/policy population");
		return index;
	}

	private ArrayList<Change> readChanges(final ObjectIntMap<String> inputIndex, final ObjectIntMap<String> policyIndex){
		assert(inputIndex!=null && policyIndex!=null) : "Changed loci require exact native sidecar bindings";
		final ArrayList<Change> result=new ArrayList<Change>();
		try(Input in=new Input(changes)){
			RrnaResourceIO.header(in.reader, "input_id\tpolicy\tchange\treference_id\tseqid\tstart\tstop\tstrand", changes);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				p.set(row);
				require(p.terms()==8, "Changed reference row needs eight fields");
				final String id=text(0), policy=text(1);
				final int input=inputIndex.get(id), pi=policyIndex.get(policy);
				require(input>=0 && pi>0, "Changed locus must belong to an observed nonbaseline policy: "+id+" "+policy);
				require(p.termEquals("LOST", 2) || p.termEquals("GAINED", 2), "Unknown paired change kind");
				final int start=p.parseInt(5), stop=p.parseInt(6);
				require(start>0 && stop>=start && (p.termEquals("+", 7) || p.termEquals("-", 7)), "Reference interval must be positive closed and stranded");
				p.setBounds(4);
				final String seqid=decoder.decode(row, p.a(), p.b());
				result.add(new Change(input, pi, line, row.clone(), seqid, start, stop));
			}
		}
		return result;
	}

	/** All target-family events on the same contig, including opposite strands and pre-DP accepts. */
	private void stages(final Output out, final Batch batch, final String policy, final String path, final long[] rows){
		assert(batch!=null && !batch.loci.isEmpty()) : "Each stage read is bound to saved changes from one original input/policy";
		assert(rows.length==batch.loci.size()) : "Every inspected locus/policy context needs its own evidence-row counter";
		try(Input in=new Input(path)){
			RrnaResourceIO.header(in.reader, NcrnaStageDiagLogSink.HEADER, path);
			long line=1;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				p.set(row);
				if(p.terms()!=21){
					throw new IllegalArgumentException("Stage sidecar must retain 21 fields: "+path+" line "+line);
				}
				if(!p.termStartsWith(family, 1)){continue;}
				if(!includeModels && p.termEquals("MODEL", 0)){continue;}
				final int strand=p.parseInt(3), length=p.parseInt(4);
				require((strand==0 || strand==1) && length>0, "Stage rows need a strand and complete contig length");
				final boolean seed=p.termEquals("SEED", 0), pathRow=p.termEquals("PATH", 0);
				int start=-1, stop=-1, os=-1, oe=-1;
				if(!seed){
					require(pathRow || p.termEquals("WINDOW", 0) || p.termEquals("MODEL", 0), "Unknown stage event");
					final int a=p.parseInt(pathRow ? 16 : 6), b=p.parseInt(pathRow ? 17 : 7);
					start=genomicStart(a, b, strand, length);
					stop=genomicStop(a, b, strand, length);
					if(!p.termEquals("MODEL", 0)){
						final int oa=p.parseInt(16), ob=p.parseInt(17);
						if(oa>=0 && ob>=oa){
							os=genomicStart(oa, ob, strand, length);
							oe=genomicStop(oa, ob, strand, length);
						}
					}
				}
				for(int i=0; i<batch.loci.size(); i++){
					final Change locus=batch.loci.get(i);
					if(contig(p, 2, locus.contig) && (seed || near(locus, start, stop))){
						require(locus.stop<=length, "Reference interval must fit its own input's declared contig length");
						emit(out, batch, locus, policy, path, line, start, stop, os, oe, row);
						rows[i]++;
					}
				}
			}
		}
	}

	private void gff(final Output out, final Batch batch, final String policy, final String path, final long[] rows){
		assert(batch!=null && !batch.loci.isEmpty()) : "Every GFF is bound to the same original input/policy as its stage sidecar";
		assert(rows.length==batch.loci.size()) : "Stage and GFF readers share one counter per original locus/policy context";
		try(Input in=new Input(path)){
			long line=0;
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				line++;
				if(row.length==0 || row[0]=='#'){continue;}
				p.set(row);
				require(p.terms()==9, "Expected nine-field final GFF");
				final int start=p.parseInt(3), stop=p.parseInt(4);
				require(start>0 && stop>=start, "Final GFF coordinates must be positive closed");
				for(int i=0; i<batch.loci.size(); i++){
					final Change locus=batch.loci.get(i);
					if(contig(p, 0, locus.contig) && near(locus, start, stop)){
						emit(out, batch, locus, policy, path, line, start, stop, start, stop, row);
						rows[i]++;
					}
				}
			}
		}
	}

	/** Reference padding is inspection context only; no grading decision is made. */
	private boolean near(final Change locus, final int start, final int stop){
		assert(start>0 && stop>=start) : "Evidence interval conversion must precede the closed-overlap context test";
		return start<=(long)locus.stop+padding && stop>=(long)locus.start-padding;
	}

	private void emit(final Output out, final Batch batch, final Change locus, final String policy, final String path, final long line,
			final int start, final int stop, final int os, final int oe, final byte[] row){
		assert(row==null ? path==null && line==0 : path!=null && line>0)
			: "Saved rows need a physical path/line; explicit empty-context witnesses have neither";
		buffer.clear().append(locus.line).tab().append(locus.row).tab();
		BuildNcrnaReference.field(buffer, policy).tab();
		BuildNcrnaReference.field(buffer, batch.stagePath).tab();
		BuildNcrnaReference.field(buffer, batch.gffPath).tab();
		BuildNcrnaReference.field(buffer, path==null ? "." : path).tab().append(line).tab().append(start).tab().append(stop).tab().append(os).tab().append(oe).tab();
		buffer.append(row==null ? "NO_SAVED_EVIDENCE" : "SAVED_ROW").tab();
		if(row==null){buffer.append('.');}
		else{BuildNcrnaReference.field(buffer, row, 0, row.length);}
		buffer.nl();
		out.writer.print(buffer);
	}

	/*--------------------------------------------------------------*/
	/*----------------      Coordinate Helpers      ----------------*/
	/*--------------------------------------------------------------*/

	/** Zero-based inclusive oriented coordinates to one-based inclusive genomic bounds. */
	public static int genomicStart(final int a, final int b, final int strand, final int length){
		require(a>=0 && b>=a && b<length && (strand==0 || strand==1), "Oriented closed interval must fit the declared contig");
		return strand==0 ? a+1 : length-b;
	}
	public static int genomicStop(final int a, final int b, final int strand, final int length){
		require(a>=0 && b>=a && b<length && (strand==0 || strand==1), "Oriented closed interval must fit the declared contig");
		return strand==0 ? b+1 : length-a;
	}

	/** Stage names may retain a description; GFF uses the first whitespace token. */
	private static boolean contig(final LineParser1 parser, final int column, final byte[] expected){
		assert(expected!=null && expected.length>0) : "Changed reference needs its literal first-token contig";
		parser.setBounds(column);
		final byte[] row=parser.line();
		int stop=parser.a();
		while(stop<parser.b() && (row[stop]&255)>32){stop++;}
		if(stop-parser.a()!=expected.length){return false;}
		for(int i=0; i<expected.length; i++){if(row[parser.a()+i]!=expected[i]){return false;}}
		return true;
	}

	/** Native metadata keys are UTF-8 and percent-framed by BuildNcrnaReference.field. */
	private String text(final int column){
		p.setBounds(column);
		return decoder.decode(p.line(), p.a(), p.b());
	}

	/** Compare the writer's encoded field bytes without decoding per-policy rows. */
	private static boolean matches(final LineParser1 parser, final int column, final byte[] token){
		assert(token!=null) : "Encoded metadata tokens must be established before comparison";
		parser.setBounds(column);
		if(parser.b()-parser.a()!=token.length){return false;}
		final byte[] row=parser.line();
		for(int i=0; i<token.length; i++){if(row[parser.a()+i]!=token[i]){return false;}}
		return true;
	}
	private static void require(final boolean condition, final String message){
		if(!condition){throw new IllegalArgumentException(message);}
	}

	private static final class Change{
		Change(final int input_, final int policy_, final long line_, final byte[] row_, final String contig_, final int start_, final int stop_){
			assert(input_>=0 && policy_>0 && start_>0 && stop_>=start_) : "Change identity and reference span were validated before retention";
			input=input_;
			policy=policy_;
			line=line_;
			row=row_;
			contig=contig_.getBytes(StandardCharsets.UTF_8);
			start=start_;
			stop=stop_;
		}
		final int input, policy, start, stop;
		final long line;
		final byte[] row, contig;
		/** Total saved rows across this change's baseline and candidate contexts. */
		long evidenceRows=0;
	}

	/** One retained original input/policy file pair, with all its requested contexts. */
	private static final class Batch{
		Batch(final int input_, final int policy_, final String prefix){
			assert(input_>=0 && policy_>=0) : "Native file ordinals must come from validated input and policy bindings";
			input=input_;
			policy=policy_;
			final String stem=prefix+".input"+input+".policy"+policy;
			stagePath=stem+".stages.tsv";
			gffPath=stem+".gff";
		}
		final int input, policy;
		final String stagePath, gffPath;
		final ArrayList<Change> loci=new ArrayList<Change>();
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	private final LineParser1 p=new LineParser1('\t');
	private final NcrnaHistoricalEvidence.Attributes decoder=new NcrnaHistoricalEvidence.Attributes();
	private final ByteBuilder buffer=new ByteBuilder();
	private PrintStream outstream;
	private String changes, inputs, policies, calls, family, output;
	private int padding=0;
	private boolean includeModels=false;
}

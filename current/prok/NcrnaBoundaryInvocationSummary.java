package prok;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;

import map.IntHashSet;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Tools;
import structures.ByteBuilder;

/** Summarizes observed custom-refiner calls, without scoring, matching or replay.
 * Forked from template.A_SampleByteFile with the shared native text owners.
 * One reader reuses its parser and primitive counters; only distinct family names
 * are retained. Candidate events, returned events and final emissions stay separate.
 * Unknown provenance is never interpreted as a known non-invocation.
 * @author Brian Bushnell, Raiden
 */
public final class NcrnaBoundaryInvocationSummary{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, NcrnaBoundaryInvocationSummary.class, false);
		try{
			final Parser parser=new Parser();
			parser.overwrite=false;
			int expected=-1;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT), b=split.length>1 ? split[1] : null;
				if(a.equals("expectedinputs")){expected=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown invocation-summary argument: "+arg);}
			}
			require(parser.in2==null && parser.out2==null && !parser.overwrite && !parser.append
				&& parser.maxReads<0 && parser.samplerate==1, "Invocation counts require an unfiltered complete trace and fresh output");
			Shared.setThreads(1);
			summarize(parser.in1, parser.out1, expected);
			pp.outstream.println("NCRNA_BOUNDARY_INVOCATION_SUMMARY_PASS");
		}finally{Shared.closeStream(pp.outstream);}
	}

	/** Validates a complete invocation-only stream before opening the fresh summary. */
	public static void summarize(final String in, final String out, final int expectedInputs){
		require(in!=null && out!=null && expectedInputs>0
			&& Tools.testInputFiles(false, true, in) && Tools.testOutputFiles(false, false, false, out)
			&& Tools.testForDuplicateFiles(true, in, out), "Supply in=invocation_trace out=fresh_summary expectedinputs=N");
		final NcrnaBoundaryInvocationSummary tool=new NcrnaBoundaryInvocationSummary();
		tool.read(in);
		require(tool.inputs.size()==expectedInputs, "Invocation trace input count differs: expected="+expectedInputs+", observed="+tool.inputs.size());
		for(int i=0; i<tool.observed.length; i++){
			require(tool.observed[i]==tool.declared[i], "Candidate/contig conservation differs at field "+i+": "+tool.observed[i]+" versus "+tool.declared[i]);
		}
		tool.rows.sort(Comparator.comparing(row->row.name));
		try(Output output=new Output(out)){
			output.writer.println(HEADER);
			final ByteBuilder b=new ByteBuilder();
			for(final Member row:tool.rows){
				b.clear().append(row.name);
				for(final long count:row.counts){b.tab().append(count);}
				output.writer.print(b.nl());
			}
		}
	}

	private NcrnaBoundaryInvocationSummary(){}

	/*--------------------------------------------------------------*/
	/*----------------           Reading            ----------------*/
	/*--------------------------------------------------------------*/

	/** Borrow the producer's explicit input IDs and compare candidate totals with T rows. */
	private void read(final String path){
		try(Input input=new Input(path)){
			final byte[] version=input.reader.nextLine();
			require(version!=null, "Invocation trace requires a version header");
			p.set(version);
			require(p.terms()==2 && p.termEquals("#ncrna_boundary_invocations", 0)
				&& (p.termEquals("1", 1) || p.termEquals("2", 1) || p.termEquals("3", 1)), "Unsupported invocation-trace version");
			final boolean legacyMaximum=!p.termEquals("3", 1);
			for(byte[] line; (line=input.reader.nextLine())!=null;){
				if(line.length==0 || line[0]=='#'){continue;}
				p.set(line);
				final int ordinal=p.parseInt(1);
				require(ordinal>=0, "Trace input ordinals must be nonnegative");
				if(p.termEquals("I", 0)){
					require(p.terms()==3 && inputs.add(ordinal), "Each trace input identity must be declared exactly once");
				}else{
					require(inputs.contains(ordinal), "Candidate/contig trace needs its preceding input declaration: "+ordinal);
					if(p.termEquals("C", 0)){candidate(legacyMaximum);}
					else if(p.termEquals("T", 0)){
						require(p.terms()==11, "Invocation contig totals require the declared eleven fields");
						for(int i=0; i<declared.length; i++){
							final long count=p.parseLong(i+6);
							require(count>=0, "Contig conservation counts cannot be negative");
							declared[i]=Math.addExact(declared[i], count);
						}
						require(p.parseLong(10)==0, "Invocation-only mode does not report score-window availability");
					}else{throw new IllegalArgumentException("Invocation-only trace cannot contain scores or unknown record types");}
				}
			}
		}
	}

	/** Counts actual call flags; no coordinate movement or neural eligibility inference. */
	private void candidate(final boolean legacyMaximum){
		require(p.terms()==37 && p.parseInt(29)==0 && p.parseInt(30)==0, "Invocation candidates have37fields and no score-replay radius/mask");
		require(NcrnaTextIO.maximumLength(p, 22, legacyMaximum)>=p.parseInt(21), "Invocation length bounds must remain ordered");
		final int known=flag(14), invoked=flag(16), returned=flag(31), selected=flag(32);
		require(invoked<=known && selected<=returned, "Unknown provenance cannot claim invocation; every selected call must have returned");
		require(p.termEquals(selected==1 ? "SELECTED" : returned==1 ? "PATH_REJECTED" : "SCAVENGER_REJECTED", 33),
			"Candidate disposition must agree with the returned/selected flags");
		Member row=null;
		for(final Member old:rows){if(p.termEquals(old.name, 7)){row=old;break;}}
		if(row==null){
			final String name=p.parseString(7);
			NcrnaFamilyRegistry.sharedBoundaryMember(name);
			row=new Member(name);
			rows.add(row);
		}
		final long[] c=row.counts;
		c[0]++;
		c[1]+=known;
		c[2]+=invoked;
		c[3]+=invoked*returned;
		c[4]+=invoked*selected;
		c[5]+=selected;
		c[6]+=(1-known)*selected;
		observed[0]++;
		observed[1]+=returned;
		observed[2]+=selected;
		observed[3]+=1-known;
	}

	private int flag(final int column){
		final int value=p.parseInt(column);
		require(p.termEquals("0", column) || p.termEquals("1", column),
			"Invocation trace flags must be literal0or1 at column "+column);
		return value;
	}
	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------       Retained Metadata      ----------------*/
	/*--------------------------------------------------------------*/

	private static final class Member{
		Member(final String name_){name=name_;}
		final String name;
		final long[] counts=new long[7];
	}
	private final LineParser1 p=new LineParser1('\t');
	private final IntHashSet inputs=new IntHashSet();
	private final ArrayList<Member> rows=new ArrayList<Member>();
	/** Candidate, returned, selected and unknown, independently declared by contig totals. */
	private final long[] observed=new long[4], declared=new long[4];
	public static final String HEADER="member\tcandidates\tknown\tcustom_refiner_invocations\treturned_invocations\tselected_invocations\tselected_calls\tunknown_selected";
}

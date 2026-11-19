package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;

import fileIO.ReadWrite;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import shared.Timer;
import shared.Tools;
import structures.ByteBuilder;

/** Checks actual compact caller resources against accepted saved parent ranks.
 * Derived from template.A_SampleByteFile's parse/validate/process lifecycle.
 * Existing byte-stream owners bind the complete original cache and model roster.
 * Uses the caller loader/scorer; no member voting, new identities or alignments.
 * @author Brian, Raiden
 */
public final class VerifyNcrnaDdlProfiles{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final Timer timer=new Timer();
		final PreParser pp=new PreParser(args, VerifyNcrnaDdlProfiles.class, false);
		final int threads=Shared.threads();
		try{
			Shared.setThreads(1);
			final Config c=parse(pp.args);
			require(Shared.threads()==1, "Saved resource parity is a single-owner pass; use t=1");
			run(c);
			timer.stop();
			pp.outstream.println("Time: "+timer);
			pp.outstream.println("NCRNA_DDL_RESOURCE_SAVED_RANK_PARITY_PASS");
		}finally{
			Shared.setThreads(threads);
			Shared.closeStream(pp.outstream);
		}
	}

	private static Config parse(final String[] args){
		assert(args!=null) : "PreParser must supply the expanded native argument list";
		final Config c=new Config();
		final Parser parser=new Parser();
		parser.overwrite=false;
		for(final String arg:args){
			final String[] split=Parse.splitOnFirst(arg, '=');
			final String a=split[0].toLowerCase(Locale.ROOT);
			final String b=split.length>1 ? split[1] : null;
			if(a.equals("panels")){c.panels=b;}
			else if(a.equals("selection")){c.selection=b;}
			else if(a.equals("ref")){c.references=b;}
			else if(a.equals("ranks")){c.ranks=b;}
			else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
			else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
			else if(a.equals("indexparity")){
				throw new IllegalArgumentException("indexparity was a temporary migration diagnostic and has been retired; use ordinary saved-rank verification");
			}
			else{require(parser.parse(arg, a, b), "Unknown resource-parity parameter: "+arg);}
		}
		require(parser.out2==null && !parser.overwrite && !parser.append
			&& parser.maxReads<0 && parser.samplerate==1,
			"Supply complete resources/cache/ranks without sampling, overwrite or append");
		c.profile512=parser.in1;
		c.profile1024=parser.in2;
		c.output=parser.out1;
		return c;
	}

	/** Verifies every original parent before opening the small result file. */
	public static void run(final Config c){
		require(c!=null && c.profile512!=null && c.profile1024!=null && c.panels!=null
			&& c.selection!=null && c.references!=null && c.ranks!=null && c.output!=null
			&& c.models>0 && c.groups>0,
			"Bind both compact resources, original panels/selection/reference, saved ranks and positive expected counts");
		require(Tools.testInputFiles(false, true, c.profile512, c.profile1024, c.panels, c.selection, c.references, c.ranks)
			&& Tools.testOutputFiles(false, false, false, c.output), "Parity inputs must be readable and output fresh");
		require(Tools.testForDuplicateFiles(true, c.profile512, c.profile1024, c.panels, c.selection, c.references, c.ranks, c.output),
			"Distinct resource/cache/rank inputs cannot overwrite or alias the result path");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions, "Resource parity requires assertions in the native cache/model readers");
		final int zip=ReadWrite.MAX_ZIP_THREADS();
		final boolean zipSet=ReadWrite.SET_ZIP_THREADS, trim=Shared.TRIM_READ_DESCRIPTION;
		try{
			ReadWrite.setZipThreads(1);
			Shared.TRIM_READ_DESCRIPTION=false;
			new VerifyNcrnaDdlProfiles(c).process();
		}finally{
			ReadWrite.setZipThreads(zip);
			ReadWrite.SET_ZIP_THREADS=zipSet;
			Shared.TRIM_READ_DESCRIPTION=trim;
		}
	}

	private VerifyNcrnaDdlProfiles(final Config c_){c=c_;}

	/*--------------------------------------------------------------*/
	/*----------------          Processing          ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		final ArrayList<String> names=new ArrayList<String>();
		final ArrayList<byte[]> sequences=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			names.add(read.id);
			sequences.add(read.bases);
		});
		require(names.size()==c.models, "The reference must retain every original cached model");
		final String[] headers=names.toArray(new String[names.size()]);
		final byte[][] bases=sequences.toArray(new byte[sequences.size()][]);
		final NcrnaDdlProfileResource[] resources={
			NcrnaDdlProfileResource.load(c.profile512, headers, bases),
			NcrnaDdlProfileResource.load(c.profile1024, headers, bases)};
		final NcrnaDdlProfileResource.Scorer[] scorers=new NcrnaDdlProfileResource.Scorer[2];
		final int[][] ranks=new int[2][data.parents.size()];
		final float[][] scores=new float[2][c.models];
		for(int variant=0; variant<2; variant++){
			require(resources[variant].compact() && resources[variant].buckets()==WIDTHS[variant],
				"Resource slots must bind selected compact512 then compact1024, not exchanged or historical profiles");
			scorers[variant]=resources[variant].newScorer();
			Arrays.fill(ranks[variant], -1);
		}
		final boolean[] acceptable=new boolean[c.models];
		final int[] tuple=new int[4];
		for(final Window window:data.windows){
			for(int variant=0; variant<2; variant++){
				scorers[variant].score(window.bases, scores[variant]);
			}
			for(final Edge edge:window.edges){
				if(edge.unique==0){continue;}
				for(int model=0; model<c.models; model++){acceptable[model]=(edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0;}
				for(int variant=0; variant<2; variant++){
					NcrnaDdlRanks.best(scores[variant], acceptable, tuple);
					final int parent=edge.parent.index;
					ranks[variant][parent]=NcrnaDdlSavedData.minimum(ranks[variant][parent], tuple[0]);
				}
			}
		}
		compare(ranks);
		write(ranks);
	}

	/** Accepted report rows follow data.parents order within each method.
	 * Bind UTF-8 bytes directly; String termEquals is ASCII-oriented. */
	private void compare(final int[][] ranks){
		assert(ranks.length==2 && ranks[0].length==data.parents.size() && ranks[1].length==data.parents.size())
			: "Both loaded resource widths must retain the complete original target roster";
		final int[] next=new int[2];
		final LineParser1 p=new LineParser1('\t');
		try(Input in=new Input(c.ranks)){
			RrnaResourceIO.header(in.reader, RANK_HEADER, c.ranks);
			for(byte[] row; (row=in.reader.nextLine())!=null;){
				p.set(row);
				require(p.terms()==10, "Saved compact locus rows require their complete ten-column schema");
				final int variant=p.termEquals(METHODS[0], 0) ? 0 : p.termEquals(METHODS[1], 0) ? 1 : -1;
				if(variant<0){continue;}
				final int index=next[variant]++;
				require(index<data.parents.size(), "Saved compact method has extra or duplicate parent rows: "+METHODS[variant]);
				final Parent parent=data.parents.get(index);
				require(equal(p, 1, parent.input.panel.nameBytes) && p.parseInt(2)==parent.input.slot
					&& equal(p, 3, parent.input.idBytes) && equal(p, 4, parent.idBytes)
					&& p.termEquals(NcrnaDdlSavedData.STRATA[parent.stratum], 5)
					&& p.termEquals(parent.nativeRank>0 ? "true" : "false", 7) && p.termEquals(parent.reason, 9),
					"Saved rank rows must bind original panel/input/parent order and dispositions: method="+METHODS[variant]+", parent="+index);
				final int expected=p.parseInt(8), observed=ranks[variant][index];
				require(expected==observed, "Caller resource rank differs from saved rank: method="+METHODS[variant]
					+", parent="+index+", expected="+expected+", observed="+observed);
			}
		}
		for(int variant=0; variant<2; variant++){
			require(next[variant]==data.parents.size(), "Saved compact method omits original parents: "+METHODS[variant]);
		}
	}

	private static boolean equal(final LineParser1 p, final int term, final byte[] expected){
		assert(expected!=null && term>=0 && term<p.terms()) : "Byte comparison binds initialized original metadata to an existing parsed term";
		p.setBounds(term);
		if(p.b()-p.a()!=expected.length){return false;}
		final byte[] row=p.line();
		for(int i=0; i<expected.length; i++){
			if(row[p.a()+i]!=expected[i]){return false;}
		}
		return true;
	}

	private void write(final int[][] ranks){
		assert(ranks.length==2) : "Only the two selected compact caller widths have passed comparison";
		try(Output out=new Output(c.output)){
			out.writer.println("total_bins\tmodels\tinputs\twindows\tparents\tmatched_ranks\tavailable_ranks\tunavailable_ranks");
			final ByteBuilder b=new ByteBuilder();
			for(int variant=0; variant<2; variant++){
				int available=0;
				for(final int rank:ranks[variant]){if(rank>0){available++;}}
				out.writer.print(b.clear().append(WIDTHS[variant]).tab().append(c.models).tab().append(data.selected.size())
					.tab().append(data.windows.size()).tab().append(data.parents.size()).tab().append(ranks[variant].length)
					.tab().append(available).tab().append(ranks[variant].length-available).nl());
			}
		}
	}

	private static void require(final boolean valid, final String why){if(!valid){throw new IllegalArgumentException(why);}}

	/*--------------------------------------------------------------*/
	/*----------------             Fields           ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String profile512, profile1024, panels, selection, references, ranks, output;
		public int models=508, groups=3153;
	}
	private final Config c;
	private NcrnaDdlSavedData data;
	private static final int[] WIDTHS={512,1024};
	private static final String[] METHODS={"partition_plain7_masked10_512", "partition_plain7_masked10_1024"};
	private static final String RANK_HEADER="method\tpanel\tinput_slot\tinput_id\tparent\tseen_stratum\tidentity_bin\tnative_common\tbest_rank\tsource_disposition";
}

package prok;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Locale;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;

import cardinality.DynamicDemiLog;
import map.ObjectMap;
import parse.LineParser1;
import parse.Parse;
import parse.Parser;
import parse.PreParser;
import prok.NcrnaDdlSavedData.Edge;
import prok.NcrnaDdlSavedData.Parent;
import prok.NcrnaDdlSavedData.Selected;
import prok.NcrnaDdlSavedData.Window;
import prok.NcrnaTextIO.Input;
import prok.NcrnaTextIO.Output;
import shared.Shared;
import structures.ByteBuilder;
import template.ThreadWaiter;

/** Native-first duplicate-skipping hybrid on saved windows and verifier outcomes.
 * Native model counts are read from the original model-complete oracle table.
 * Only the selected sketch and original canonical controls are reconstructed.
 * Native/canonical event tuples and all original control minima must match;
 * the selected new sketch matches its saved locus minima, not an unavailable
 * independent event-tuple cache. No genome, caller, verifier or alignment runs.
 * Uses native ByteFile/Streamer input owners and the template worker lifecycle.
 * @author Brian, Raiden
 */
public final class HybridSavedNcrnaWindows{

	/*--------------------------------------------------------------*/
	/*----------------        Initialization        ----------------*/
	/*--------------------------------------------------------------*/

	public static void main(final String[] args){
		final PreParser pp=new PreParser(args, HybridSavedNcrnaWindows.class, false);
		final int previous=DynamicDemiLog.exponentBits();
		try{
			final Config c=new Config();
			final Parser parser=new Parser();
			parser.overwrite=false;
			for(final String arg:pp.args){
				final String[] split=Parse.splitOnFirst(arg, '=');
				final String a=split[0].toLowerCase(Locale.ROOT);
				final String b=split.length>1 ? split[1] : null;
				if(a.equals("panels")){c.panels=b;}
				else if(a.equals("selection")){c.selection=b;}
				else if(a.equals("grid")){c.grid=b;}
				else if(a.equals("method")){c.method=Parse.parseIntKMG(b);}
				else if(a.equals("members")){c.members=b;}
				else if(a.equals("ledger")){c.ledger=b;}
				else if(a.equals("assignments")){c.assignments=b;}
				else if(a.equals("census")){c.census=b;}
				else if(a.equals("summary")){c.summary=b;}
				else if(a.equals("ref")){c.references=b;}
				else if(a.equals("expectedgroups")){c.groups=Parse.parseIntKMG(b);}
				else if(a.equals("expectedmodels")){c.models=Parse.parseIntKMG(b);}
				else if(a.equals("t") || a.equals("threads")){c.threads=Parse.parseIntKMG(b);}
				else{require(parser.parse(arg, a, b), "Unknown saved-hybrid parameter: "+arg);}
			}
			require(parser.out1!=null && parser.in1==null && parser.in2==null && parser.out2==null
				&& !parser.overwrite && !parser.append && parser.maxReads<0 && parser.samplerate==1,
				"Use grid= panels= selection= out=fresh_prefix; no genomes, overwrite, append or sampling");
			c.output=parser.out1;
			DynamicDemiLog.setExponent(4);
			run(c);
			pp.outstream.println("NCRNA_SAVED_NATIVE_FIRST_HYBRID_COMPLETE_NOT_CALLER_ADOPTED");
		}finally{
			DynamicDemiLog.setExponent(previous);
			Shared.closeStream(pp.outstream);
		}
	}

	public static void run(final Config c){
		require(c!=null && c.panels!=null && c.selection!=null && c.grid!=null && c.output!=null
			&& c.members!=null && c.ledger!=null && c.assignments!=null && c.census!=null
			&& c.summary!=null && c.references!=null && c.method>0 && c.models>0 && c.groups>0 && c.threads>0,
			"Declare the selected saved method, complete cache/member bindings and positive model/group/worker counts");
		boolean assertions=false;
		assert(assertions=true);
		require(assertions && DynamicDemiLog.exponentBits()==4 && !DynamicDemiLog.blacklistExists(),
			"Saved hybrids require assertions, native exponent4 and no DDL blacklist");
		RrnaResourceIO.fresh(c.output, ".methods.tsv", ".summary.tsv", ".loci.tsv.gz", ".input_summary.tsv", ".last_losses.tsv", ".work.tsv");
		Shared.TRIM_READ_DESCRIPTION=false;
		new HybridSavedNcrnaWindows(c).process();
	}
	private HybridSavedNcrnaWindows(final Config c_){c=c_;}

	/*--------------------------------------------------------------*/
	/*----------------        Owner Prepasses       ----------------*/
	/*--------------------------------------------------------------*/

	private void process(){
		method();
		data=NcrnaDdlSavedData.load(c.panels, c.selection, c.models, c.groups);
		selectedMinima();
		final NcrnaDdlMembers members=members();
		data.readNativeCounts(names);
		windows=data.windows.toArray(new Window[data.windows.size()]);
		profiles=new NcrnaDdlProfile[9][];
		profiles[0]=NcrnaDdlModelProfiles.build(members, k, width, false, middle);
		for(int g=1; g<profiles.length; g++){
			profiles[g]=NcrnaDdlModelProfiles.build(members, g<=4 ? 9 : 11, WIDTHS[(g-1)&3], true, false);
		}
		final int workers=Math.min(c.threads, Math.max(1, windows.length));
		final ArrayList<ProcessThread> threads=new ArrayList<ProcessThread>(workers);
		for(int i=0; i<workers; i++){threads.add(new ProcessThread(this, i));}
		next.set(Math.min(workers, windows.length));
		ThreadWaiter.startAndWait(threads);
		if(failure.get()!=null){throw new IllegalStateException("Saved hybrid failed before result publication", failure.get());}
		minima=new int[data.parents.size()*TOTAL];
		Arrays.fill(minima, -1);
		for(final ProcessThread worker:threads){
			require(worker.success, "Every hybrid worker must finish before its minima are used");
			for(int i=0; i<minima.length; i++){minima[i]=NcrnaDdlSavedData.minimum(minima[i], worker.minima[i]);}
			scored+=worker.scored;
			nativeTuples+=worker.nativeTuples;
			canonicalTuples+=worker.canonicalTuples;
			if(worker.scored>0){scoringWorkers++;}
		}
		for(final Parent parent:data.parents){
			final int offset=parent.index*TOTAL;
			sameMinimum(parent, minima[offset], parent.nativeRank, "Cached native event tuples must conserve every original locus minimum");
			sameMinimum(parent, minima[offset+1], selected[parent.index], "The selected new sketch must reproduce every saved grid locus minimum");
			for(int control=0; control<NcrnaDdlSavedData.CONTROLS; control++){
				sameMinimum(parent, minima[offset+METHODS+control], parent.controlMin[control], "Original canonical event controls must conserve all saved locus minima");
			}
		}
		report(workers);
	}

	/** Parameters come from the actual selected row, preventing nominal flags
	 * from disagreeing with the method whose saved minima are being checked. */
	private void method(){
		boolean found=false, baseline=false;
		try(Input in=new Input(c.grid+".methods.tsv")){
			NcrnaDdlSavedData.header(in, "method\tk\tbuckets\torientation\tmiddle_mask\tscore\trole");
			final LineParser1 p=new LineParser1('\t');
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==7, "Selected grid methods require the original seven-column schema");
				if(p.parseInt(0)==0){
					require(!baseline && p.parseInt(2)==0 && p.termEquals("FORWARD_NATIVE_INDEX", 3)
						&& p.termEquals("false", 4) && p.termEquals("POSITIVE_SHARED_KEYS", 5) && p.termEquals("SAVED_BASELINE", 6),
						"The saved comparison needs exactly one native positive-count baseline");
					baseline=true;
					nativeK=p.parseInt(1);
				}
				if(p.parseInt(0)!=c.method){continue;}
				require(!found && p.termEquals("FORWARD", 3) && p.termEquals("NEW_GRID", 6)
					&& (p.termEquals("true", 4) || p.termEquals("false", 4)), "Select exactly one new forward grid method");
				found=true;
				k=p.parseInt(1);
				width=p.parseInt(2);
				middle=p.termEquals("true", 4);
				score=p.termEquals("SSU", 5) ? NcrnaDdlRanks.SSU : p.termEquals("RAW_ALL_QUERY_KEYS", 5) ? NcrnaDdlRanks.ALL_KEYS : -1;
			}
		}
		require(found && baseline && nativeK>0 && nativeK<=15 && k>=7 && k<=11 && width>=256 && width<=2048
			&& Integer.bitCount(width)==1 && score>=0, "The chosen method must be one supported forward k7–11/width256–2048 raw or SSU geometry");
	}

	private void selectedMinima(){
		selected=new int[data.parents.size()];
		final boolean[][] seen=new boolean[2][selected.length];
		final ObjectMap<String, Parent> parents=new ObjectMap<String, Parent>(String.class, Parent.class);
		for(final Parent parent:data.parents){parents.put(parent.input.panel.name+'\t'+new String(parent.idBytes, java.nio.charset.StandardCharsets.UTF_8), parent);}
		try(Input in=new Input(c.grid+".loci.tsv.gz")){
			NcrnaDdlSavedData.header(in, "panel\tinput_slot\tinput_id\tparent\tseen_stratum\tmethod\tbest_rank\tsource_disposition\tinterpretation");
			final LineParser1 p=new LineParser1('\t');
			for(byte[] row=in.reader.nextLine(); row!=null; row=in.reader.nextLine()){
				p.set(row);
				require(p.terms()==9, "Saved selected-method minima require the complete grid locus schema");
				final int method=p.parseInt(5);
				if(method!=0 && method!=c.method){continue;}
				final Parent parent=parents.get(NcrnaDdlSavedData.text(p, 0)+'\t'+NcrnaDdlSavedData.text(p, 3));
				require(parent!=null && p.parseInt(1)==parent.input.slot
					&& NcrnaDdlSavedData.text(p, 2).equals(parent.input.id)
					&& p.termEquals(NcrnaDdlSavedData.STRATA[parent.stratum], 4) && p.termEquals(parent.reason, 7)
					&& p.termEquals("SAVED_LOCUS_MINIMUM_NOT_CALLER_WORK", 8), "Grid minima must bind the exact original input, parent, seen label and disposition");
				final int column=method==0 ? 0 : 1, rank=p.parseInt(6);
				require(!seen[column][parent.index] && (rank==-1 || rank>0 && rank<=c.models), "Each selected parent needs one valid native and chosen-method minimum");
				seen[column][parent.index]=true;
				if(column==0){require(rank==parent.nativeRank, "The grid native minimum cannot exchange the original baseline");}
				else{selected[parent.index]=rank;}
			}
		}
		for(int parent=0; parent<selected.length; parent++){
			require(seen[0][parent] && seen[1][parent], "Saved grid minima cannot omit unrankable or zero-native selected targets");
		}
	}

	private NcrnaDdlMembers members(){
		final ArrayList<String> labels=new ArrayList<String>();
		final ArrayList<byte[]> bases=new ArrayList<byte[]>();
		NcrnaDdlSavedData.sequences(c.references, read->{
			labels.add(read.id);
			bases.add(read.bases);
		});
		require(labels.size()==c.models, "The cached native model ordinals need the complete unchanged reference roster");
		names=labels.toArray(new String[labels.size()]);
		return NcrnaDdlMembers.load(c.members, c.ledger, c.assignments, c.census, c.summary, names, bases.toArray(new byte[bases.size()][]));
	}

	/*--------------------------------------------------------------*/
	/*----------------         Query Workers        ----------------*/
	/*--------------------------------------------------------------*/

	static final class ProcessThread extends Thread{
		ProcessThread(final HybridSavedNcrnaWindows owner_, final int index_){
			assert(owner_.profiles.length==9 && owner_.profiles[0].length==owner_.c.models)
				: "Workers index one selected geometry and eight canonical controls in the complete model order";
			owner=owner_;
			index=index_;
			minima=new int[owner.data.parents.size()*TOTAL];
			Arrays.fill(minima, -1);
			ranks[0]=new NcrnaDdlRanks(owner.profiles[0], owner.k, false, owner.middle);
			for(int g=1; g<ranks.length; g++){ranks[g]=new NcrnaDdlRanks(owner.profiles[g], g<=4 ? 9 : 11);}
			order=new NcrnaDdlInterleave(owner.c.models);
			scores=new float[owner.c.models];
			acceptable=new boolean[owner.c.models];
		}
		@Override
		public void run(){
			try{
				if(index<owner.windows.length){window(owner.windows[index]);}
				for(int first=owner.next.getAndAdd(8); first<owner.windows.length; first=owner.next.getAndAdd(8)){
					if(owner.failure.get()!=null){return;}
					for(int i=first; i<Math.min(first+8, owner.windows.length); i++){window(owner.windows[i]);}
				}
				success=true;
			}catch(RuntimeException | Error invalid){owner.failure.compareAndSet(null, invalid);}
		}
		private void window(final Window window){
			assert(window.nativeCounts!=null && window.nativeCounts.length==acceptable.length)
				: "The native-count prepass must finish every selected window before query workers start";
			boolean useful=false;
			for(final Edge edge:window.edges){useful|=edge.unique>0;}
			if(useful){
				for(final NcrnaDdlRanks rank:ranks){rank.score(window.bases);}
				for(int model=0; model<scores.length; model++){scores[model]=ranks[0].score(owner.score, model);}
				order.build(window.nativeCounts, scores);
				scored++;
			}
			for(final Edge edge:window.edges){
				final int offset=edge.parent.index*TOTAL;
				for(int model=0; model<acceptable.length; model++){
					acceptable[model]=(edge.uniqueModels[model>>>6]&(1L<<(model&63)))!=0;
				}
				if(edge.unique>0){order.nativeBest(acceptable, tuple);}else{Arrays.fill(tuple, -1);}
				if(!Arrays.equals(tuple, edge.nativeTuple)){
					throw new IllegalArgumentException("Cached native rank/tie/model parity failed: panel="+window.input.panel.name
						+", event="+window.key+", parent="+new String(edge.parent.idBytes, java.nio.charset.StandardCharsets.UTF_8)
						+", expected="+Arrays.toString(edge.nativeTuple)+", observed="+Arrays.toString(tuple));
				}
				nativeTuples++;
				minima[offset]=NcrnaDdlSavedData.minimum(minima[offset], tuple[0]);
				if(edge.unique>0){order.sketchBest(acceptable, tuple);}else{Arrays.fill(tuple, -1);}
				minima[offset+1]=NcrnaDdlSavedData.minimum(minima[offset+1], tuple[0]);
				if(edge.unique>0){order.hybridBest(acceptable, tuple);}else{Arrays.fill(tuple, -1);}
				minima[offset+2]=NcrnaDdlSavedData.minimum(minima[offset+2], tuple[0]);
				for(int g=1; g<ranks.length; g++){
					for(int s=0; s<2; s++){
						final int control=(g-1)*2+s;
						if(edge.unique>0){ranks[g].best(SCORES[s], acceptable, tuple);}else{Arrays.fill(tuple, -1);}
						for(int field=0; field<4; field++){
							if(tuple[field]!=edge.controls[control*4+field]){
							throw new IllegalArgumentException("Canonical hybrid-control rank/tie/model parity failed: panel="+window.input.panel.name
								+", event="+window.key+", control="+control+", field="+field+", expected="+edge.controls[control*4+field]+", observed="+tuple[field]);
							}
						}
						canonicalTuples++;
						minima[offset+METHODS+control]=NcrnaDdlSavedData.minimum(minima[offset+METHODS+control], tuple[0]);
					}
				}
			}
		}
		private final HybridSavedNcrnaWindows owner;
		private final int index;
		private final NcrnaDdlRanks[] ranks=new NcrnaDdlRanks[9];
		private final NcrnaDdlInterleave order;
		private final float[] scores;
		private final boolean[] acceptable;
		private final int[] tuple=new int[4];
		final int[] minima;
		boolean success;
		long scored, nativeTuples, canonicalTuples;
	}

	/*--------------------------------------------------------------*/
	/*----------------            Reports           ----------------*/
	/*--------------------------------------------------------------*/

	private void report(final int workers){
		final long[][][] hist=new long[3][METHODS][c.models+1];
		final long[][] sums=new long[3][METHODS];
		final int[][] better=new int[3][METHODS], tied=new int[3][METHODS], worse=new int[3][METHODS];
		final int[] common=new int[3];
		final boolean[][] reachable=new boolean[METHODS][data.selected.size()];
		try(Output loci=new Output(c.output+".loci.tsv.gz"); Output stats=new Output(c.output+".summary.tsv");
			Output losses=new Output(c.output+".last_losses.tsv"); Output inputs=new Output(c.output+".input_summary.tsv");
			Output methods=new Output(c.output+".methods.tsv"); Output work=new Output(c.output+".work.tsv")){
			final ByteBuilder b=new ByteBuilder();
			methods.writer.println("method\tk\tbuckets\torientation\tmiddle_mask\tscore\tinterpretation");
			methods.writer.print(b.clear().append("NATIVE\t").append(nativeK).append("\t0\tFORWARD_NATIVE_INDEX\tfalse\tPOSITIVE_SHARED_KEYS\tSAVED_BASELINE\n"));
			for(int method=1; method<METHODS; method++){
				methods.writer.print(b.clear().append(METHOD_NAMES[method]).tab().append(k).tab().append(width)
					.append("\tFORWARD\t").append(middle).tab().append(score==NcrnaDdlRanks.SSU ? "SSU" : "RAW_ALL_QUERY_KEYS")
					.tab().append(method==1 ? "SAVED_GRID_SELECTED_METHOD" : "NATIVE_FIRST_PROPOSAL_ZIP_SKIP_REPEATS").nl());
			}
			loci.writer.println("panel\tinput_slot\tinput_id\tparent\tseen_stratum\tnative_rank\tsketch_rank\thybrid_rank\tsource_disposition");
			for(final Parent parent:data.parents){
				final int offset=parent.index*TOTAL;
				if(parent.nativeRank>0){common[parent.stratum]++;}
				for(int method=0; method<METHODS; method++){
					final int rank=minima[offset+method];
					if(rank>0 && rank<=50){reachable[method][parent.input.index]=true;}
					if(parent.nativeRank<0){continue;}
					require(rank>0 && rank<=c.models, "A common locus must retain an actual model position in all three complete orders");
					hist[parent.stratum][method][rank]++;
					sums[parent.stratum][method]+=rank;
					if(rank<parent.nativeRank){better[parent.stratum][method]++;}
					else if(rank==parent.nativeRank){tied[parent.stratum][method]++;}
					else{worse[parent.stratum][method]++;}
				}
				loci.writer.print(b.clear().append(parent.input.panel.nameBytes).tab().append(parent.input.slot).tab()
					.append(parent.input.idBytes).tab().append(parent.idBytes).tab().append(NcrnaDdlSavedData.STRATA[parent.stratum])
					.tab().append(minima[offset]).tab().append(minima[offset+1]).tab().append(minima[offset+2]).tab().append(parent.reason).nl());
			}
			stats.writer.println("seen_stratum\tmethod\tcommon_loci\trank_sum\tmean_rank\twithin10\twithin50\twithin100\tp90\tp95\tp99\tbetter\ttied\tworse\tinterpretation");
			for(int stratum=0; stratum<3; stratum++){
				for(int method=0; method<METHODS; method++){
					final long[] h=hist[stratum][method];
					b.clear().append(NcrnaDdlSavedData.STRATA[stratum]).tab().append(METHOD_NAMES[method]).tab().append(common[stratum])
						.tab().append(sums[stratum][method]).tab();
					if(common[stratum]==0){b.append("NA");}else{b.append((double)sums[stratum][method]/common[stratum], 9);}
					stats.writer.print(b.tab().append(covered(h, 10)).tab().append(covered(h, 50)).tab().append(covered(h, 100))
						.tab().append(RankNcrnaWindows.quantile(h, .9)).tab().append(RankNcrnaWindows.quantile(h, .95))
						.tab().append(RankNcrnaWindows.quantile(h, .99)).tab().append(better[stratum][method])
						.tab().append(tied[stratum][method]).tab().append(worse[stratum][method]).append("\tSAVED_LOCUS_MINIMA_NOT_CALLER_WORK\n"));
				}
			}
			losses.writer.println("method\tpanel\tinput_slot\tinput_id\tfasta\tinterpretation");
			inputs.writer.println("method\tselected_input_records\tnative_reachable_input_records_at50\treachable_input_records_at50\tlast_locus_losses_at50\tinterpretation");
			int baseline=0;
			for(final boolean reached:reachable[0]){if(reached){baseline++;}}
			for(int method=0; method<METHODS; method++){
				int lost=0, reached=0;
				for(final Selected input:data.selected){
					if(reachable[method][input.index]){reached++;}
					if(!reachable[0][input.index] || reachable[method][input.index]){continue;}
					lost++;
					losses.writer.print(b.clear().append(METHOD_NAMES[method]).tab().append(input.panel.nameBytes).tab().append(input.slot)
						.tab().append(input.idBytes).tab().append(input.pathBytes).append("\tLAST_SAVED_REACHABLE_LOCUS_LOST_AT50_NOT_CALLER_DETECTION\n"));
				}
				inputs.writer.print(b.clear().append(METHOD_NAMES[method]).tab().append(data.selected.size()).tab().append(baseline)
					.tab().append(reached).tab().append(lost).append("\tINPUT_RECORDS_NOT_ESTABLISHED_INDEPENDENT_SPECIES\n"));
			}
			work.writer.print(b.clear().append("item\tvalue\nselected_inputs\t").append(data.selected.size())
				.append("\nselected_targets\t").append(data.parents.size()).append("\nselected_windows\t").append(windows.length)
				.append("\ncached_native_model_rows\t").append((long)windows.length*c.models)
				.append("\nscored_windows\t").append(scored).append("\nworkers_started\t").append(workers)
				.append("\nworkers_with_scored_windows\t").append(scoringWorkers).append("\nselected_grid_method\t").append(c.method)
				.append("\nnative_event_tuples_checked\t").append(nativeTuples).append("\ncanonical_event_tuples_checked\t").append(canonicalTuples)
				.append("\nselected_sketch_parity\tSAVED_LOCUS_MINIMA_ONLY_NO_INDEPENDENT_NEW_EVENT_TUPLE_CACHE")
				.append("\nhybrid_rule\tNATIVE_I_THEN_SKETCH_I_DUPLICATE_PROPOSALS_DO_NOT_EMIT_ANOTHER_RANK")
				.append("\nresult\tSAVED_ORDER_DIAGNOSTIC_NOT_CALLER_DETECTION_SPEED_OR_ADOPTION\n"));
		}
	}
	private static long covered(final long[] h, final int limit){
		assert(h!=null && h.length>0 && limit>0) : "Rank coverage scans positive histogram bins for the declared cutoff "+limit;
		long count=0;
		for(int rank=1; rank<=limit && rank<h.length; rank++){count+=h[rank];}
		return count;
	}
	private static void require(final boolean valid, final String why){NcrnaDdlSavedData.require(valid, why);}
	private static void sameMinimum(final Parent parent, final int observed, final int expected, final String why){
		if(observed!=expected){
			throw new IllegalArgumentException(why+": panel="+parent.input.panel.name+", parent="
				+new String(parent.idBytes, java.nio.charset.StandardCharsets.UTF_8)+", expected="+expected+", observed="+observed);
		}
	}

	/*--------------------------------------------------------------*/
	/*----------------            Fields            ----------------*/
	/*--------------------------------------------------------------*/

	public static final class Config{
		public String panels, selection, grid, output, members, ledger, assignments, census, summary, references;
		public int method=-1, models=508, groups=3153, threads=1;
	}
	private final Config c;
	private NcrnaDdlSavedData data;
	private Window[] windows;
	private NcrnaDdlProfile[][] profiles;
	private String[] names;
	private int[] selected, minima;
	private int k, width, score, nativeK, scoringWorkers;
	private boolean middle;
	private long scored, nativeTuples, canonicalTuples;
	private final AtomicInteger next=new AtomicInteger();
	private final AtomicReference<Throwable> failure=new AtomicReference<Throwable>();
	private static final int METHODS=3, TOTAL=METHODS+NcrnaDdlSavedData.CONTROLS;
	private static final int[] WIDTHS={256,512,1024,2048}, SCORES={NcrnaDdlRanks.SSU,NcrnaDdlRanks.ALL_KEYS};
	private static final String[] METHOD_NAMES={"NATIVE","SKETCH","HYBRID"};
}
